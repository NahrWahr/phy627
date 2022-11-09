#include "simulation.h"
#include "openacc.h"
#include "openacc_curand.h"
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>

void Simulation::init()
{
    t = 0;

    // Setup RGen distributions
    rgen.seed(p.seed);
    srand(p.seed * p.seed);
    std::uniform_real_distribution<real> RLx(0.0, p.Lx);
    std::uniform_real_distribution<real> RLy(0.0, p.Ly);
    std::uniform_real_distribution<real> RLz(0.0, p.Lz);

    std::normal_distribution<real> RVx(p.xVel, 0.5);
    std::normal_distribution<real> RVy(0.0, 0.5);
    std::normal_distribution<real> RVz(0.0, 0.5);

    RXi = std::normal_distribution<real>(0.0, 1.0);

    // Setup initial pos, vel and acc.
    for (int i = 0; i < p.numParticles; i++)
    {
        posx[i] = RLx(rgen);
        posy[i] = RLy(rgen);
        posz[i] = 0.0;

        velx[i] = RVx(rgen);
        vely[i] = RVy(rgen);
        velz[i] = 0.0;

        // vel[i].x = 0.0;
        // vel[i].y = 0.0;
        // vel[i].z = 0.0;

        accx[i] = 0.0;
        accy[i] = 0.0;
        accz[i] = 0.0;

        acc_oldx[i] = 0.0;
        acc_oldy[i] = 0.0;
        acc_oldz[i] = 0.0;

        pos_bufferx[i] = posx[i];
        pos_buffery[i] = posy[i];
        pos_bufferz[i] = posz[i];

        vel_bufferx[i] = velx[i];
        vel_buffery[i] = vely[i];
        vel_bufferz[i] = velz[i];
    }

// Update the data on the GPU:
#pragma acc update device(posx[:p.numParticles])
#pragma acc update device(posy[:p.numParticles])
#pragma acc update device(posz[:p.numParticles])

#pragma acc update device(velx[:p.numParticles])
#pragma acc update device(vely[:p.numParticles])
#pragma acc update device(velz[:p.numParticles])

#pragma acc update device(accx[:p.numParticles])
#pragma acc update device(accy[:p.numParticles])
#pragma acc update device(accz[:p.numParticles])

#pragma acc update device(acc_oldx[:p.numParticles])
#pragma acc update device(acc_oldy[:p.numParticles])
#pragma acc update device(acc_oldz[:p.numParticles])

#pragma acc update device(pos_bufferx[:p.numParticles * p.bufferSize])
#pragma acc update device(pos_buffery[:p.numParticles * p.bufferSize])
#pragma acc update device(pos_bufferz[:p.numParticles * p.bufferSize])

#pragma acc update device(vel_bufferx[:p.numParticles])
#pragma acc update device(vel_buffery[:p.numParticles])
#pragma acc update device(vel_bufferz[:p.numParticles])

    // Setup output streams
    osPos.open("outPos.csv");
    osVel.open("outVel.csv");
    int i = 0;
    osPos << std::setprecision(5) << std::fixed << "(" << posx[i] << " " << posy[i] << " " << posz[i] << ")";
    osVel << std::setprecision(5) << std::fixed << "(" << velx[i] << " " << vely[i] << " " << velz[i] << ")";
    for (i = 1; i < p.numParticles; i++)
    {
        osPos << ","
              << "(" << posx[i] << " " << posy[i] << " " << posz[i] << ")";
        osVel << ","
              << "(" << velx[i] << " " << vely[i] << " " << velz[i] << ")";
    }
}
/**
 * @brief The simulation step.
 * The basic steps are:
 * 1) Update r
 * 2) Update v_tilda
 * 3) Calculate F
 * 4) Update v
 * 
 */
void Simulation::step()
{
    // Define local variable, increment time.
    int n = p.numParticles;
    t++;
#pragma acc update device(t)

    // Setup random numbers
#if GPUACC
    int ng = 1024;
    long long *seed;
    seed = (long long *)malloc(ng * sizeof(long long));
    for (int i = 0; i < ng; ++i)
        seed[i] = rand();
#else
    int ng = 12;
    long long *seed;
    seed = (long long *)malloc(ng * sizeof(long long));
    for (int i = 0; i < ng; ++i)
        seed[i] = rand();
#endif

// Update r and calculate v tilde
#pragma acc parallel loop present(posx[:n], posy[:n], posz[:n], velx[:n], vely[:n], velz[:n], accx[:n], accy[:n], accz[:n], vel_tildex[:n], vel_tildey[:n], vel_tildez[:n], acc_oldx[:n], acc_oldy[:n], acc_oldz[:n], p)
    for (int i = 0; i < n; i++)
    {
        posx[i] += p.dt * (velx[i] + 0.5 * p.dt * accx[i]);
        posy[i] += p.dt * (vely[i] + 0.5 * p.dt * accy[i]);
        posz[i] += p.dt * (velz[i] + 0.5 * p.dt * accz[i]);

        vel_tildex[i] = velx[i] + p.lambda * p.dt * accx[i];
        vel_tildey[i] = vely[i] + p.lambda * p.dt * accy[i];
        vel_tildez[i] = velz[i] + p.lambda * p.dt * accz[i];

        acc_oldx[i] = accx[i];
        acc_oldy[i] = accy[i];
        acc_oldz[i] = accz[i];

        accx[i] = 0.0;
        accy[i] = 0.0;
        accz[i] = 0.0;
        
        
    }

// Calculate Force with v tilde on GPU
#if GPUACC
#pragma acc parallel loop gang vector present(posx[:n], posy[:n], posz[:n], vel_tildex[:n], vel_tildey[:n], vel_tildez[:n], accx[:n], accy[:n], accz[:n], p) copyin(seed[:ng])
    for (int i = 0; i < n; i++)
    {
        int gangidx = __pgi_gangidx();
        curandState_t state;
        curand_init(seed[gangidx] + i + 1, 0ULL, 0ULL, &state);
#pragma acc loop seq
        for (int j = i+1; j < n; j++)
        {
            interactionForceGPU(i, j, state);
        }
        
        // Add external force
        accx[i] += p.Fe;

        // Check wall colissions:
        if (!p.isPeriodic && posx[i] < p.r_c / 3)
        {
            real wposx, wposy, wposz;
            wposx = 0.0 - 2 * p.r_c / 3.0;
            wposy = posy[i];
            wposz = posz[i];
            interactionForceWallGPU(i, wposx, wposy, wposz, state);
        }
        if (posy[i] < p.r_c / 3)
        {
            real wposx, wposy, wposz;
            wposx = posx[i];
            wposy = 0.0 - 2 * p.r_c / 3.0;
            wposz = posz[i];
            interactionForceWallGPU(i, wposx, wposy, wposz, state);
        }
        if (!p.isPeriodic && posx[i] > p.Lx - p.r_c / 3)
        {
            real wposx, wposy, wposz;
            wposx = p.Lx + 2 * p.r_c / 3.0;
            wposy = posy[i];
            wposz = posz[i];
            interactionForceWallGPU(i, wposx, wposy, wposz, state);
        }
        if (posy[i] > p.Ly - p.r_c / 3)
        {
            real wposx, wposy, wposz;
            wposx = posx[i];
            wposy = p.Ly + 2 * p.r_c / 3.0;
            wposz = posz[i];
            interactionForceWallGPU(i, wposx, wposy, wposz, state);
        }
    }

// Calculate Force with v tilde on CPU
#else
#pragma acc parallel loop gang vector num_gangs(ng) present(posx[:n], posy[:n], posz[:n], vel_tildex[:n], vel_tildey[:n], vel_tildez[:n], accx[:n], accy[:n], accz[:n], p) copyin(seed[:ng])
    for (int i = 0; i < n; i++)
    {
        int gangidx = __pgi_gangidx();
        std::default_random_engine gen(seed[gangidx] + i + 1);
        std::normal_distribution<real> dist(0.0, 1.0);
#pragma acc loop seq
        for (int j = i+1; j < n; j++)
        {
            interactionForceCPU(i, j, gen, dist);
        }

        // Add external force
        accx[i] += p.Fe;
        
        // Check wall colissions:
        if (!p.isPeriodic && posx[i] < p.r_c / 3.0)
        {
            real wposx, wposy, wposz;
            wposx = 0.0 - 2 * p.r_c / 3.0;
            wposy = posy[i];
            wposz = posz[i];
            interactionForceWallCPU(i, wposx, wposy, wposz, gen, dist);
        }
        if (posy[i] < p.r_c / 3)
        {
            real wposx, wposy, wposz;
            wposx = posx[i];
            wposy = 0.0 - 2 * p.r_c / 3.0;
            wposz = posz[i];
            interactionForceWallCPU(i, wposx, wposy, wposz, gen, dist);
        }
        if (!p.isPeriodic && posx[i] > p.Lx - p.r_c / 3.0)
        {
            real wposx, wposy, wposz;
            wposx = p.Lx + 2 * p.r_c / 3.0;
            wposy = posy[i];
            wposz = posz[i];
            interactionForceWallCPU(i, wposx, wposy, wposz, gen, dist);
        }
        if (posy[i] > p.Ly - p.r_c / 3.0)
        {
            real wposx, wposy, wposz;
            wposx = posx[i];
            wposy = p.Ly + 2 * p.r_c / 3.0;
            wposz = posz[i];
            interactionForceWallCPU(i, wposx, wposy, wposz, gen, dist);
        }
        
    }
#endif

// Update v and reset F
#pragma acc parallel loop present(posx[:n], posy[:n], posz[:n], acc_oldx[:n], acc_oldy[:n], acc_oldz[:n], accx[:n], accy[:n], accz[:n], p)
    for (int i = 0; i < n; i++)
    {
        velx[i] = velx[i] + 0.5 * p.dt * (accx[i] + acc_oldx[i]);
        vely[i] = vely[i] + 0.5 * p.dt * (accy[i] + acc_oldy[i]);
        velz[i] = velz[i] + 0.5 * p.dt * (accz[i] + acc_oldz[i]);
        
        if (p.isSlipBoundary)
        {
            if (posx[i] < 0.0)
            {
                if (p.isPeriodic)
                    posx[i] += p.Lx;
                else
                {
                    posx[i] = -posx[i];
                    velx[i] = -velx[i];
                }
            }
            if (posy[i] < 0.0)
            {
                posy[i] = -posy[i];
                vely[i] = -vely[i];
            }
            if (posx[i] > p.Lx)
            {
                if (p.isPeriodic)
                    posx[i] -= p.Lx;
                else
                {
                    posx[i] = 2 * p.Lx - posx[i];
                    velx[i] = -velx[i];
                }
            }
            if (posy[i] > p.Ly)
            {
                posy[i] = 2 * p.Ly - posy[i];
                vely[i] = -vely[i];
            }
        }
        else
        {
            if (posx[i] < 0.0)
            {
                if (p.isPeriodic)
                    posx[i] += p.Lx;
                else
                {
                    posx[i] = -posx[i];
                    velx[i] = -velx[i];
                    vely[i] = -vely[i];
                    velz[i] = -velz[i];
                }
            }
            if (posy[i] < 0.0)
            {
                posy[i] = -posy[i];
                velx[i] = -velx[i];
                vely[i] = -vely[i];
                velz[i] = -velz[i];
            }
            if (posx[i] > p.Lx)
            {
                if (p.isPeriodic)
                    posx[i] -= p.Lx;
                else
                {
                    posx[i] = 2 * p.Lx - posx[i];
                    velx[i] = -velx[i];
                    vely[i] = -vely[i];
                    velz[i] = -velz[i];
                }
            }
            if (posy[i] > p.Ly)
            {
                posy[i] = 2 * p.Ly - posy[i];
                velx[i] = -velx[i];
                vely[i] = -vely[i];
                velz[i] = -velz[i];
            }
        }
    }

    // Copy data to CPU for asyncrohous print.

    if (p.outputEnabled)
    {
        int shift = t%p.bufferSize;
        int max_size = p.bufferSize - 1;
#pragma acc parallel loop present(pos_bufferx[:n], pos_buffery[:n], pos_bufferz[:n], vel_bufferx[:n], vel_buffery[:n], vel_bufferz[:n], posx[:n], posy[:n], posz[:n], velx[:n], vely[:n], velz[:n], p)
        for (int i = 0; i < n; i++)
        {
            int x = i + shift*p.numParticles;
            pos_bufferx[x] = posx[i];
            pos_buffery[x] = posy[i];
            pos_bufferz[x] = posz[i];
            
            vel_bufferx[x] = velx[i];
            vel_buffery[x] = vely[i];
            vel_bufferz[x] = velz[i];
        }
        if (shift == max_size)
        {
            #pragma acc update host(pos_bufferx[:p.numParticles*p.bufferSize], pos_buffery[:p.numParticles*p.bufferSize], pos_bufferz[:p.numParticles*p.bufferSize], vel_bufferx[:p.numParticles*p.bufferSize], vel_buffery[:p.numParticles*p.bufferSize], vel_bufferz[:p.numParticles*p.bufferSize])
            output();
        }
    }
#pragma acc wait
}

#if GPUACC
#pragma acc routine seq nohost
inline void Simulation::interactionForceGPU(int i, int j, curandState_t state)
{
    real rx, ry, rz, dist_sqr;
    rx = posx[i] - posx[j];
    ry = posy[i] - posy[j];
    rz = posz[i] - posz[j];

    dist_sqr = rx * rx + ry * ry + rz * rz;

    if (dist_sqr > p.r_c)
    {
        if (!p.isPeriodic)
            return;
        if (rx > 0)
            rx = rx - p.Lx;
        else
            rx = p.Lx - rx;
        dist_sqr = rx * rx + ry * ry + rz * rz;
        if (dist_sqr > p.r_c)
            return;
    }

    real dist, rnx, rny, rnz;
    real omegaC, omegaD, omegaR;
    real FC, FD, FR, F;
    real vx, vy, vz;

    vx = vel_tildex[i] - vel_tildex[j];
    vy = vel_tildey[i] - vel_tildey[j];
    vz = vel_tildez[i] - vel_tildez[j];

    dist = sqrt(dist_sqr);
    rnx = rx / dist;
    rny = ry / dist;
    rnz = rz / dist;

    omegaC = 1 - dist / p.r_c;
    omegaD = sqrt(omegaC);
    omegaR = sqrt(omegaD);

    real xi_temp;
    xi_temp = curand_normal_double(&state);

    FC = p.a * omegaC;
    FD = -p.gamma * omegaD * (vx * rnx + vy * rny + vz * rnz);
    FR = -p.sigma * omegaR * xi_temp * p.invSqrtDt;

    F = FC + FD + FR;

#pragma acc atomic update
    accx[i] += F * rnx;
#pragma acc atomic update
    accy[i] += F * rny;
#pragma acc atomic update
    accz[i] += F * rnz;
#pragma acc atomic update
    accx[j] -= F * rnx;
#pragma acc atomic update
    accy[j] -= F * rny;
#pragma acc atomic update
    accz[j] -= F * rnz;
}
#pragma acc routine seq nohost
void Simulation::interactionForceWallGPU(int i, real wposx, real wposy, real wposz, curandState_t state)
{
    real rx, ry, rz, dist_sqr;
    rx = posx[i] - wposx;
    ry = posy[i] - wposy;
    rz = posz[i] - wposz;

    dist_sqr = rx * rx + ry * ry + rz * rz;

    real dist, rnx, rny, rnz;
    real omegaC, omegaD, omegaR;
    real FC, FD, FR, F;
    real vx, vy, vz;

    vx = velx[i];
    vy = vely[i];
    vz = velz[i];

    dist = sqrt(dist_sqr);
    rnx = rx / dist;
    rny = ry / dist;
    rnz = rz / dist;

    omegaC = 1 - dist / p.r_c;
    omegaD = sqrt(omegaC);
    omegaR = sqrt(omegaD);

    real xi_temp;
    xi_temp = curand_normal_double(&state);

    FC = p.a * omegaC;
    FD = -p.gamma * omegaD * (vx * rnx + vy * rny + vz * rnz);
    FR = -p.sigma * omegaR * xi_temp * p.invSqrtDt;

    F = FC + FD + FR;

#pragma acc atomic update
    accx[i] += F * rnx;
#pragma acc atomic update
    accy[i] += F * rny;
#pragma acc atomic update
    accz[i] += F * rnz;
}

#else

void Simulation::interactionForceCPU(int i, int j,
                                     std::default_random_engine gen, std::normal_distribution<real> normal_dist)
{
    real rx, ry, rz, dist_sqr;
    rx = posx[i] - posx[j];
    ry = posy[i] - posy[j];
    rz = posz[i] - posz[j];

    dist_sqr = rx * rx + ry * ry + rz * rz;

    if (dist_sqr > p.r_c)
    {
        if (!p.isPeriodic)
            return;
        if (rx > 0)
            rx = rx - p.Lx;
        else
            rx = p.Lx - rx;
        dist_sqr = rx * rx + ry * ry + rz * rz;
        if (dist_sqr > p.r_c)
            return;
    }

    real dist, rnx, rny, rnz;
    real omegaC, omegaD, omegaR;
    real FC, FD, FR, F;
    real vx, vy, vz;

    vx = vel_tildex[i] - vel_tildex[j];
    vy = vel_tildey[i] - vel_tildey[j];
    vz = vel_tildez[i] - vel_tildez[j];

    dist = sqrt(dist_sqr);
    rnx = rx / dist;
    rny = ry / dist;
    rnz = rz / dist;

    omegaC = 1 - dist / p.r_c;
    omegaD = sqrt(omegaC);
    omegaR = sqrt(omegaD);

    real xi_temp;
    xi_temp = normal_dist(gen);

    FC = p.a * omegaC;
    FD = -p.gamma * omegaD * (vx * rnx + vy * rny + vz * rnz);
    FR = -p.sigma * omegaR * xi_temp * p.invSqrtDt;

    F = FC + FD + FR;

#pragma acc atomic update
    accx[i] += F * rnx;
#pragma acc atomic update
    accy[i] += F * rny;
#pragma acc atomic update
    accz[i] += F * rnz;
#pragma acc atomic update
    accx[j] -= F * rnx;
#pragma acc atomic update
    accy[j] -= F * rny;
#pragma acc atomic update
    accz[j] -= F * rnz;
}

void Simulation::interactionForceWallCPU(int i, real wposx, real wposy, real wposz, std::default_random_engine gen, std::normal_distribution<real> normal_dist)
{
    real rx, ry, rz, dist_sqr;
    rx = posx[i] - wposx;
    ry = posy[i] - wposy;
    rz = posz[i] - wposz;

    dist_sqr = rx * rx + ry * ry + rz * rz;

    real dist, rnx, rny, rnz;
    real omegaC, omegaD, omegaR;
    real FC, FD, FR, F;
    real vx, vy, vz;

    vx = velx[i];
    vy = vely[i];
    vz = velz[i];

    dist = sqrt(dist_sqr);
    rnx = rx / dist;
    rny = ry / dist;
    rnz = rz / dist;

    omegaC = 1 - dist / p.r_c;
    omegaD = sqrt(omegaC);
    omegaR = sqrt(omegaD);

    real xi_temp;
    xi_temp = normal_dist(gen);

    FC = p.a * omegaC;
    FD = -p.gamma * omegaD * (vx * rnx + vy * rny + vz * rnz);
    FR = -p.sigma * omegaR * xi_temp * p.invSqrtDt;

    F = FC + FD + FR;

#pragma acc atomic update
    accx[i] += F * rnx;
#pragma acc atomic update
    accy[i] += F * rny;
#pragma acc atomic update
    accz[i] += F * rnz;
}

#endif

void Simulation::output()
{
    int shift = p.bufferSize;
    for (int i = 0; i < p.bufferSize*p.numParticles; i++)
    {
        if (i%p.numParticles==0){

            osPos << std::endl;
            osVel << std::endl;
            osPos << "(" << pos_bufferx[i] << " " << pos_buffery[i] << " " << pos_bufferz[i] << ")";
            osVel << "(" << vel_bufferx[i] << " " << vel_buffery[i] << " " << vel_bufferz[i] << ")";
        }
        else {
            osPos << ","
                  << "(" << pos_bufferx[i] << " " << pos_buffery[i] << " " << pos_bufferz[i] << ")";
            osVel << ","
                  << "(" << vel_bufferx[i] << " " << vel_buffery[i] << " " << vel_bufferz[i] << ")";
        }
    }
}

Simulation::Simulation()
{
}

/**
 * @brief Construct a new Simulation:: Simulation object
 * 
 * @param params 
 */
Simulation::Simulation(Parameters params)
{
    p = params;

    // Setup memory
    posx = new real[p.numParticles];
    posy = new real[p.numParticles];
    posz = new real[p.numParticles];

    velx = new real[p.numParticles];
    vely = new real[p.numParticles];
    velz = new real[p.numParticles];

    vel_tildex = new real[p.numParticles];
    vel_tildey = new real[p.numParticles];
    vel_tildez = new real[p.numParticles];

    accx = new real[p.numParticles];
    accy = new real[p.numParticles];
    accz = new real[p.numParticles];

    acc_oldx = new real[p.numParticles];
    acc_oldy = new real[p.numParticles];
    acc_oldz = new real[p.numParticles];

    pos_bufferx = new real[p.numParticles * p.bufferSize];
    pos_buffery = new real[p.numParticles * p.bufferSize];
    pos_bufferz = new real[p.numParticles * p.bufferSize];

    vel_bufferx = new real[p.numParticles * p.bufferSize];
    vel_buffery = new real[p.numParticles * p.bufferSize];
    vel_bufferz = new real[p.numParticles * p.bufferSize];

#pragma acc enter data create(this)
#pragma acc update device(this)
#pragma acc enter data create(p)

#pragma acc enter data create(posx[:p.numParticles])
#pragma acc enter data create(posy[:p.numParticles])
#pragma acc enter data create(posz[:p.numParticles])

#pragma acc enter data create(velx[:p.numParticles])
#pragma acc enter data create(vely[:p.numParticles])
#pragma acc enter data create(velz[:p.numParticles])

#pragma acc enter data create(vel_tildex[:p.numParticles])
#pragma acc enter data create(vel_tildey[:p.numParticles])
#pragma acc enter data create(vel_tildez[:p.numParticles])

#pragma acc enter data create(accx[:p.numParticles])
#pragma acc enter data create(accy[:p.numParticles])
#pragma acc enter data create(accz[:p.numParticles])

#pragma acc enter data create(acc_oldx[:p.numParticles])
#pragma acc enter data create(acc_oldy[:p.numParticles])
#pragma acc enter data create(acc_oldz[:p.numParticles])

#pragma acc enter data create(pos_bufferx[:p.numParticles * p.bufferSize])
#pragma acc enter data create(pos_buffery[:p.numParticles * p.bufferSize])
#pragma acc enter data create(pos_bufferz[:p.numParticles * p.bufferSize])

#pragma acc enter data create(vel_bufferx[:p.numParticles * p.bufferSize])
#pragma acc enter data create(vel_buffery[:p.numParticles * p.bufferSize])
#pragma acc enter data create(vel_bufferz[:p.numParticles * p.bufferSize])
}

Simulation::~Simulation()
{
    // Free memory
    delete (posx);
    delete (posy);
    delete (posz);

    delete (velx);
    delete (vely);
    delete (velz);

    delete (vel_tildex);
    delete (vel_tildey);
    delete (vel_tildez);

    delete (accx);
    delete (accy);
    delete (accz);

    delete (acc_oldx);
    delete (acc_oldy);
    delete (acc_oldz);

    delete (pos_bufferx);
    delete (pos_buffery);
    delete (pos_bufferz);

    delete (vel_bufferx);
    delete (vel_buffery);
    delete (vel_bufferz);

#pragma acc exit data delete (posx[:p.numParticles])
#pragma acc exit data delete (posy[:p.numParticles])
#pragma acc exit data delete (posz[:p.numParticles])

#pragma acc exit data delete (velx[:p.numParticles])
#pragma acc exit data delete (vely[:p.numParticles])
#pragma acc exit data delete (velz[:p.numParticles])

#pragma acc exit data delete (vel_tildex[:p.numParticles])
#pragma acc exit data delete (vel_tildey[:p.numParticles])
#pragma acc exit data delete (vel_tildez[:p.numParticles])

#pragma acc exit data delete (accx[:p.numParticles])
#pragma acc exit data delete (accy[:p.numParticles])
#pragma acc exit data delete (accz[:p.numParticles])

#pragma acc exit data delete (acc_oldx[:p.numParticles])
#pragma acc exit data delete (acc_oldy[:p.numParticles])
#pragma acc exit data delete (acc_oldz[:p.numParticles])

#pragma acc exit data delete (pos_bufferx[:p.numParticles * p.bufferSize])
#pragma acc exit data delete (pos_buffery[:p.numParticles * p.bufferSize])
#pragma acc exit data delete (pos_bufferz[:p.numParticles * p.bufferSize])

#pragma acc exit data delete (vel_bufferx[:p.numParticles * p.bufferSize])
#pragma acc exit data delete (vel_buffery[:p.numParticles * p.bufferSize])
#pragma acc exit data delete (vel_bufferz[:p.numParticles * p.bufferSize])

#pragma acc exit data delete (p)
#pragma acc exit data delete (this)
}
