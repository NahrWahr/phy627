#include "../include/simulation.h"
#include "../include/base.h"
#include <nvtx3/nvToolsExt.h>

/**
 * @brief Construct a new Simulation:: Simulation object
 * 
 * No params mentioned. Init everything to zero.
 */
Simulation::Simulation()
{
    Parameters params;
    params.numSoParticles = 0;
    params.numSiParticles = 0;
    params.numWParticles = 0;
    params.numWParticles = 0;
    params.numParticles = 0;
    params.seed = 2021;
    this->p = params;

    // Init random variables
    rgen.seed(p.seed);
    std::uniform_real_distribution<double> RLx = std::uniform_real_distribution<double>(0, p.Lx);
    std::uniform_real_distribution<double> RLy = std::uniform_real_distribution<double>(0, p.Ly);
    std::uniform_real_distribution<double> RLz = std::uniform_real_distribution<double>(0, p.Lz);
    std::normal_distribution<double> xi = std::normal_distribution<double>(0.0, 1.0);

    initWParticles();
    initSoParticles();
    initCParticles();
    initSiParticles();

    time = 0;

    position.open("positions.csv");
}
/**
 * @brief Construct a new Simulation:: Simulation object
 * 
 * @param params : use params to initialize a simulation.
 */
Simulation::Simulation(Parameters params)
{
    this->p = params;

    // Init random variables
    rgen.seed(p.seed);
    RLx = *new std::uniform_real_distribution<double>(0, p.Lx);
    RLy = *new std::uniform_real_distribution<double>(0, p.Ly);
    RLz = *new std::uniform_real_distribution<double>(0, p.Lz);
    xi = *new std::normal_distribution<double>(0.0, 1.0);

    initWParticles();
    initSoParticles();
    initCParticles();
    initSiParticles();

    time = 0;

    position.open("positions.csv");
    velocity.open("velocities.csv");
}

/**
 * @brief Return set of parameters
 * 
 * @return Parameters 
 */
Parameters Simulation::getParams()
{
    return p;
}

/**
 * @brief Initialize the So type particles and add them to the particles vector.
 * 
 */
void Simulation::initSoParticles()
{
    for (int i = 0; i < p.numSoParticles; i++)
    {
        Particle *x = new Particle("So");
        // x->pos = Vec3(RLx(rgen), RLy(rgen), RLz(rgen));
        x->pos = Vec3(RLx(rgen), RLy(rgen), 0);
        x->vel = Vec3(xi(rgen) + 1.0, xi(rgen), 0);
        particles.push_back(*x);
    }
}

/**
 * @brief Initialize the Si type particles and add them to the particles vector.
 * 
 */
void Simulation::initSiParticles()
{
    for (int i = 0; i < p.numSiParticles; i++)
    {
        Particle *x = new Particle("Si");
        particles.push_back(*x);
    }
}

/**
 * @brief Initialize the W type particles and add them to the particles vector.
 * 
 */
void Simulation::initWParticles()
{
    for (int i = 0; i < p.numWParticles; i++)
    {
        Particle *x = new Particle("W");

        particles.push_back(*x);
    }
}

/**
 * @brief Initialize the C type particles and add them to the particles vector.
 * 
 */
void Simulation::initCParticles()
{
    for (int i = 0; i < p.numCParticles; i++)
    {
        Particle *x = new Particle("C");
        particles.push_back(*x);
    }
}

/**
 * @brief Make the simulation take a single step (dt) in time.
 * Involves calculating forces, doing first half of the interactions, 
 * calculating more forces, doing second half of the interaction.
 * 
 */
void Simulation::step()
{
    {
        // Calculating forces.
        nvtxRangePushA("CalculateForce1");
        for (int i = 0; i < p.numParticles; i++)
        {
            calculateForces(i);
        }
        nvtxRangePop();

        // Do first half of interactions
        nvtxRangePushA("firstHalfInteractions");
        for (int i = 0; i < p.numParticles; i++)
        {
            firstHalfInteractions(i);
        }
        nvtxRangePop();

        // Calculate forces again
        nvtxRangePushA("CalculateForce2");
        for (int i = 0; i < p.numParticles; i++)
        {
            calculateForces(i);
        }
        nvtxRangePop();

        // Calculate second half of interactions
        nvtxRangePushA("secondHalfInteractions");
        for (int i = 0; i < p.numParticles; i++)
        {
            secondHalfInteractions(i);
        }
        nvtxRangePop();
    }
    time++;
}

/**
 * @brief Calculate all the forces on a particle i.
 * 
 * @param i Which ith particle to calculate the forces on.
 */
void Simulation::calculateForces(int i)
{
    double x1, x2, x3;
    x1 = particles[i].pos.x;
    x2 = particles[i].pos.y;
    x3 = particles[i].pos.z;
    {
        if (x1 < p.r_c)
        {
            if (!p.periodicBoundary)
            {
                applyWallForce(i, particles[i].pos + Vec3(-particles[i].pos.x, 0.0, 0.0));
            }
        }
        if (x2 < p.r_c)
        {
            applyWallForce(i, particles[i].pos + Vec3(0.0, -particles[i].pos.y, 0.0));
        }
        // if (x3 < p.r_c)
        // {
        //     applyWallForce(i, particles[i].pos + Vec3(0.0, 0.0, -particles[i].pos.z));
        // }
        if (x1 > p.Lx - p.r_c)
        {
            if (!p.periodicBoundary)
            {
                applyWallForce(i, particles[i].pos + Vec3(p.Lx - particles[i].pos.x, 0.0, 0.0));
            }
        }
        if (x2 > p.Ly - p.r_c)
        {
            applyWallForce(i, particles[i].pos + Vec3(0.0, p.Ly - particles[i].pos.y, 0.0));
        }
        // if (x3 > p.Lz - p.r_c)
        // {
        //     applyWallForce(i, particles[i].pos + Vec3(0.0, 0.0, p.Lz - particles[i].pos.z));
        // }
    }

    // Interaction with other particles
    // ! Note: I have counter is running from i+1 to n_particles. Make sure to observe while setting up parallel code.
    for (int j = i + 1; j < p.numParticles; j++)
    {
        Vec3 r_i, r_j, r_ij;
        Vec3 v_i, v_j, v_ij;
        bool flag = false;

        r_i = particles[i].pos;
        r_j = particles[j].pos;
        r_ij = r_i - r_j;

        double dist = r_ij.norm();
        if (dist < p.r_c)
            flag = true;
        else if (p.periodicBoundary)
        {
            Vec3 boundary = Vec3(p.Lx, 0, 0);
            // cout<<(boundary);
            // double periodic_distance;
            if (r_i.x > p.Lx - 1)
            {
                r_ij = r_i - r_j - boundary;
                // cout<<(r_ij);
                dist = r_ij.norm();
                if (dist < p.r_c)
                    // cout<< dist;
                    flag = true;
            }
            else if (r_j.x > p.Lx - 1)
            {
                r_ij = r_i - r_j + boundary;
                dist = r_ij.norm();
                if (dist < p.r_c)
            //         cout<<dist;
                    flag = true;
            }
        }
        // Check if within distance
        if (!flag)
            continue;

        // Calculate force

        static double temp, omegaD, omegaR, F_C, F_D, F_R;
        static Vec3 DeltaF;

        v_i = particles[i].vel;
        v_j = particles[j].vel;
        v_ij = v_i - v_j;

        temp = (1 - r_ij.norm() / p.r_c);
        omegaD = sqrt(temp);
        omegaR = sqrt(omegaD);

        F_C = p.a * temp;
        F_D = -p.gamma * omegaD * (v_ij.dot(r_ij.normalize()));
        F_R = -p.sigma * omegaR * (p.invSqrtDt) * (double)xi(rgen);
        // std::cout<<xi(rgen);
        DeltaF = r_ij * (F_C + F_D + F_R);

        particles[i].acc += DeltaF;
        particles[j].acc -= DeltaF;
    }
}

void Simulation::applyWallForce(int i, Vec3 wallPos)
{
    static Vec3 r_ij;
    static Vec3 v_ij;

    r_ij = particles[i].pos - wallPos;
    v_ij = particles[i].vel;

    static double temp, omegaD, omegaR, F_C, F_D, F_R;
    static Vec3 DeltaF;

    temp = (1 - r_ij.norm() / p.r_c);
    omegaD = sqrt(temp);
    omegaR = sqrt(omegaD);

    F_C = p.a * temp;
    F_D = -p.gamma * omegaD * (v_ij.dot(r_ij.normalize()));
    F_R = -p.sigma * omegaR * (p.invSqrtDt) * (double)xi(rgen);
    DeltaF = r_ij * (F_C + F_D + F_R);
    // std::cout <<F_C<<F_D<<F_R;

    particles[i].acc += DeltaF;
    // std::cout << DeltaF;
}

void Simulation::firstHalfInteractions(int i)
{
    particles[i].pos = particles[i].pos + particles[i].vel * p.dt + particles[i].acc * p.dt * p.dt * 0.5;
    specularReflection(i);
    particles[i].velOld = particles[i].vel * 1.0;
    particles[i].accOld = particles[i].acc * 1.0;
    particles[i].vel = particles[i].vel + particles[i].acc * p.lambda * p.dt;
    particles[i].acc = Vec3(0, 0, 0);
}

void Simulation::secondHalfInteractions(int i)
{
    particles[i].vel = particles[i].velOld + (particles[i].acc + particles[i].accOld) * p.dt * 0.5;
    particles[i].acc = Vec3(0, 0, 0);
}

void Simulation::specularReflection(int i)
{
    double x1, x2, x3;
    x1 = particles[i].pos.x;
    x2 = particles[i].pos.y;
    x3 = particles[i].pos.z;
    if (x1 < 0)
    {
        if (p.periodicBoundary)
        {
            particles[i].pos.x += p.Lx;
        }
        else
        {
            particles[i].pos.x = -particles[i].pos.x;
            particles[i].vel.x = -particles[i].vel.x;
        }
    }
    if (x2 < 0)
    {
        particles[i].pos.y = -particles[i].pos.y;
        particles[i].vel.y = -particles[i].vel.y;
    }
    if (x3 < 0)
    {
        particles[i].pos.z = -particles[i].pos.z;
        particles[i].vel.z = -particles[i].vel.z;
    }
    if (x1 > p.Lx)
    {
        if (p.periodicBoundary)
        {
            particles[i].pos.x -= p.Lx;
        }
        else
        {
            particles[i].pos.x = 2 * p.Lx - particles[i].pos.x;
            particles[i].vel = Vec3(-particles[i].vel.x, 0,0);
        }
    }
    if (x2 > p.Ly)
    {
        particles[i].pos.y = 2 * p.Ly - particles[i].pos.y;
        particles[i].vel = Vec3(0,-particles[i].vel.y,0);
    }
    if (x3 > p.Ly)
    {
        particles[i].pos.z = 2 * p.Lz - particles[i].pos.z;
        particles[i].vel = Vec3(0,0,-particles[i].vel.z);
    }
}

/**
 * @brief Function to output the position data
 * 
 */
void Simulation::outputPositionData()
{
    // Output first line
    static bool flag = true;
    if (flag)
    {
        position << particles[0].type;
        for (int i = 1; i < p.numParticles; i++)
        {

            position << "," << particles[i].type;
            // position << "," << particles[i].type;
            // position << "," << particles[i].type;
        }
        flag = false;
    }

    // Output positions:
    position << std::endl;
    position << particles[0].pos;

    for (int i = 1; i < p.numParticles; i++)
    {
        position << "," << particles[i].pos;
    }
}

void Simulation::outputVelocityData() 
{
    // Output first line
    static bool flag = true;
    if (flag)
    {
        velocity << particles[0].type;
        for (int i = 1; i < p.numParticles; i++)
        {

            velocity << "," << particles[i].type;
            // position << "," << particles[i].type;
            // position << "," << particles[i].type;
        }
        flag = false;
    }

    // Output velocities:
    velocity << std::endl;
    velocity << particles[0].vel;

    for (int i = 1; i < p.numParticles; i++)
    {
        velocity << "," << particles[i].vel;
    }
}
