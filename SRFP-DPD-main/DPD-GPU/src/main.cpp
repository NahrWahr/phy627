/**
 * @file main.cpp
 * @author Debaditya Bhattacharya (debbh922@gmail.com)
 * @brief Simulation driver. Parameters (in order):
 * @param numSoParticles Number of particles
 * @param Lx Box length 
 * @param Ly Box width
 * @param nt Number of steps
 * @version 0.1
 * @date 2021-07-17
 */

#include <iostream>
#include <chrono>
#include "simulation.h"
#include <openacc.h>
#include <fstream>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

Parameters setDefaultParams();

int main(int argc, char *argv[])
{
    Parameters p;
    p = setDefaultParams();
    p.outputEnabled = true;
    int x = 1;
    if(argc > 1){
        try{
            p.numSoParticles = atoi(argv[x++]);
            p.Lx             = atof(argv[x++]);
            p.Ly             = atof(argv[x++]);
            p.Lz             = atof(argv[x++]);
            p.Fe             = atof(argv[x++]);
            p.xVel           = atof(argv[x++]);
            p.nt             = atoi(argv[x++]);
            p.outputEnabled  = atoi(argv[x++]);
            p.isSlipBoundary = atoi(argv[x++]);
            p.isPeriodic     = atoi(argv[x++]);
            p.seed           = atoi(argv[x++]);
        }
        catch (int x) {
            std::cout<<"Invalid argument please check parameter " <<x;;
            exit(-1);
        }
    }
    p.numParticles   = p.numCParticles + p.numSiParticles + p.numSoParticles + p.numWParticles;
    std::ofstream time;
    time.open("time.txt");

    auto t1 = high_resolution_clock::now();
    Simulation sim(p);
    sim.init();
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    std::cout<<"Time taken for init: "<<ms_double.count()<<"ms \n";
    time<<ms_double.count()<<std::endl;

    t1 = high_resolution_clock::now();
    std::cout<<"Processing step 0";
    for (int t = 0; t < sim.p.nt; t++)
    {
        sim.step();
        if (t%100==0){
            std::cout<<"\rProcessing step "<<t;
            fflush(stdout);
        }
    }
    std::cout<<std::endl;
    t2 = high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout<<"Time taken for simulation of "<<p.numParticles <<" particles for "<<p.nt<<" steps: "<<ms_double.count()<<"ms \n";
    time<<ms_double.count()<<std::endl;
    return 0;
}

Parameters setDefaultParams()
{
    Parameters p;
    // Particle Count
    p.numSoParticles = 3200;
    p.numSiParticles = 0;
    p.numWParticles  = 0;
    p.numCParticles  = 0;
    p.numParticles   = p.numCParticles + p.numSiParticles + p.numSoParticles + p.numWParticles;

    // Simulation Paramters
    p.gpuAcclerated  = true;
    p.isSlipBoundary = false;
    p.isPeriodic     = true;
    p.seed           = 2020;
    p.nt             = 2000;
    p.bufferSize     = 10000;
    p.Lx             = 40;
    p.Ly             = 20;
    p.Lz             = 0;
    p.Fe             = 0.01;

    p.r_c            = 1;
    p.a              = 20;
    p.gamma          = 18;
    p.sigma          = sqrt((double)2 * p.gamma * 1);
    p.lambda         = 0.65;
    p.dt             = 0.01;
    p.invSqrtDt      = 1 / sqrt(p.dt);

    return p;
}
