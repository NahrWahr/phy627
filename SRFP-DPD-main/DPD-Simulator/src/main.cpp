#include <iostream>
#include "../include/base.h"
#include "../include/simulation.h"
#include <nvtx3/nvToolsExt.h>

Parameters setParams();

int main(int argc, char *argv[])
{
	Parameters p = setParams();

    // Initialize the simulation
    nvtxRangePushA("init");
    Simulation sim = Simulation(p);
    nvtxRangePop();

    // Output the parameters for the first time
    nvtxRangePushA("outputParams");
    sim.p.outputParams();
    nvtxRangePop();

    // Output postiion data for the first time.
    nvtxRangePushA("OutputPositionData");
    sim.outputPositionData();
    sim.outputVelocityData();
    nvtxRangePop();

    // Start iterations.
    nvtxRangePushA("Steps");
    for (int i = 0; i < sim.p.steps; i++)
    {
        // Single step
        nvtxRangePushA("Step");
        sim.step();
        nvtxRangePop();

        // After taking step, output data.
        nvtxRangePushA("OutputPositionData");
        sim.outputPositionData();
        sim.outputVelocityData();
        nvtxRangePop();
    }
    nvtxRangePop();
    
}

/**
 * @brief Set the Params object with data. Modify at compile time.
 * 
 * @return Parameters 
 */
Parameters setParams(){
    Parameters p;
    p.periodicBoundary = true;
    p.seed = 2021;
    p.steps = 2000;
    p.numSoParticles = 1000;
    p.numWParticles = 0;
    p.numSiParticles = 0;
    p.numCParticles = 0;
    p.numParticles = p.numSiParticles + p.numCParticles + p.numSoParticles + p.numWParticles;
    
    p.Ly = p.Lz = 10.0;
    p.Lx = 40.0;
    p.r_c = 1.0;
    p.gamma = 18;
    p.a = 20;
    p.dt = 0.01;
    p.invSqrtDt = 1/sqrt(p.dt);
    p.sigma = sqrt((double)2*p.gamma*1); //310K * K_b = x j
	p.lambda = 0.65;
    return p;
}
