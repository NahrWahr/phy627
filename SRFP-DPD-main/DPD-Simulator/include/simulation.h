#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <fstream>
#include <random>
#include <vector>
#include "base.h"

class Simulation
{
public:
    // IO streams
    std::ofstream position;
    std::ofstream velocity;
    std::ofstream acceleration;

    // Data
    Parameters p;
    std::vector<Particle> particles;
    std::default_random_engine rgen;
    std::uniform_real_distribution<double> RLx;
    std::uniform_real_distribution<double> RLy;
    std::uniform_real_distribution<double> RLz;
    std::normal_distribution<double> xi;
    int time;
    // Constructors
    Simulation();
    Simulation(Parameters params);

    // Methods
    Parameters getParams();

    void initSoParticles();
    void initSiParticles();
    void initWParticles();
    void initCParticles();

    void step();
    void calculateForces(int i);
    void applyWallForce(int i, Vec3 wallPos);
    void applyForces(int i);
    void firstHalfInteractions(int i);
    void secondHalfInteractions(int i);
    void specularReflection(int i);

    void outputPositionData();
    void outputVelocityData();
};

#endif //__SIMULATION_H__
