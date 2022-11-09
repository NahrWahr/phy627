#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <iostream>
#include <fstream>
#include <random>
#include "openacc_curand.h"

typedef double real;
typedef struct
{
    real x;
    real y;
    real z;
} real3;

class Parameters
{
public:
    int numSoParticles;
    int numSiParticles;
    int numWParticles;
    int numCParticles;
    int numParticles;

    int seed;
    int nt;
    bool isPeriodic;
    bool isSlipBoundary;
    bool gpuAcclerated;
    bool outputEnabled;
    int bufferSize;

    double Lx;
    double Ly;
    double Lz;
    double Fe;
    double xVel;
    double r_c;
    double a;
    double sigma;
    double gamma;
    double lambda;
    double dt;
    double invSqrtDt;

    void outputParams();
};

class Simulation
{
public:
    //IO Streams
    std::ofstream osPos;
    std::ofstream osVel;

    //Data
    Parameters p;
    real *posx;
    real *posy;
    real *posz;
    real *velx;
    real *vely;
    real *velz;
    real *vel_tildex;
    real *vel_tildey;
    real *vel_tildez;
    real *accx;
    real *accy;
    real *accz;
    real *acc_oldx;
    real *acc_oldy;
    real *acc_oldz;
    real *pos_bufferx;
    real *pos_buffery;
    real *pos_bufferz;
    real *vel_bufferx;
    real *vel_buffery;
    real *vel_bufferz;
    int t;

    //RGen
    std::default_random_engine rgen;
    std::normal_distribution<real> RXi;

    //Main Methods
    void init();
    void step();

    //Auxilary Methods
    void interactionForceGPU(int i, int j, curandState_t state);
    void interactionForceCPU(int i, int j, std::default_random_engine gen, std::normal_distribution<real> normal_dist);
    void interactionForceWallGPU(int i, real wposx, real wposy, real wposz, curandState_t state);
    void interactionForceWallCPU(int i, real wposx, real wposy, real wposz, std::default_random_engine gen, std::normal_distribution<real> normal_dist);

    //Helper Methods
    void output();
    //Data Methods

    //Constructors
    Simulation();
    Simulation(Parameters params);
    ~Simulation();
};

#endif // __SIMULATION_H__
