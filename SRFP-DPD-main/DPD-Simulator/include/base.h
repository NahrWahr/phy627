#ifndef __BASE_H__
#define __BASE_H__

#include <string>
#include <iostream>
#include <math.h>

using namespace std;
typedef double vecType;

class Vec3
{
public:
    vecType x, y, z;

    Vec3(); // Default constructor;
    Vec3(vecType x1, vecType y1, vecType z1 = 0);
    Vec3(const Vec3 &vec); //copy constructor

    //Operator overloads
    Vec3 operator+(const Vec3 &vec);
    Vec3 &operator+=(const Vec3 &vec);
    Vec3 operator-(const Vec3 &vec);
    Vec3 &operator-=(const Vec3 &vec);
    Vec3 operator*(vecType value);
    Vec3 &operator*=(vecType value);
    Vec3 operator/(vecType value);
    Vec3 &operator/=(vecType value);
    Vec3 &operator=(const Vec3 &vec);
    friend ostream &operator<<(ostream &os, const Vec3 &vec);

    //Functions
    vecType dot(const Vec3 &vec);
    Vec3 cross(const Vec3 &vec);
    vecType norm();
    Vec3 normalize();
    vecType distance(const Vec3 &vec);
};

class Particle
{
public:
    Vec3 pos;
    Vec3 vel;
    Vec3 velOld;
    Vec3 acc;
    Vec3 accOld;
    string type;
    Particle(string type1);
};

class Parameters
{
public:
    int numSoParticles;
    int numSiParticles;
    int numWParticles;
    int numCParticles;
    int numParticles;
    uint32_t seed;
    int steps;
    bool periodicBoundary;

    double Lx;
    double Ly;
    double Lz;
    double r_c;
    double a;
    double sigma;
    double gamma;
    double lambda;
    double dt;
    double invSqrtDt;

    void outputParams();
    
};

#endif // __BASE_H__