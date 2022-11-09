#include "../include/base.h"
#include <math.h>
#include <fstream>

Vec3::Vec3()
{
    x = 0;
    y = 0;
    z = 0;
}

Vec3::Vec3(vecType x1, vecType y1, vecType z1)
{
    x = x1;
    y = y1;
    z = z1;
}

Vec3::Vec3(const Vec3 &vec)
{
    x = vec.x;
    y = vec.y;
    z = vec.z;
}

Vec3 Vec3::operator-(const Vec3 &vec)
{
    return Vec3(x - vec.x, y - vec.y, z - vec.z);
}

Vec3 &Vec3::operator-=(const Vec3 &vec)
{
    this->x -= vec.x;
    this->y -= vec.y;
    this->z -= vec.z;
    return *this;
}

Vec3 Vec3::operator/(vecType value)
{
    return Vec3(x / value, y / value, z / value);
}

Vec3 &Vec3::operator/=(vecType value)
{
    this->x /= value;
    this->y /= value;
    this->z /= value;
    return *this;
}

Vec3 &Vec3::operator=(const Vec3 &vec)
{
    this->x = vec.x;
    this->y = vec.y;
    this->z = vec.z;
    return *this;
}

Vec3 &Vec3::operator*=(vecType value)
{
    this->x *= value;
    this->y *= value;
    this->z *= value;
    return *this;
}

Vec3 Vec3::operator*(vecType value)
{
    return Vec3(this->x * value, this->y * value, this->z * value);
}

Vec3 &Vec3::operator+=(const Vec3 &vec)
{
    this->x += vec.x;
    this->y += vec.y;
    this->z += vec.z;
    return *this;
}

Vec3 Vec3::operator+(const Vec3 &vec)
{
    return Vec3(this->x + vec.x, this->y + vec.y, this->z + vec.z);
}

vecType Vec3::dot(const Vec3 &vec)
{
    return this->x * vec.x + this->y * vec.y + this->z * vec.z;
}

Vec3 Vec3::cross(const Vec3 &vec)
{
    vecType x1 = this->y * vec.z - this->z * vec.y;
    vecType y1 = this->z * vec.x - this->x * vec.z;
    vecType z1 = this->x * vec.y - this->y * vec.x;
    return Vec3(x1, y1, z1);
}

vecType Vec3::norm()
{
    return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

Vec3 Vec3::normalize()
{
    vecType n = this->norm();
    return *this / n;
}

vecType Vec3::distance(const Vec3 &vec)
{
    return sqrt((this->x - vec.x) * (this->x - vec.x) +
                (this->y - vec.y) * (this->y - vec.y) +
                (this->z - vec.z) * (this->z - vec.z));
}

Particle::Particle(string type1) 
{
    this->type = type1;
}

void Parameters::outputParams() 
{
    std::ofstream output;
    output.open("params.csv");
    output << "Parameter," << "Value"<<std::endl;

    output << "Num So Particles," << numSoParticles <<std::endl;
    output << "Num Si Particles," << numSiParticles <<std::endl;
    output << "Num W Particles," << numWParticles <<std::endl;
    output << "Num C Particles," << numCParticles <<std::endl;
    output << "Num Particles," << numParticles <<std::endl;

    output << "Lx," << Lx << std::endl;
    output << "Ly," << Ly << std::endl;
    output << "Lz," << Lz << std::endl;
    output << "r_c," << r_c << std::endl;
    output << "a," << a << std::endl;
    output << "sigma," << sigma << std::endl;
    output << "gamma," << gamma << std::endl;
    output << "lambda," << lambda << std::endl;
    output << "dt," << dt;
    
    output.close(); 
    
}

ostream &operator<<(ostream &os, const Vec3 &vec)
{
    os<<"("<<vec.x<<" "<<vec.y<<" "<<vec.z<<")";
    return os;
}
