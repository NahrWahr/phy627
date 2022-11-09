#include<bits/stdc++.h>
using namespace std;

class FluidSim
{
public:
  struct Particle
  {
    //vector<double> r = vector<double>(3);
    vector<double> r{3};
    vector<double> v{3};

    //vector<int> intpos = vector<int>(3);
    vector<int> intpos{3};
  };

  struct Accel
  {
    vector<double> a{3};
  };

private:
  Particle *P, *Pend;
  int np;
  
  Accel *A;
  double dt;
  double L;
  double cutoff_r;

  int npcf;
  double pcf_rmax;
  double *pcf_sum;
  double *pcf_sum2;
  int pcf_samples;

  double r2ir;
  int irmax;

public:
  // Some simulation statistics: 
  int nskip0;
  int nskip1;
  int ncalc;
  int nborder;
  double r_min;
		
private:
  static inline double potential_1r2(double r2)
  {
    double r6=r2*r2*r2;
    return(4.0*r6*(r6-1.0));
  }

  static inline double potential(double r)
  {  return(potential_1r2(1.0/(r*r)));  }
		
  static inline double potentialD_1r2(double r2)
  {
    double r6=r2*r2*r2;
    return(48.0*r6*r2*(0.5-r6));
  }
		
  // Calc r+dv/dr = -24*(2*r^-12-r^-6). Argument r. 
  static inline double potential_w(double r)
  {
    double r6=1.0/pow(r,6);  // r^-6
    return(24*r6*(1-2*r6));
  }

};

int main()
{
  return 0;
}
