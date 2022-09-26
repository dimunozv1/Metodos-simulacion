#include <iostream>
#include <cmath>
using namespace std;

const double Lambda = 1;

double f1(double r, double x1, double x2)
{
    return -x1/r-Lambda*Lambda*x2;
}
double f2(double r, double x1, double x2)
{
    return x1;
}

void pasoRK4(double& t, double& x1, double& x2, double dt)
{
    double dx11, dx21, dx31, dx41;                                        double dx12, dx22, dx32, dx42;
    dx11 = dt * f1(t,          x1,            x2);                        dx12 = dt * f2(t, x1, x2);
    dx21 = dt * f1(t + dt / 2, x1 + dx11 / 2, x2 + dx12 / 2);            dx22 = dt * f2(t + dt / 2, x1 + dx11 /2 , x2 + dx12 / 2);
    dx31 = dt * f1(t + dt / 2, x1 + dx21 / 2, x2 + dx22 / 2);            dx32 = dt * f2(t + dt / 2, x1 + dx21 /2 , x2 + dx22 / 2);
    dx41 = dt * f1(t + dt, x1 + dx31, x2 + dx32);                        dx42 = dt * f2(t + dt, x1 + dx31, x2 + dx32);


     x1 += (dx11 + 2 * (dx21 + dx31) + dx41) / 6;                            x2 += (dx12 + 2 * (dx22 + dx32) + dx42) / 6;
    t += dt;
   
}

int main()
{
  double x1,x2, r, dr = 0.01;

  for (r = 0.01, x1=0, x2 = 1; r < 10.0;)
    {
      
      cout << r << "\t" << x2 <<  endl;
      
        pasoRK4(r, x1, x2, dr);

    }
    return 0;
}
