#include <iostream>
#include <cmath>
using namespace std;

const double w = 1;

double f1(double r, double x1, double x2)
{
    return -w*w*x2;
}
double f2(double r, double x1, double x2)
{
    return x1;
}

void pasoRK4(double& t, double& x1, double& x2, double dt)
{
  double dx11, dx21, dx31, dx41, dx1,dx2;                                        double dx12, dx22, dx32, dx42;
    dx11 = dt * f1(t,          x1,            x2);                        dx12 = dt * f2(t, x1, x2);
    dx21 = dt * f1(t + dt / 2, x1 + dx11 / 2, x2 + dx12 / 2);            dx22 = dt * f2(t + dt / 2, x1 + dx11 /2 , x2 + dx12 / 2);
    dx31 = dt * f1(t + dt / 2, x1 + dx21 / 2, x2 + dx22 / 2);            dx32 = dt * f2(t + dt / 2, x1 + dx21 /2 , x2 + dx22 / 2);
    dx41 = dt * f1(t + dt, x1 + dx31, x2 + dx32);                        dx42 = dt * f2(t + dt, x1 + dx31, x2 + dx32);

   
     dx1= (dx11 + 2 * (dx21 + dx31) + dx41) / 6;                           dx2 = (dx12 + 2 * (dx22 + dx32) + dx42) / 6;
     // cout<<dx41<<"\t"<<dx42<<endl;
     x1+=dx1;                                                              x2+=dx2;
    t += dt;
    

}

int main()
{
  double x1,x2, t, dt = 0.1;

  for (t = 0, x1=0, x2 =1 ; t < 2;)
    {
      
       cout << t << "\t" << x2 <<"\t"<<cos(t)<<  endl;
      
       pasoRK4(t, x1, x2, dt);

    }
    return 0;
}
