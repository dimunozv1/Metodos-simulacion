#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

//const double Beta = 0.35;
const double Gamma = 0.008;
double f1(double t, double s, double i,double Beta)
{
    return -Beta * s * i;
}
double f2(double t, double s, double i, double Beta)
{
    return Beta * s * i - Gamma * i;
}

void pasoRK4(double& t, double& x1, double& x2, double dt,double Beta)
{
    double dx11, dx21, dx31, dx41;                                        double dx12, dx22, dx32, dx42;
    dx11 = dt * f1(t,          x1,            x2,Beta);                        dx12 = dt * f2(t, x1, x2,Beta);
    dx21 = dt * f1(t + dt / 2, x1 + dx11 / 2, x2 + dx12 / 2,Beta);            dx22 = dt * f2(t + dt / 2, x1 + dx11 /2 , x2 + dx12 / 2,Beta);
    dx31 = dt * f1(t + dt / 2, x1 + dx21 / 2, x2 + dx22 / 2,Beta);            dx32 = dt * f2(t + dt / 2, x1 + dx21 /2 , x2 + dx22 / 2,Beta);
    dx41 = dt * f1(t + dt, x1 + dx31, x2 + dx32,Beta);                        dx42 = dt * f2(t + dt, x1 + dx31, x2 + dx32,Beta);


     x1 += (dx11 + 2 * (dx21 + dx31) + dx41) / 6;                            x2 += (dx12 + 2 * (dx22 + dx32) + dx42) / 6;
    t += dt;
    //cout<<"hola   " <<f1(t,0.999,0)<<"\t"<<f2(t,0.999,0)<<endl;

}

int main()
{
  cout.precision(7);
  cout.setf(ios::scientific);
  double t, s, i, r, dt = 0.1;
  t = 0;
  s = 0.999;
  i = 0.001;
  r=1-s-i;
  double sold=0;
  double Err=1.0e-09;
  double Beta;

  for(Beta=0.01;Beta<1;Beta+=0.001){
       while(abs(sold-s)>Err)
    {
      sold=s;
      //cout << t << "\t" << s << "\t" << i <<"\t"<< r <<  endl;
      
      pasoRK4(t, s, i, dt,Beta);
	r=1-s-i;
    }
       cout<<s<<"\t"<<Beta/Gamma<<endl;
       s=0.999;
       t=0;
       i=0.001;
       r=0;}
    return 0;
}
