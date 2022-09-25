#include <cmath>
#include <iostream>

using namespace std;
double f(double t,double x){
  
  return t;
}

void PasoEuler(double & t,double & x, double dt){
  double dx;
  dx=dt*f(t,x);
  x+=dx; t+=dt;
}
  
const double Err= 1.0e-07;

  int main(void) {

    double t,x;
    double dt=0.1;

    for(t=0,x=0;t<2;){
      cout<<t<<"\t"<<x<<"\t"<<t*t/2<<endl;
      PasoEuler(t,x,dt);

    }



    return 0;
} 
