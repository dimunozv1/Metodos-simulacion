#include <cmath>
#include <iostream>

double f(double x){
  
  return sin(x)/x;
}
const double Err= 1.0e-07;

  int main(void) {

  double x;
  double a=2;
  double b=4;
  double m=(b+a)/2;
  while((b-a)/2>Err){
    m=(b+a)/2;
    if(f(a)*f(m)>0){
      a=m;
    }
    if(f(b)*f(m)>0){
      b=m;
    }
     }
   std::cout<< m <<std::endl;
  for (x =0.1; x <= 10; x += 0.1) {
    // std::cout << x << " " << sin(x) / x << std::endl;
  }
//std::cout << x << std::endl;
  return 0;
} 
