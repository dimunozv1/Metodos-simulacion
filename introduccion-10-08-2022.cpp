#include <cmath>
#include <iostream>

double f(double x){
  
  return sin(x)/x;
}  

  int main(void) {

  double x;
  double a;
  double b;  
  for (x =0.1; x <= 10; x += 0.1) {
    std::cout << x << " " << sin(x) / x << std::endl;
  }}
  std::cout << x << std::endl;
  return 0;
} 
