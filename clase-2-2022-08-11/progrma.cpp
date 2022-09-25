#include <iostream>
#include <cmath>

using namespace std;

double f(double x){
  return sin(x)/x;
}

double bisection(double a, double b);
const double ErrMax=1e-7;

int main(){
  double a=2, b=4;

  cout<<"El cero es "<<bisection(a,b)<<endl;
  
  return 0;
}


double bisection(double a, double b){ 
  double fa=f(a);
  double fm,m;
  
  while(b-a >= ErrMax){
    m = (b+a)/2; fm=f(m);
      if(f(a)*f(m)>0)
	{a=m; fa=fm;}
      else
      b=m;
    }

  double zero=(a+b)/2;

  return zero;
}
