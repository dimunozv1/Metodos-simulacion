#include <iostream>
 #include <cmath>

using namespace std;

double f(double t,int m, double x){
  return cos(m*t-x*sin(t));
}

double Integral_Simpson(double a, double b, int N, int m, double x);

double Bessel(int m, double x);
double bisection(double a, double b, int k);

const double ErrMax=1e-7;

int main(){
  double m=0, a=2,b=4;
  cout<<bisection(a,b,m)<<endl;
    
  
  
  return 0;
}


double Integral_Simpson(double a, double b, int N, int m, double x){ 
  
 int n=2*N;//numero de pasos o divisiones que se van a usar
 double suma=0;//aqui se van a guardar las sumas de la integral
 double h=(b-a)/n; //ancho de cada seccion
 double t;
    for(int i=0;i<=n;i++){
      t=a+i*h;
      if(i==0){
	suma+=f(t,m,x);
      }
      if(i==n){
	suma+=f(t,m,x);
      }
      else if(i%2==0){
	suma+=2*f(t,m,x);

      }
      else{
	suma+=4*f(t,m,x);
      }

    }
 double integral=suma*h/3;

  

  return integral;
}

double Bessel(int m, double x){

  double a =0;
  double b=M_PI;
  int N=50;

  return Integral_Simpson(a,  b, N, m, x)/M_PI;

}
 double bisection(double a, double b, int alpha){ 
  double fa=Bessel(alpha,a);
  double fm,m;
  
  while(b-a >= ErrMax){
    m = (b+a)/2; fm=Bessel(alpha,m);
      if(fa*fm>0)
	{a=m; fa=fm;}
      else
      b=m;
    }

  double zero=(a+b)/2;

  return zero;
}

