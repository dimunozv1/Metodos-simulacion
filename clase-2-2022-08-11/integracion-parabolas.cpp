#include <iostream>
 #include <cmath>

using namespace std;

double f(double x){
  return cos(x);
}

double Integral_Simpson(double a, double b, int N);
const double ErrMax=1e-7;

int main(){
  double a=0, b=M_PI/2;
  int N=50; //En este metodo el numero de particiones debe ser par,
  //por esto se multiplico por 2 en la funcion
  
 
  cout<<"La integral, por Simpson, es "<<Integral_Simpson(a,b,N)<<endl;
  
  return 0;
}


double Integral_Simpson(double a, double b, int N){ 
  
 int n=2*N;//numero de pasos o divisiones que se van a usar
 double suma=0;//aqui se van a guardar las sumas de la integral
 double h=(b-a)/n; //ancho de cada seccion 
    for(int i=0;i<=n;i++){
      if(i==0){
	suma+=f(a);
      }
      if(i==n){
	suma+=f(a+i*h);
      }
      else if(i%2==0){
	suma+=2*f(a+i*h);

      }
      else{
	suma+=4*f(a+i*h);
      }

    }
 double integral=suma*h/3;

  

  return integral;
}
