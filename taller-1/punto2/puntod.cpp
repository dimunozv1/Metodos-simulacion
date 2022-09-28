#include <iostream>
#include <cmath>
#include <cstdlib>


using namespace std;


const double Err= 1.0e-011;

double f(double t,int m, double r){
   return cos(m*t-r*sin(t));
  // return cos(t);
}


double integral_besel(double x, double alpha, int n,double a,double b)
{
	double t, h, sum;
	int i;
	n = n * 2;
	h = (b - a) / n;
	for (i = 0, sum = 0; i <= n; i += 1)
	{
		t = a + i * h;
		if (i == 0 || i == n)
		{
			sum += f(t,alpha,x);
		}
		else if (i % 2 == 0)
		{
			sum += 2 * f(t,alpha,x);
		}
		else
		{
		  sum += 4 * f(t,alpha,x);
		}
	}
	return sum * h / 3;

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

double Bessel(int m, double r){
  

  return (1/M_PI)* integral_besel(r, m, 100,0,M_PI);

}
double Biseccion(double a,double b){


  double lambda=(b+a)/2;
  while((b-a)>Err){
    lambda=(b+a)/2;
    
    if(Bessel(0,a)*Bessel(0,lambda)>=0){
      a=lambda;
    }
    if(Bessel(0,b)*Bessel(0,lambda)>=0){
      b=lambda;
    }
     }
  

     return lambda;
}
 int main(int argc,char** argv)
{
   int a=atoi(argv[1]);
   if(a==0){
   for(double r=0.1;r<16;r+=0.1){
     cout<<r<<"\t"<< Bessel(0,  r)<<endl;
   }} //Datos para graficar la funcion de bessel de orden 0
  double rango0[]={2.0,4.0,4.0,6.0,7.0,10.0,10.0,12.0};
  double ceros[4];
  int j=0;
  if(a==1){  //imprimir los ceros de la funcion de bessel de orden 0
    // cout<<"Los ceros de la funcion de bessel son:"<<endl;
   for(int i =0;i<7;i++){
     if(i%2==0){
       ceros[j]=Biseccion(rango0[i],rango0[i+1]);
       cout<<ceros[j]<<"\t"<<0<<endl;
       j+=1;
     }
   }
    }
 
    return 0;
}
