#include <iostream>
#include <cmath>
#include <cstdlib>


using namespace std;


const double Err= 1.0e-09;
double f1(double r, double x1, double x2,double Lambda)
{
    return -x1/r-Lambda*Lambda*x2;
}
double f2(double r, double x1, double x2)
{
    return x1;
}

void pasoRK4(double& t, double& x1, double& x2, double dt,double Lambda)
{
    double dx11, dx21, dx31, dx41;                                        double dx12, dx22, dx32, dx42;
    dx11 = dt * f1(t,          x1,            x2,Lambda);                        dx12 = dt * f2(t, x1, x2);
    dx21 = dt * f1(t + dt / 2, x1 + dx11 / 2, x2 + dx12 / 2,Lambda);            dx22 = dt * f2(t + dt / 2, x1 + dx11 /2 , x2 + dx12 / 2);
    dx31 = dt * f1(t + dt / 2, x1 + dx21 / 2, x2 + dx22 / 2,Lambda);            dx32 = dt * f2(t + dt / 2, x1 + dx21 /2 , x2 + dx22 / 2);
    dx41 = dt * f1(t + dt, x1 + dx31, x2 + dx32,Lambda);                        dx42 = dt * f2(t + dt, x1 + dx31, x2 + dx32);


     x1 += (dx11 + 2 * (dx21 + dx31) + dx41) / 6;                            x2 += (dx12 + 2 * (dx22 + dx32) + dx42) / 6;
    t += dt;
   

}
double f(double Lambda,int flag)
{
  double x1,x2, r, dr = 0.01;
  x1=0;
   x2=1;
  
    for (r = 0.01; r < 1.0;)
    {
      
      if(flag==1){cout << r << "\t" << x2 <<  endl;}
      
      pasoRK4(r, x1, x2, dr,Lambda);

    }
    return x2;
      // cout << Lambda << "\t" << x2 <<  endl;;
}
double Biseccion(double a,double b){


  double lambda=(b+a)/2;
  while((b-a)>Err){
    lambda=(b+a)/2;
    
    if(f(a,0)*f(lambda,0)>=0){
      a=lambda;
    }
    if(f(b,0)*f(lambda,0)>=0){
      b=lambda;
    }
     }
  

     return lambda;
}
 int main(int argc,char** argv)
{
 
  double abLambda[]={2.0,4.0,4.0,6.0,7.0,10.0,10.0,12.0};
  double Lambda0[4];
  int j=0;
  // cout<<Biseccion(-4,4)<<endl;
   for(int i =0;i<7;i++){
     if(i%2==0){
       Lambda0[j]=Biseccion(abLambda[i],abLambda[i+1]);
       j+=1;
     }
   }
  int a=atoi(argv[1]);
  if(a<0){
    cout<< "Use un valor entero entre 0 y 3"<<endl;
    return 0;}
   if(a>3){
    cout<< "Use un valor entero entre 0 y 3"<<endl;
    return 0;}
 
  f(Lambda0[a],1);
    
 
    return 0;
}
