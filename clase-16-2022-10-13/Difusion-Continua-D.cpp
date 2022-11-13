#include <cmath>
#include <iostream>
#include "Random64.h"
using namespace std;

const int Lx=1024;
const double p=0.5;

const int Q=2; //Numero de flechas

//-------Clase LatticeGas--------

class LatticeGas{

private:
  int V[Q]; //Q=0 (derecha) Q=1 (izquierda), esto me da las direcciones
  double f[Lx][Q], fnew[Lx][Q];    //es 2 porque pueden haber dos direcciones
public:
  LatticeGas(void);
  void Inicie(int N,double mu,double sigma);
  void Colisione(void);
  void Adveccione(void);
  double rho(int x);
  double Varianza(void);
  void GrafiqueRho(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1;
  V[1]=-1;
}
void LatticeGas::Inicie(int N,double mu,double sigma){
  for(int ix=0;ix<Lx;ix++){
    double rho=(N/(sigma*sqrt(2*M_PI)))*exp(-0.5*pow((ix-mu)/sigma,2.0));
  for(int i=0;i<Q;i++){
      f[ix][i]=rho/Q;
  }
  }
}

void LatticeGas::Colisione(void){
  int j;
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      j=(i+1)%Q;
      fnew[ix][i]=f[ix][i]+(1-p)*(f[ix][j]-f[ix][i]);
    }
   
      
  }
}
void LatticeGas::Adveccione(void){
   for(int ix=0;ix<Lx;ix++){
     for(int i=0;i<Q;i++){
       f[(ix+V[i]+Lx)%Lx][i]=fnew[ix][i];}
   }

}

double LatticeGas::rho(int ix){
  double suma=0; 
  for(int i=0;i<Q;i++){
    suma+=f[ix][i];}
  return suma;
}

double LatticeGas::Varianza(void){

  int ix;
  double N,Xprom,Sigma2;

  for(N=0,ix=0;ix<Lx;ix++){
    N+=rho(ix);
  }
  for(Xprom=0,ix=0;ix<Lx;ix++){
    Xprom+=ix*rho(ix);
  }
  
  Xprom/=N;

  for(Sigma2=0,ix=0;ix<Lx;ix++){
    Sigma2+=pow(ix-Xprom,2.0)*rho(ix);}

  Sigma2/=N;

  return Sigma2;
      
    

}

void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++){
    cout<<ix<<"\t"<<rho(ix)<<endl;}
}

//--------Programa Principal--------

int main(void){
  LatticeGas Difusion;
  Crandom ran64(1);
  int N=400;
  double mu=Lx/2;
  double sigma=Lx/8;
  int tmax=400;

  
 
  Difusion.Inicie(N,mu,sigma);
 
  for(int t=0;t<tmax;t++){   
    
    cout<<t<<"\t"<<Difusion.Varianza()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
   
  }


  return 0;
}
