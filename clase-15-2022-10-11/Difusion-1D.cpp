#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=1024;
const double p=0.5;

const int Q=2; //Numero de flechas

//-------Clase LatticeGas--------

class LatticeGas{

private:
  int V[Q]; //Q=0 (derecha) Q=1 (izquierda), esto me da las direcciones
  int n[Lx][Q], nnew[Lx][Q];    //es 2 porque pueden haber dos direcciones
public:
  LatticeGas(void);
  void Inicie(int N,double mu,double sigma,Crandom &ran64);
  void Borrese(void);
  void Show(void);
  void Colisione(Crandom &ran64);
  void Adveccione(void);
  double rho(int x);
  double Varianza(void);
  void GrafiqueRho(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1;
  V[1]=-1;
}
void LatticeGas::Inicie(int N,double mu,double sigma,Crandom &ran64){
  int ix,i;
  while(N>0){
    ix=(int)ran64.gauss(mu,sigma);
    if(ix<0){
      ix=0;}
    if(ix>Lx-1){
      ix=Lx-1;}
    i=(int)Q*ran64.r();//Escoger direccion al azar
    
    if(n[ix][i]==0){
      n[ix][i]=1;N--;
     
     
    }
      }
}
void LatticeGas::Borrese(void){
  for(int ix=0;ix<Lx;ix++){
    for(int i=0;i<Q;i++){
      n[ix][i]=nnew[ix][i]=0;
    }
  }


}
void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){
    for(int ix=0;ix<Lx;ix++){
      cout<<n[ix][i];
    }
    cout<<endl; //Queremos primero imprirmir todos los que va hacia un lado y luego hacia otro
  }
}

void LatticeGas::Colisione(Crandom &ran64){
  for(int ix=0;ix<Lx;ix++){
    if(ran64.r()>p){
      nnew[ix][0]=n[ix][1];
      nnew[ix][1]=n[ix][0];
    }
    else{
      nnew[ix][0]=n[ix][0];
      nnew[ix][1]=n[ix][1];
    }
      
  }
}
void LatticeGas::Adveccione(void){
   for(int ix=0;ix<Lx;ix++){
     for(int i=0;i<Q;i++){
       n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i];}
   }

}

double LatticeGas::rho(int x){
  return n[x][0]+n[x][1];
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

  Sigma2/=(N-1);

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

  Difusion.Borrese();
 
  Difusion.Inicie(N,mu,sigma,ran64);
 
  for(int t=0;t<tmax;t++){   
    
    cout<<t<<"\t"<<Difusion.Varianza()<<endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
   
  }


  return 0;
}
