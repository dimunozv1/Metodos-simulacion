#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes globales

const int N=4;
const double g=980;  //CGS
const double K=1e9;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  double theta,omega,tau;
  double m,R,l,x;
public:
  void Inicie(double theta0,
	      double omega0,double x0,double l0,double m0,double R0);
  void BorreTorque(void){tau=0;};
  void SumeTorque(double tau0);
  void Mueva_theta(double dt,double coeficiente);
  void Mueva_omega(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return x+l*sin(theta);}; //Inline
  double Gety(void){return -l*cos(theta);}; //Inline
  double Gettau(void){return tau;}; //Inline
  double Gettheta(void){return theta;};
  double Getomega(void){return omega;};
  friend class Colisionador;
};
void Cuerpo::Inicie(double theta0,
	      double omega0,double x0,double l0,double m0,double R0){
  omega=omega0;
  theta=theta0;
  l=l0;
  x=x0;
  m=m0;
  R=R0;
}
void Cuerpo::SumeTorque(double tau0){
  tau+=tau0;
}
  //cout<<"SumaTorque tau0"<<tau0<<" tau "<<tau<<endl;
void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(dt*coeficiente);
  // cout<<"Mueva_theta "<<theta<< " omega " << omega<<endl;
  
}
void Cuerpo::Mueva_omega(double dt,double coeficiente){
  omega+=tau*(dt*coeficiente/(m*l*l));
  //cout<<"Mueva_omega "<<omega<<"tau "<<tau<<endl;
  
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)  ";
  cout<<" , "<<x<<"+"<<l/7<<"*t*sin("<<theta<<"),-"<<l/7<<"*t*cos("<<theta<<")";
  
  
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeTorque(Cuerpo* Pendulo);
  void CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);    
};
void Colisionador::CalculeTorque(Cuerpo*Pendulo){
  double tau0;
  int i,j;
  //Borrar torques
  for(i=0;i<N;i++){
    Pendulo[i].BorreTorque();
    tau0=-Pendulo[i].l*Pendulo[i].m*g*sin(Pendulo[i].theta);
    Pendulo[i].SumeTorque(tau0);
  }
    //cout<<"tau"<<Pendulo[i].Gettau()<<endl;
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=N-1;i>0;i--){
    
      CalculeTorqueEntre(Pendulo[i],Pendulo[i-1]);
      // cout<<i<<endl;
  }
}
void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  double s=(Pendulo2.Getx()+Pendulo2.R)-(Pendulo1.Getx()-Pendulo2.R);
  double F=0;
  if(s>0){
    F=K*pow(s,1.5);}
  Pendulo1.SumeTorque(F*Pendulo1.l);
  Pendulo2.SumeTorque(-F*Pendulo1.l);
}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Pendulo.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-5:12]"<<endl;
  cout<<"set yrange[-15:2]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(){
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=50, m1=50;
  double l0=12;
  double R=1; //radio para dibujarlos
  // double M=m0+m1, x0=-m1*r/M, x1=m0*r/M;
  double T=2*M_PI/sqrt(g/l0);
  double tdibujo,tcuadro=T/100;
  int i;
  double dt=0.001;
  double tmax=3*T;
  double t;
  //---------(theta0,omega0,x0,l0,m0,R0)
  Pendulo[0].Inicie(-0.2, 0,0,l0,m0,R);
  for(i=1;i<N;i++){
    Pendulo[i].Inicie(0,0,2*R*i ,l0,m0,R);

  }
  
  
  // InicieAnimacion();
  
  for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
    
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<N;i++) Pendulo[i].Dibujese();
      TermineCuadro();
      tdibujo=0;
    } 
    // if(abs(Pendulo[1].Gettau())<1e05){
    //  cout<<t<<"\t"<<Pendulo[1].Gettau()<<endl;}
    
    //cout<<Pendulo[1].Getx()<<" "<<Pendulo[1].Gety()<<endl;
    // Mover por PEFRL
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);
    Newton.CalculeTorque(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorque(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Coeficiente2);
    Newton.CalculeTorque(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorque(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);   
  }
  
  return 0;
}
