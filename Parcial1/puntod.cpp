// Simular el movimiento de N granos con gravedad, fricción y choques inelásticos
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double g=9.8,  Gamma=M_PI/9,Ca=0.5,p=1.224,Cm=1.0;
const int N=1;



const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F,omega,tau,theta; double m,R; double I;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0,
	      double theta0x, double theta0y,double theta0z, double omega0x,double omega0y,double omega0z);
  void BorreFuerza(){F.load(0,0,0); tau.load(0,0,0);};
  void AdicioneFuerza(vector3D F0,vector3D tau0){F+=F0;tau+=tau0;};
  void AdicioneTorque(vector3D tau0){tau+=tau0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double GetthetaX(void){return theta.x();}; //inline
  friend class Colisionador;
};
void Cuerpo:: Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0,
	      double theta0x, double theta0y,double theta0z, double omega0x,double omega0y,double omega0z){
  r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0);  m=m0;  R=R0;
  theta.load(theta0x,theta0y,theta0z);
  omega.load(omega0x,omega0y,omega0z); I=2.0/3*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);  omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<r.x()<<"+"<<R*cos(theta.z())/7<<"*t,"<<r.y()<<"+"<<R*sin(theta.z())/7<<"*t";
}

//--- clase Colisionador ----
class Colisionador{
private:

public:
  void Inicie(void);
  void CalculeFuerzas(Cuerpo * Balon,double dt);
 
 
};




void Colisionador::CalculeFuerzas(Cuerpo * Balon,double dt){
  int i,j; vector3D Fg,Fa,Fm,wvec,Vvec,taux; double Vnorm,zita;
  taux.load(0,0,0);
  double A=M_PI*pow(0.22,2);
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++){
    Balon[i].BorreFuerza();}
  //--- Sumar el peso ---
  for(i=0;i<N;i++){
    Fg.load(0,0,-Balon[i].m*g);
    Balon[i].AdicioneFuerza(Fg,taux);
  }

  //--Fuerza de arrastre
 for(i=0;i<N;i++){
   Vnorm=Balon[i].V.norm();
   Vvec=Balon[i].V;
   Fa= Vvec*Vnorm*Ca*p*A*-0.5;
   Balon[i].AdicioneFuerza(Fa,taux);}


 //--Fuerza de Magnus
 for(i=0;i<N;i++){
   wvec=Balon[i].omega;
   Vvec=Balon[i].V;
   zita=angle(wvec,Vvec);
   
   Fm= wvec^Vvec*-0.5*Cm*p*A*0.22/sin(zita);
   Balon[i].AdicioneFuerza(Fm,taux);}
}


 




//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Puntoc.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1:14]"<<endl;
  cout<<"set yrange[-7:2]"<<endl;
  cout<<"set title \"Tiro parabolico con Fg, Fa y Fm visto desde arriba\" "<<endl; 
  cout<<"set xlabel \"x\" "<<endl;
  cout<<"set ylabel \"y\" "<<endl; 
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

//-----------  Programa Principal --------------  
int main(void){
  
  Cuerpo Balon[N];
  Colisionador Newton;
  double m0=0.43, R0=0.22, V0=20;
  double omega0[6]={1.0*M_PI/5,2*M_PI/5,M_PI,2*M_PI,4*M_PI,10*M_PI};
  int i;
  double Err=5.0e-03;
  
  double t,tdibujo,tmax=2*V0*sin(Gamma)/g,tcuadro=tmax/200,dt=1e-3;
 
  
  // InicieAnimacion(); //Dibujar

  //Inicializar el Balon
  for(int j=0; j<6;j++){
      //----------(   x0,y0,z0, Vx0, Vy0,    Vz0,      m0,R0,theta0x,theta0y,theta0z,omega0x...)
  Balon[0].Inicie(0,0,0,V0*cos(Gamma),0,V0*sin(Gamma), m0,R0,0,0,0,0,0,omega0[j]);
    
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(abs(Balon[0].Getx()-5.0)<Err){
      
      cout<<abs(Balon[0].Gety())<<"\t"<<omega0[j]/(2*M_PI)<<endl;
    }

    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Balon[i].Mueva_r(dt,epsilon);
    Newton.CalculeFuerzas(Balon,dt);
    for(i=0;i<N;i++)Balon[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Balon[i].Mueva_r(dt,chi);
    Newton.CalculeFuerzas(Balon,dt);
    for(i=0;i<N;i++)Balon[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Balon[i].Mueva_r(dt,chiepsilon);
    Newton.CalculeFuerzas(Balon,dt);
    for(i=0;i<N;i++)Balon[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Balon[i].Mueva_r(dt,chi);
    Newton.CalculeFuerzas(Balon,dt);
    for(i=0;i<N;i++)Balon[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Balon[i].Mueva_r(dt,epsilon);  

  }   
  }
  
  return 0;
}

  
