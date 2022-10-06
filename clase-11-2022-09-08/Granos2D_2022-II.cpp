// Simular el movimiento de N granos con gravedad, fricción y choques inelásticos
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double g=9.8, K=1.0e4, Gamma=50, Kcundall=500, MU=0.4;
const double Lx=100, Ly=100;
const int Nx=2, Ny=2, N=Nx*Ny;

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
  vector3D r,V,F; double m,R; double theta,omega,tau; double I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
	      double theta0,double omega0);
  void BorreFuerza(){F.load(0,0,0); tau=0;};
  void AdicioneFuerza(vector3D F0,double tau0){F+=F0;tau+=tau0;};
  void AdicioneTorque(double tau0){tau+=tau0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double Gettheta(void){return theta;}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
		    double theta0,double omega0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);  omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";
}

//--- clase Colisionador ----
class Colisionador{
private:
  double xCundall[N+4][N+4];
  double sold[N+4][N+4];
public:
  void Inicie(void);
  void CalculeFuerzas(Cuerpo * Grano,double dt);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2, double &x_cundall,double &s_old,double t) ;
 
};

void Colisionador::Inicie(void){
  for(int i=0;i<N+4;i++){
    for(int j=0; j<N+4;j++){
      xCundall[i][j]=0;
      sold[i][j]=0;
    }
  }

}


void Colisionador::CalculeFuerzas(Cuerpo * Grano,double dt){
  int i,j; vector3D Fg;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N+4;i++)
    Grano[i].BorreFuerza();
  //--- Sumar el peso ---
  for(i=0;i<N;i++){
    Fg.load(0,-Grano[i].m*g,0);
    Grano[i].AdicioneFuerza(Fg,0);
  }
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++)
      CalculeFuerzaEntre(Grano[i], Grano[j],xCundall[i][j],sold[i][j],dt);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,double &x_Cundall,double &s_old,double dt){
  vector3D r21=Grano2.r-Grano1.r;
  double d=r21.norm();
  double s=Grano1.R+Grano2.R-d;
  double m12=(Grano1.m*Grano2.m)/(Grano1.m*Grano2.m);
  if(s>0){
    vector3D RW;
   
   
    //Vectores Unitarios
    vector3D n=r21*(1.0/d);
    vector3D t;
    t.load(n.y(),-n.x(),0);
    
    //Velocidad relativa
    vector3D V21=Grano2.V-Grano1.V;

    RW.load(0,0,Grano1.R*Grano1.omega+Grano2.R*Grano2.omega);
    vector3D Vc=V21-RW^n;
    double Vn=Vc*n;
    double Vt=Vc*t;
    
    //Fuerza de Hertz-Kuwara-Kono   
    double Fn=(K*pow(s,1.5))-Gamma*m12*sqrt(s)*Vn;
   

    //Fuerza de Cundall
    x_Cundall+=Vt*dt;
    double Ft;
    double Fk=-Kcundall*x_Cundall;
    double Fr=-MU*Fn;
    if(fabs(Fk)>fabs(Fr)){
      Ft=Fk;}
    else{
      Ft=Fk;
    }

    //Calcule y cargue las fuerzas
    
    vector3D tau1,k,F1,F2,tau2;
    k.load(0,0,1);
    F2=n*Fn+t*Ft;
    F1=(-1)*F2;
    tau2=(-Grano2.R*n)^F2;
    tau1=(Grano1.R*n)^F1;
    Grano1.AdicioneFuerza(F1,tau1*k);   Grano2.AdicioneFuerza(F2,tau2*k);
    
  }
  if(s_old>0 && s<0){
    x_Cundall=0;
  }
  s_old=s;
}


 




//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  // cout<<"set output 'Gas2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Grano[N+4];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1, R0=2, kT=10, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,tmax=10*(Lx/V0),tcuadro=tmax/1000,dt=1e-3;
  double dx=Lx/(Nx+1), dy=Ly/(Ny+1);
  double Theta,OmegaMax=1.0;
  
  InicieAnimacion(); //Dibujar

  //Inicializar las paredes
  double Rpared=100*Lx, Mpared=100*m0;
  //---------------(  x0,       y0,Vx0,Vy0,    m0,    R0, theta0,omega0) 
  Grano[N+0].Inicie(Lx/2,Ly+Rpared,  0,  0,Mpared,Rpared,      0,     0); //Pared de arriba
  Grano[N+1].Inicie(Lx/2,  -Rpared,  0,  0,Mpared,Rpared,      0,     0); //Pared de abajo
  Grano[N+2].Inicie(Lx+Rpared,Ly/2,  0,  0,Mpared,Rpared,      0,     0); //Pared derecha
  Grano[N+3].Inicie(  -Rpared,Ly/2,  0,  0,Mpared,Rpared,      0,     0); //Pared izquierda
  //Inicializar las moléculas
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //--------------------(   x0,   y0,          Vx0,          Vy0, m0,R0,theta0,omega0)
      Grano[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy,V0*cos(Theta),V0*sin(Theta), m0,R0,0,OmegaMax);//OJO
    }
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }

    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);  

  }   

  
  return 0;
}

  
