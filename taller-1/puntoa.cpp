#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//------Constantes Globales------
const double g=9.8;
const double Lx=60;
const double Ly=60;  //Tamanio de la caja
const int Nx=1;
const int Ny=1;
const int N=Nx*Ny;
const double K=1.0e04;
const double r0=10;
const double epsi=1.0;



const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//----Funciones----
void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);

//-------Declarar clases------

class Cuerpo;
class Colisionador;
class Cuerpo{

private:
  vector3D r,F,v;
  double m,R;

public:
  void Inicie(double x0,double y0,double vx0,double vy0,double m0,double R0);
  void SumeFuerza(vector3D F0){F+=F0;};
  void BorreFuerza(void){F.load(0,0,0);};
  void Mueva_r(double dt,double coeficiente);
  void Mueva_v(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double Getvx(void){return v.x();};
  double Getvy(void){return v.y();};
  friend class Colisionador;

};

class Colisionador{

private:

public:
  void CalculeFuerza(Cuerpo *Molecula); //porque sabemos que se va a utilizar un array de moleculas, se usa un pointer
  void CalculeFuerzaEntre(Cuerpo &Molecula1,Cuerpo  &Molecula2);
    };
int main(void){

  double m0=1.0;
  double R0=2.5;
  double x0=10;
  double y0=0;
  double kT=0.5;
  double V0=sqrt(2*kT/m0);
  double Vx0=V0;
  double Vy0=0;
  int tmin=0;
  int tmax=10*(Lx/V0);
  double dt=1e-3;
  double Rpared=100*Lx;
  double Mpared=100*m0;
  double tcuadro=tmax/1000;
  int dib=0;
  //cout<<"tcuadro "<<tcuadro<<endl;
    
  double t,tdibujo;
  int i;
  Cuerpo Molecula[N+4];
  Colisionador LJ;
  //Inicializar paredes
  Molecula[N+0].Inicie( -Rpared,Ly/2,0,0,100*m0,R0*100);//izquierda
  Molecula[N+1].Inicie(Lx+Rpared,Ly/2,0,0,100*m0,R0*100);//derecha
  Molecula[N+2].Inicie(Lx/2, -Rpared,0,0,100*m0,R0*100);//abajo  
  Molecula[N+3].Inicie(Lx/2,Ly+Rpared,0,0,100*m0,R0*100);
  //inicializar moleculas
 
  for(int i=0;i<N;i++){//x0, y0, vx0, vy0, m0, R0
    Molecula[i].Inicie(0,Ly/2,V0,0,m0,R0);
    // cout<<"inicio "<<Molecula[i].Getx()<<endl;
  }
  if (dib==1){
  InicieAnimacion();}

  for(t=0,tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if (dib==1){
    if(tdibujo>tcuadro){
         InicieCuadro();
      for(int i=0;i<N;i++){
	      Molecula[i].Dibujese();
        }
       TermineCuadro();
      tdibujo=0;
     
    }
    }
     //cout<<" x "<<Molecula[0].Getx()<<endl;
    //---  PEFRL ---
     for(i=0;i<N;i++){
       Molecula[i].Mueva_r(dt,epsilon);}
    //cout<<" 1x "<<Molecula[0].Getx()<<"      1v   "<<Molecula[0].Getvx()<<endl;
    LJ.CalculeFuerza(Molecula);
    // cout<<" x  after fuerza"<<Molecula[0].Getx()<<endl;
    for(i=0;i<N;i++){
      Molecula[i].Mueva_v(dt,lambda2);}
    //cout<<" 2x after v "<<Molecula[0].Getx()<<"      2v   "<<Molecula[0].Getvx() <<endl;
    for(i=0;i<N;i++){
      Molecula[i].Mueva_r(dt,chi);}
    //cout<<" 3x "<<Molecula[0].Getx()<<"      3v   "<<Molecula[0].Getvx()<<endl;
    LJ.CalculeFuerza(Molecula);
    for(i=0;i<N;i++){
      Molecula[i].Mueva_v(dt,lambda);}
    for(i=0;i<N;i++){
      Molecula[i].Mueva_r(dt,chiepsilon);}
    // cout<<" x "<<Molecula[0].Getx()<<endl;
    LJ.CalculeFuerza(Molecula);
    for(i=0;i<N;i++){
      Molecula[i].Mueva_v(dt,lambda);}
    for(i=0;i<N;i++){
      Molecula[i].Mueva_r(dt,chi);}
    //  cout<<" x "<<Molecula[0].Getx()<<endl;
    LJ.CalculeFuerza(Molecula);
    for(i=0;i<N;i++){
      Molecula[i].Mueva_v(dt,lambda2);}
    for(i=0;i<N;i++){
      Molecula[i].Mueva_r(dt,epsilon);}
    // cout<<" x "<<Molecula[0].Getx()<<endl;
      


    }





  return 0;

}

void Cuerpo::Inicie(double x0,double y0,double vx0,double vy0,double m0,double R0){
  r.load(x0,y0,0);
  v.load(vx0,vy0,0);
  m=m0;
  R=R0;    
}
void Cuerpo::Mueva_r(double dt,double coeficiente){

  r+=v*(coeficiente*dt);
}
void Cuerpo::Mueva_v(double dt,double coeficiente){

  v+=F*(coeficiente*dt/m);
  //cout<<"Fuerza : "<< F.x() <<endl;
}
void Cuerpo::Dibujese(void){
  cout<<","<< r.x() << "+" << R <<"*cos(t),"<<r.y() << "+" << R <<"*sin(t)";
}

//------Colisionador-----------

void Colisionador::CalculeFuerza(Cuerpo *Molecula){
  //Borrar Fuerza
  for(int i=0;i<N+4;i++){
    Molecula[i].BorreFuerza();
  }

  //Calcular fuerza
  for(int i=0;i<N;i++){
    for(int j=i+1;j<N+4;j++)
    {
      CalculeFuerzaEntre(Molecula[i],Molecula[j]);
    }
  }
}

void Colisionador::CalculeFuerzaEntre(Cuerpo &Molecula1,Cuerpo  &Molecula2){
   vector3D r21=Molecula2.r-Molecula1.r;
   double d21=r21.norm();
  
  double s=Molecula1.R+Molecula2.R-d21;
  cout<<s<<endl;
  if(s>0){
    vector3D n=r21*(1.0/d21);
    //Fuerza de Lennard-Jones
    
    // double F=12*epsi/d21*(pow(r0,12)/pow(d21,12)-pow(r0,6)/pow(d21,6));
    vector3D F2=n*(K*pow(s,1.5));
    Molecula2.SumeFuerza(F2);   Molecula1.SumeFuerza(F2*(-1));
  }

}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Gas2D.gif'"<<endl;
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
