#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;
//Constantes globales
const double G = 1.0;
const int N=2;

//constante PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);


//Declaración de las clases
class Cuerpo;
class Colisionador;
//---------- Clase Cuerpo --------------
class Cuerpo {
private:
  vector3D r, V, F;
    double  m, R;
public:
    void Inicie(double x0, double y0, double z0, double Vx0, double Vy0,
		double Vz0, double m0, double R0);
  void BorreFuerza(void){F.load(0,0,0);};
  void SumeFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt,double coeficiente);
    void Mueva_V(double dt,double coeficiente);
    void Dibujese(void);

     void Getr(void) { r.show(); }; //Inline
   void GetF(void) {  F.show(); };
    double Getx(void) { return r.x(); }; //Inline
    double Gety(void) { return r.y(); }; //Inline
    double Getz(void) { return r.z(); }; //Inline
    friend class Colisionador;
};



void Cuerpo::Inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0) {
  
    r.load(x0, y0, z0);
    V.load(Vx0, Vy0, Vz0);
    m = m0; R = R0;
}


void Cuerpo:: Mueva_r(double dt,double coeficiente){
  r+=V*(dt*coeficiente);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(dt*coeficiente/m);
}

void Cuerpo::Dibujese(void){
cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}


//----------Clase Colisionador------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo* Planeta);
  void CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo& Planeta2);
};

void Colisionador::CalculeFuerzas(Cuerpo* Planeta){

  int i,j;
  //Borrar fuerzas
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
	}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo& Planeta2){
  vector3D r21;
  double d21;
  vector3D n;
  r21=Planeta2.r-Planeta1.r;
  d21=r21.norm();
  // cout<<d21<<endl;
  n=r21/d21;
  double F=G*Planeta1.m*Planeta2.m*pow(d21,-2.0);
  vector3D F1=F*n;
  Planeta1.SumeFuerza(F1);
  Planeta2.SumeFuerza((-1)*F1);}

//-------------------------- Funciones de Animacion -------------------

void InicieAnimacion(void){
cout<<"set terminal gif animate"<<endl;
cout<<"set output 'Na.gif'"<<endl;
cout<<"unset key"<<endl;
cout<<"set xrange[-1500:1500]"<<endl;
cout<<"set yrange[-1500:1500]"<<endl;
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

//----------------------- Programa Principal ------------------------
  

  


int main(void) {
    Cuerpo Planeta[N];
    double m0 = 1047, m1 = 1,r=1000;
    double M=m0+m1,x0=-m1*r/M,x1=m0*r/M;
    double omega, V0, T,V1;
    double t,tmax, dt = 0.1;
    Colisionador Newton;
    
    
    omega = sqrt(G*M*pow(r,-3.0));
    V0 = omega * x0;
    T = 2 * M_PI / omega;
    V1=omega*x1;
    tmax=20*T;
    dt=0.1;
    double tdibujo,tcuadro=T/200;
    int i;

    //------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
    //double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0
    Planeta[0].Inicie(x0, 0, 0, 0, V0,0, m0, 100);
    Planeta[1].Inicie(x1, 0, 0, 0, V1 ,0, m1, 40 );
    //Comandos gnuplot
    InicieAnimacion();
    int flag=1;
   
    for (t = 0,tdibujo=0; t < tmax; t += dt,tdibujo+=dt) {
      
     
      if(flag==1){
      if(tdibujo>tcuadro){
	InicieCuadro();
	for(int i=0;i<N;i++){
	  Planeta[i].Dibujese();}
	TermineCuadro();
	tdibujo=0;
      }
      }
      // cout<<Planeta[0].GetF()<<"\t"<<Planeta[1].GetF()<<endl;
      //Mover por PEFRL
	
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_r(dt,Zeta);}

	Newton.CalculeFuerzas(Planeta);
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_V(dt,Coeficiente1);}
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_r(dt,Chi);}
	Newton.CalculeFuerzas(Planeta);
	for(i=0;i<N;i++){	
	  Planeta[i].Mueva_V(dt,Lambda);}
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_r(dt,Coeficiente2);}
        Newton.CalculeFuerzas(Planeta);
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_V(dt,Lambda);}
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_r(dt,Chi);}
	Newton.CalculeFuerzas(Planeta);
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_V(dt,Coeficiente1);}
	for(i=0;i<N;i++){
	  Planeta[i].Mueva_r(dt,Zeta);}
    }

    return 0;
}


