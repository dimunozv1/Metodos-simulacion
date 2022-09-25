#include <cmath>
#include <iostream>

using namespace std;
//constantes globales

const double g=9.8; //la gravedad del asunto
const double Gm=1;

//Declaracion de clases
class Cuerpo{  //datos internos que almacena, deben ser privados idealmente
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;

  //la parte publica son las ordenes que le pueden dar a las clases, todos deben
  //saber que ordenes le pueden dar al cuerpo

public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0); //solo da valores iniciales, no necesita devolver nada
  void CalculeFuerza(void);
  void Muevase(double dt); //para que se mueva si necesitamos el tiempo
  double Getx(void){return x;}; //Esto se conoce como una funcion inline, es mas eficiente porque no toca buscar la info en el stack
  double Gety(void){return y;};//Inline

};
//---------Funciones Globales--------
//se menciona a que objeto pertenece la orden
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  x=x0;
  y=y0;
  Vx=Vx0;
  Vy=Vy0;
  m=m0;
  R=R0;
  }

void Cuerpo:: CalculeFuerza(void){
  Fx=-Gm*m*x/pow(x*x+y*y,1.5); Fy=-Gm*m*y/pow(x*x+y*y,1.5);

}

void Cuerpo::Muevase(double dt){
  x+=Vx*dt; y+=Vy*dt;
  Vx+=Fx/m*dt; Vy+=Fy/m*dt;
}

//Al usar esta orden, se inicializa el balon con valores iniciales

  int main(void) {

    Cuerpo Planeta; //Es un instance, que es un ejemplar de la clase cuerpo
    double omega;
    double t,dt=0.0001;
    double r=10;
    double m=1;
    omega=sqrt(Gm/pow(r,3));
    double V0=omega*r;
    double T=2*M_PI/omega;

    Planeta.Inicie(r,0,0,V0/2,m,r);

    for(t=0;t<1.1*T;t+=dt){
      Planeta.CalculeFuerza();
      Planeta.Muevase(dt);
      cout<<Planeta.Getx()<<"\t"<<Planeta.Gety()<<endl;
    }

    
    

    

 
  return 0;
} 
