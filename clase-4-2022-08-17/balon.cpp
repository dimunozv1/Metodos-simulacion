#include <cmath>
#include <iostream>

using namespace std;
//constantes globales

const double g=9.8; //la gravedad del asunto

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
  Vy=y0;
  m=m0;
  R=R0;
  }

void Cuerpo:: CalculeFuerza(void){
  Fx=0; Fy=-m*g;

}

void Cuerpo::Muevase(double dt){
  x+=Vx*dt; y+=Vy*dt;
  Vx+=Fx/m*dt; Vy+=Fy/m*dt;
}

//Al usar esta orden, se inicializa el balon con valores iniciales

  int main(void) {

    Cuerpo Balon; //Es un instance, que es un ejemplar de la clase cuerpo
    double t,dt=0.1;

    Balon.Inicie(0,0,20,15,0.453,0.15);

    for(t=0;t<3.2;t+=dt){
      Balon.CalculeFuerza();
      Balon.Muevase(dt);
      cout<<Balon.Getx()<<"\t"<<Balon.Gety()<<endl;
    }

    

 
  return 0;
} 
