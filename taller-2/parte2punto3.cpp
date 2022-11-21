#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <exception>
#include <memory>


using namespace std;

const int Lx=256*2;
const int Ly=64;

const int Q=9; // q es el numero de flechas
const int N=800;


const double Ufan0=0.1;
const double tau=1.55;
const double utau=1.0/tau;
const double um_utau=1-utau;
const double eta=(tau-0.5)/3.0;

const double omega=2*M_PI/1000;
const double P=1;



//----------------Clase de la simulacion LB---------------//

class LatticeBoltzmann
{
private:
double Rc=8;
int ixc=128;
int iyc=32;
double w[Q];
int Vx[Q], Vy[Q];
double *f, *fnew; // where the asterisk is used for making
// f and fnew be crated throught the use
// of dynamic memory
double dA_dir[N*2];
double theta=2*M_PI/N;
double g=2*Rc*tan(theta/2);

public:
LatticeBoltzmann(void) // constructor function
{
w[0]=4.0/9;

w[1]=w[2]=w[3]=w[4]=1.0/9;


w[5]=w[6]=w[7]=w[8]=1.0/36;
w[5]=w[6]=w[7]=w[8]=1.0/36;

Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;


Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;



Vx[5]=1; Vx[6]=-1; Vx[7]=-1; Vx[8]=1;

Vy[5]=1; Vy[6]=1; Vy[7]=-1; Vy[8]=-1;



int ArraySize=Lx*Ly*Q;
f=new double [ArraySize]; fnew= new double [ArraySize];// use od dynamic memory


for(int i=0; i<N;i++)
{
dA_dir[2*i]=cos(i*theta);
dA_dir[2*i+1]=sin(i*theta);
//cout<<"x: "<<dA_dir[2*i]<<" \t"<<" y: "<<dA_dir[2*i+1]<<"\n";
}
};
~LatticeBoltzmann(void) // destructor function
{
delete[] f; delete[] fnew;
}
int n(int ix, int iy, int i)
{
return (ix*Ly+iy)*Q+i;
};
double rho(int ix, int iy, bool UseNew)
{
double sum;
int i, n0;
for (sum=0, i=0; i<Q; i++)
{
n0=n(ix,iy,i);
if(UseNew)
{sum+=fnew[n0];}
else
{sum+=f[n0];}
}

return sum;

}
double J(int ix, int iy, bool UseNew, bool Is_x)// hace el valor del momento j para una celda XY
{ // La variable Usenew reporta si usamos fnew o no
// La variable Is_x decide si la componente de J es x o y
double sum;
int i, n0;
for(sum=0,i=0; i<Q; i++)
{
n0=n(ix,iy,i);
if(Is_x)
{ if(UseNew)
{sum+=Vx[i]*fnew[n0];}
else
{sum+=Vx[i]*f[n0];}
}
else
{
if(UseNew)
{sum+=Vy[i]*fnew[n0];}
else
{sum+=Vy[i]*f[n0];}
}
}
return sum;
}
double f_eq(double rho0, double Ux0, double Uy0, int i)
{
double UdotV=Ux0*Vx[i]+ Uy0*Vy[i];
double U2=Ux0*Ux0+Uy0*Uy0;
return rho0*w[i]*(1+3*UdotV+4.5*UdotV*UdotV-1.5*U2);
}
double sigmaxx(double rho0, int ix, int iy)
{
double sum=0,Ux;
int i;
for(i=0; i<Q; i++)
{
Ux=J(ix+Vx[i],iy+Vy[i],false,true)/rho0;
sum+=w[i]*Vx[i]*Ux;
}
return -rho0/3.+eta*3*2*sum;
}
double sigmayy(double rho0, int ix, int iy)
{
double sum=0,Uy;

int i;
for(i=0; i<Q; i++)
{
Uy=J(ix+Vx[i],iy+Vy[i],false,false)/rho0;
sum+=w[i]*Vy[i]*Uy;
}
return -rho0/3.+eta*3*2*sum;

}
double sigmaxy(double rho0, int ix, int iy)
{
double sumx=0,sumy=0,Ux,Uy;
int i;
for(i=0; i<Q; i++)
{
Uy=J(ix+Vx[i],iy+Vy[i],false,false)/rho0;
Ux=J(ix+Vx[i],iy+Vy[i],false,true )/rho0;
sumy+=w[i]*Vy[i]*Uy;
sumx+=w[i]*Vx[i]*Ux;
}

return eta*3*(sumy+sumx);

}
double phi(double x, double y)
{
int ix=int(x),iy=int(y);
double u=(x-ix), v=(y-iy);
double sum=0;
double rhoxx=rho(ix,iy+1,false),rhoyy=rho(ix+1,iy,false);
double rhoxy=rho(ix+1,iy+1,false), rhoyx=rho(ix,iy,false);

sum+=sigmaxx(rhoxx,ix,iy+1)*(1-u)*v;
sum+=sigmayy(rhoyy,ix+1,iy)*u*(1-v);
sum+=sigmaxy(rhoxy,ix+1,iy+1)*u*v;
sum+=sigmaxy(rhoyx,ix,iy)*(1-u)*(1-v);
return sum;
}
void forces(void)
{ double Px,Py;
double sumx=0, sumy=0;
double stress;
double Fxp, Fyp;
int R=Rc;
double R2=R*R;
double Magnus=-0.5*P*M_PI*R2*R*Ufan0*omega;

for(int i=0;i<N;i++)
{
Px=Rc*dA_dir[2*i ]+ixc;
Py=Rc*dA_dir[2*i+1]+iyc;
stress=phi(Px,Py);
Fxp=stress*dA_dir[2*i]*g;
Fyp=stress*dA_dir[2*i+1]*g;
sumx+=Fxp;
sumy+=Fyp;

}
cout<<"Fuerzas en x: "<<sumx<<"\t"<<"Fuerzas en y:"<<sumy<<"\n";
 cout<<"Magnus"<< "\t"<< Magnus<< endl;

   

}
void Start(double rho0, double Ux0, double Uy0)
{
int ix,iy,i,n0;
for(ix=0;ix<Lx;ix++)
{for(iy=0;iy<Ly;iy++)
{for(i=0; i<Q; i++)
{
n0=n(ix,iy,i);
f[n0]=f_eq(rho0,Ux0,Uy0,i);
}
}
}
}
void Collision(void)
{
int ix, iy,i,n0;
double rho0,Ux0,Uy0;
for(ix=0;ix<Lx;ix++)
{for(iy=0;iy<Ly;iy++)
{
rho0=rho(ix,iy,false);
Ux0=J(ix,iy,false,true)/rho0;
Uy0=J(ix,iy,false,false)/rho0;
for(i=0;i<Q;i++)
{
n0=n(ix,iy,i);
fnew[n0]=um_utau*f[n0]+utau*f_eq(rho0,Ux0,Uy0,i);
}
}

}

}
void fields(double Ufan)// impose conditions for the problem
{ // Ufan es the velocity U of the fluid emitted by the fan of the wind tunel
int i,ix,iy,n0;
double lambda, rho0, Jx0,Jy0;

//--obstacle--//
 int R=Rc;
 double R2=R*R;

double Ux0,Uy0;
    
for(ix=0;ix<Lx;ix++)
{for(iy=0;iy<Ly;iy++)
{
rho0=rho(ix,iy,false);
if(ix==0) // If it's on the left border, then be a fan
{for(i=0;i<Q;i++)
{ n0=n(ix,iy,i);
fnew[n0]=f_eq(rho0,Ufan,0,i);
}

}
else if( (ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc) <= R2 )
  { Ux0=-omega*(iy-iyc);
    Uy0=omega*(ix-ixc);
    
for(i=0;i<Q;i++)
{
n0=n(ix,iy,i);
 fnew[n0]=f_eq(rho0,Ux0,Uy0,i);
}
}
else if(ix==ixc && iy==iyc+R+1) // this adds a small imperfection on the obstacle, to make it asymetric
{ // It is important for it to be asymetric,
for(i=0;i<Q;i++) // othewise we would never see the instability of Von-karman
{
n0=n(ix,iy,i);
fnew[n0]=f_eq(rho0,0,0,i);
}
}

}

}
}
void advection(void)
{
int ix,iy,i,ixnew,iynew,n0,n0new;
for(ix=0;ix<Lx;ix++)
{for(iy=0;iy<Ly;iy++)
{for(i=0;i<Q;i++)
{
ixnew=(ix+Vx[i]+Lx)%Lx;// the modulo L makes periodic boundary conditions
iynew=(iy+Vy[i]+Ly)%Ly;
n0=n(ix,iy,i);
n0new=n(ixnew,iynew,i);
f[n0new]=fnew[n0];

}

}

}

}
void Show(const char * NameFile, double Ufan)
{ ofstream salida(NameFile);
int ix,iy;
double rho0, Ux0,Uy0;

for(ix=0;ix<Lx;ix+=1)
{for(iy=0;iy<Ly;iy+=1)
{
rho0=rho(ix,iy,true);
Ux0=J(ix,iy,true,true)/rho0;
Uy0=J(ix,iy,true,false)/rho0;
salida<<ix<<" "<<iy<<" "<<Ux0/Ufan*2<<" "<<Uy0/Ufan*2<<"\n";
}
salida<<"\n";
}
salida.close();
}
};


//----------------Correr el prgrama----------------//


int main(void)
{
LatticeBoltzmann Aire;
int t,tmax=1000;//=10000;
double rho0=1.0;

//Start
Aire.Start(rho0,Ufan0,0);
//Run
for(t=0;t<tmax;t++){
Aire.Collision();
Aire.fields(Ufan0);
Aire.advection();
}
Aire.Collision();
Aire.fields(Ufan0);
Aire.advection();
Aire.forces();

//Print
Aire.Show("punto_a.dat",Ufan0);
return 0;
} 
