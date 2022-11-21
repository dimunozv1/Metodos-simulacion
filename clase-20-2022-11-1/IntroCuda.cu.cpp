#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

#define Lx 16 //esto es una orden al compilador
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//---------Programa del Device o GPU-----
//---------Kernels----------

__global__ void AddTwoVectors(float *d_a,float *d_b,float *d_c){
  //global indica que la funcion es para la GPU
  //Que tarea me toca

  int ix;
  ix=vlockIdx.x*blockDim.x+threadIdx.x;
  d_c[ix]=d_a[ix]+d_b[ix];
}

//-----Codigo del Host------

int main(void){
  int ix;
  //Declarar todas las variables por duplicado

  //----En el Host---
  float h_a[Lx],h_b[Lx],h_c[Lx];

  //----En el device----
  float *d_a; cudaMalloc((void**)&d_a, Lx*sizeof(float));
  float *d_b; cudaMalloc((void**)&d_b, Lx*sizeof(float));
  float *d_c; cudaMalloc((void**)&d_c, Lx*sizeof(float));

  //Llenar en el host los datos a procesar

  for(ix=0;ix<Lx;ix++){
    h_a[ix]=ix;
    h_b[ix]=2*x;
  }

  //Enviar datos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHosTtoDevice);
  cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHosTtoDevice);

  //Correr en Device
  dim3 ThreadsPerBlock(Nx,0,0);
  dim3 BlocksPerGrid(Mx,0,0);
  AddTwoVectors<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  //es necesario poner las flechas para correr una funcion en el device,
  //Entonces entre los signos se ponen las variables de block per grid y
  //threads per grid, y por ultimo los parentesis con los argumentos de la
  //funcion



  //DEvolver resultado al host

  cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDeviceToHost);

  //Imprimir resultados
  for(ix=0;ix<Lx;ix++)
    cout<<h_c[ix]<<endl;

  //Liberar la memoria dinamica
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_c);
return 0;

}
