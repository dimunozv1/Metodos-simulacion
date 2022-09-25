#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

int main(void){

  vector3D a,b,c;

  a.load(1,2,2);
  c.load(2,1,1);
  double k=a.operator*(c);
  b.show();
  cout<<k<<endl;


  return 0;
}
