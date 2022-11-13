#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;


int main(void){

  Crandom ran64(1);
  for(int i=0;i<5;i++){
  cout<<ran64.r()<<endl;
    }


  return 0;
}
