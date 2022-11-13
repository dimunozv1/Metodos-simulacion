#include <iostream>
#include <cmath>
using namespace std;
const double w = 1,alpha=0;


double f1(double r, double x1, double x2, double dt)
{	
	//cout << -(x1 / r) - lamda * lamda * x2 << "  return f1  " << endl;
	return  - w * w * x2;

}
double f2(double r, double x1, double x2, double dt)
{	
	//cout << x1 << "  return f2  " << endl;
	return x1;
}

void RK4_unpaso(double& r, double& x1,double& x2, double dr)
{	
  //cout << r << " " << x1 << " " << x2 << " " << dr;
	double k1p, k2p, k3p, k4p, k1v, k2v, k3v, k4v,dp,dv;
	k1p = dr * f1(r, x1, x2, dr);						k1v = dr * f2(r, x1, x2, dr);
	k2p = dr * f1(r + dr / 2, x1 + k1p / 2, x2 + k1v / 2, dr);		k2v = dr * f2(r + dr / 2, x1 + k1p / 2, x2 + k1v / 2, dr);
	k3p = dr * f1(r + dr / 2, x1 + k2p / 2, x2 + k2v / 2, dr);		k3v = dr * f2(r + dr / 2, x1 + k2p / 2, x2 + k2v / 2, dr);
	k4p = dr * f1(r + dr, x1 + k3p, x2 + k3v,dr);				k4v = dr * f2(r + dr, x1 + k3p, x2 + k3v, dr);
	
	dp = (k1p + 2 * (k2p + k3p) + k4p) / 6;					dv = (k1v + 2 * (k2v + k3v) + k4v) / 6;
	//	cout << k4p << "\t" << k4v << endl;
	r += dr; x1 += dp;	x2 += dv;


}

int main()
{	
	double r = 0,dr=0.1,x1=0,x2=1;
	for (; r < 2;)
	{
	   std::cout << r << " " << x2 <<"\t"<<cos(r)<< endl;
	  	RK4_unpaso(r, x1, x2, dr);

	}
	return 0;
}
