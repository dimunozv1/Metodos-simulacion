#include <iostream>
#include <cmath>
//#include "ConsoleApplication1.h"
using namespace std;
double funcion(double x)
{
	return sin(x);
}

double CerosBisecion(double a, double b) 
{
	double  m = 0, fa, fm;
	const double Epsilon = 1e-12;

	fa = funcion(a);

	while ((b - a) > Epsilon)
	{
		m = (a + b) / 2.0;
		fm = funcion(m);
		if (fa * fm > 0)
		{
			a = m;
			fa = funcion(a);
		}
		else
			b = m;
	}
	return (a + b) / 2;
}


double integralsin(double a, double b, int n)
{
	double x, h, suma;
	int i;
	n = n * 2;
	h = (b - a) / n;
	for (suma = 0, i = 0; i <= n; i += 1)
	{
		x = a + (i * h);
		if (i == 0 || i == n)
		{
			suma += funcion(x);

		}
		else if (i % 2 == 0)
		{
			suma += 2 * funcion(x);
		}
		else
		{
			suma += 4 * funcion(x);
		}
	}

	return suma * h / 3;
}

double funcion_interior(double x, double alpha, double t)
{
  //	return cos( (alpha * t) - (x * sin(t)) );
  return cos(t);
}

double integral_besel(double x, double alpha, int n,double a,double b)
{
	double t, h, sum;
	int i;
	n = n * 2;
	h = (b - a) / n;
	for (i = 0, sum = 0; i <= n; i += 1)
	{
		t = a + i * h;
		if (i == 0 || i == n)
		{
			sum += funcion_interior(x,alpha,t);
		}
		else if (i % 2 == 0)
		{
			sum += 2 * funcion_interior(x, alpha, t);
		}
		else
		{
			sum += 4 * funcion_interior(x, alpha, t);
		}
	}
	return sum * h / 3;

}

double CerosBisecion2(double a, double b,int alpha)
{
	double pi = CerosBisecion(a, b);
	double  m = 0, fa, fm;
	const double Epsilon = 1e-9;
	int n = 1000;

	fa = integral_besel(a, alpha, n, 0, pi) / pi;

	while ((b - a) > Epsilon)
	{
		m = (a + b) / 2.0;
		fm = integral_besel(m, alpha, n, 0, pi) / pi;
		if (fa * fm > 0)
		{
			a = m;
			fa = integral_besel(a, alpha, n, 0, pi) / pi;
		}
		else
			b = m;
		//cout << m << std::endl;

	}	
	return (a + b) / 2;

}
int main()
{
	double a = 2, b = 4,x,alpha;
	double pi = CerosBisecion(a, b);
	//double k = M_PI;
	alpha = 0;

	double ceroL = CerosBisecion2(2, 4,alpha);
	double o=integral_besel(0, 0, 50,0,M_PI/2);
	cout << "el cero esta en:" << o << std::endl;

}
