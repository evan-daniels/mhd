#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u, const doublearray1d& parameters, doublearray1d& y)
{
	double p, p_star;
	double gamma = parameters(0);
	double B1 = parameters(1);
	
	p = u(6) - 0.5*((pow(u(1),2.0) + pow(u(2),2.0) + pow(u(3),2.0))/u(0)) - 0.5*(pow(B1,2.0) + pow(u(4),2.0) + pow(u(5),2.0));
	p_star = p + 0.5*(pow(B1,2.0) + pow(u(4),2.0) + pow(u(5),2.0));
	
	y(0) = u(1);
	y(1) = (pow(u(1),2.0)/u(0)) + p_star;
	y(2) = (u(1)*u(2)/u(0)) - B1*u(4);
	y(3) = (u(1)*u(3)/u(0)) - B1*u(5);
	y(4) = (u(4)*u(1)/u(0)) - B1*u(2)/u(0); // change this to canonical vorticity
	y(5) = (u(5)*u(1)/u(0)) - B1*u(3)/u(0); // change this to canonical vorticity 
	y(6) = (u(6) + p_star)*(u(1)/u(0)) - B1*(B1*(u(1)/u(0)) + (u(2)/u(0))*u(4) + (u(3)/u(0))*u(5));
}




















