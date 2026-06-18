#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u, const doublearray1d& parameters, doublearray1d& y)
{
	double rho = u(0);
	double u1 = u(1)/rho;
	double E = u(2);
	double gamma = parameters(0);
	double p = (gamma - 1.0)*(E - 0.5*rho*pow(u1,2.0));
	
	y(0) = rho*u1;
	y(1) = rho*pow(u1,2.0) + p;
	y(2) = u1*(E + p);
}




















