#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::spectral_radius(const doublearray1d& u, const doublearray1d& parameters, double& rx)
{

	double rho, ux, p, A, c;
	double gamma = parameters(0);
	
	rho = u(0);
	ux = u(1)/rho;
	p = (gamma -1.0)*(u(2) - .5*rho*pow(ux,2.0));
	A = gamma*p/rho;
	c = sqrt(A);
	rx = fabs(ux) + c;
}
