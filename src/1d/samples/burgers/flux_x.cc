#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u_vector, const doublearray1d& parameters, doublearray1d& y)
{
	double u = u_vector(0);
	
	y(0) = 0.5*u*u;
}
