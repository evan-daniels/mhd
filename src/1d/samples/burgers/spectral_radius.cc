#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::spectral_radius(const doublearray1d& u, const doublearray1d& parameters, double& rx)
{
	rx = fabs(u(0));
}
