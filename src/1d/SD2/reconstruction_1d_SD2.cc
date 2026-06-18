////////////////////////////////////////////////////////////////////////////
// minmod polynomial reconstruction of interface values
////////////////////////////////////////////////////////////////////////////

#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::reconstruction_1d_SD2(doublearray2d& un, doublearray2d& u_E, doublearray2d& u_W, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const double& alpha)
{
	long J = un.getIndex1Size() - 4;
	long L = un.getIndex2Size();
	long j, l;

	double u, s;

	for (l = 0; l < L; l++)
	{
		for (j = 1; j < J+3; j++)
		{
			u = un(j,l);
			s = minmod3(alpha*(un(j+1,l) - un(j,l))/dx_interface(j), (un(j+1,l) - un(j-1,l))/(dx_interface(j) + dx_interface(j-1)), alpha*(un(j,l) - un(j-1,l))/dx_interface(j-1));
		
			u_E(j,l) = u + .5*dx_cell(j)*s;
			u_W(j,l) = u - .5*dx_cell(j)*s;
		}
	}
}
