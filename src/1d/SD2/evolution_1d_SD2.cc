#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::evolution_1d_SD2(doublearray2d& un, const doublearray1d& lambda, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const double& alpha, const doublearray1d& parameters, const int& id, const int& p)
{
	long j, l, J, L;
	
	MPI::Status status;
	int tag;
	int index;
		
	J = un.getIndex1Size() - 4;
	L = un.getIndex2Size();
	
	doublearray2d u_E(J+4,L);
	doublearray2d u_W(J+4,L);
	doublearray2d C0(J+4,L), C1(J+4,L);
	
	reconstruction_1d_SD2(un, u_E, u_W, dx_cell, dx_interface, alpha);
	C_flux_1d_SD2(u_E, u_W, lambda, parameters, C0);

	for (l = 0; l < L; l++)
	{
		for (j = 2; j < J+2; j++)
			un(j,l) = un(j,l) + C0(j,l);
	}
	
	boundary_conditions(un, id, p);
	
	reconstruction_1d_SD2(un, u_E, u_W, dx_cell, dx_interface, alpha);
	C_flux_1d_SD2(u_E, u_W, lambda, parameters, C1);
	
	for (l = 0; l < L; l++)
	{
		for (j = 2; j < J + 2; j++)
			un(j,l) = un(j,l) + 0.5*(C1(j,l) - C0(j,l));
	}
	
	boundary_conditions(un, id, p);
}
