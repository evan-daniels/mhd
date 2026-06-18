#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::time_step_1d(const doublearray2d& un, const doublearray1d& dx_cell, const double& cfl, double& dtp, double& t, double& t_out, double& dt_out, doublearray1d& lambda, const doublearray1d& parameters)
{
	long J = un.getIndex1Size() - 4;
	long L = un.getIndex2Size();
	long j, l;
	double rx, r_maxx;
	double dx_min = dx_cell(2);
	
	doublearray1d u_vector(L);
	
	r_maxx = 0.0;
	
	for (j = 2; j < J+2; j++)
		dx_min = min(dx_min, dx_cell(j));
	
	for (j = 2; j < J+2; j++)
	{
		for (l = 0; l < L; l++)
			u_vector(l) = un(j,l);
		
		spectral_radius(u_vector, parameters, rx);
		
		if (rx > r_maxx)
			r_maxx = rx;
	}
	
	dtp = dx_min*cfl/r_maxx;
	dtp = min(dtp, dt_out - t_out);
}
