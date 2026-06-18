#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::H_flux_1d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, const doublearray1d& parameters, doublearray1d& H)
{
	long L = u_e.getIndex1Size();
	long l;
	double re, rw, a;
	
	doublearray1d f_e(L), f_w(L);
	
	spectral_radius(u_e, parameters, re);
	spectral_radius(u_w, parameters, rw);

	if (re > rw)
		a = re;
	else
		a = rw;
	
  	flux_x(u_e, parameters, f_e);
	flux_x(u_w, parameters, f_w);
	
	for (l = 0; l < L; l++)
		H(l) = 0.5*(f_w(l) + f_e(l)) - 0.5*a*(u_w(l) - u_e(l));
}
