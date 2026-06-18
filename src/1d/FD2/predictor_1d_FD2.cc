#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::predictor_1d_FD2(const doublearray2d& un, doublearray2d& un_half, const double& alpha, const doublearray1d& lambda, const doublearray1d& parameters)
{
	long j, l, J, L;
	
	J = un.getIndex1Size() - 4;
	L = un.getIndex2Size();
	
	doublearray1d u_vector(L);
	doublearray1d f_vector(L);
	doublearray2d f(J+4,L);
	doublearray2d f_prime(J+4,L);
	
	for (j = 0; j < J+4; j++)
	{
		for (l = 0; l < L; l++)
			u_vector(l) = un(j,l);
		
		flux_x(u_vector, parameters, f_vector);
		
		for (l = 0; l < L; l++)
			f(j,l) = f_vector(l);
	}
	
	for (j = 1; j < J+3; j++)
	{
		for (l = 0; l < L; l++)
			f_prime(j,l) = minmod3(alpha*(f(j,l)-f(j-1,l)), .5*(f(j+1,l)-f(j-1,l)), alpha*(f(j+1,l)-f(j,l)));
	}
	
	for (j = 1; j < J+3; j++)
	{
		for (l = 0; l < L; l++)
			un_half(j,l) = un(j,l) - 0.5*lambda(j)*f_prime(j,l);
	}
}
