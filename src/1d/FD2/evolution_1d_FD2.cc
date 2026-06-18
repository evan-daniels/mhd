#include "centpack_1d_FD2.h"

using namespace std;

void CENTPACK::evolution_1d_FD2(doublearray2d& un, const doublearray2d& uj_half, const doublearray1d& lambda, const double& alpha, const doublearray1d& parameters, const bool& odd, const int& id, const int& p)
{

	long j, l, J, L;
	
	J = un.getIndex1Size() - 4;
	L = un.getIndex2Size();
	
	doublearray2d un_half(J+4,L);
	
	predictor_1d_FD2(un, un_half, alpha, lambda, parameters);
	corrector_1d_FD2(un, uj_half, un_half, lambda, parameters, odd);
	
	boundary_conditions(un, odd, id, p);
}
