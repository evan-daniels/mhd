#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::initial_conditions(doublearray2d& un, const doublearray1d& parameters, const doublearray1d& x)
{
	
	long J = un.getIndex1Size() - 4;
	long L = un.getIndex2Size();
	long j, l;
	long j1 = J/2 + 2;
	double dx;
    double pi = acos(-1.0);
	
	for(l = 0; l < L; l++)
	{
		for (j = 0; j < J+4; j++)
			un(j,l) = sin(2.0*pi*x(j)) + 0.5*sin(pi*x(j));
	}
}
