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
	double p_left, p_right;
	
	double gamma = parameters(0);
	
	doublearray1d u_left(L), u_right(L);
	
	p_left = 1.0;
	p_right = 0.1;
	
	u_left(0) = 1.0;
	u_left(1) = 0.0;
	u_left(2) = p_left/(gamma - 1.0);
	
	u_right(0) = 0.125; 
	u_right(1) = 0.0;
	u_right(2) = p_right/(gamma - 1.0);
	
	for(l = 0; l < L; l++)
	{
	  for (j = 0; j < J+4; j++)
	  {
		  if (x(j) < 0.0)
			  un(j,l) = u_left(l);
		  else
			  un(j,l) = u_right(l);
	  }
	}
}
