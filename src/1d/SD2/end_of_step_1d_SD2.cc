#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::end_of_step_1d_SD2(const doublearray2d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& sum_t, const double& dt_cpu, const doublearray1d& parameters, const int& id, const int& p)
{
	if (t_out == dt_out)
	{
		writeout(un, t, parameters, n, id, p);
		n++;
		
		t_out = 0.0;
		
		cout<<"output written at t = "<< t <<endl;
	}
	
	if (id == 0)
	{
		cout<<"run = "<< t <<",dt = "<< dt <<",cpu_t = "<< dt_cpu <<",t = "<< sum_t <<" in processor "<< id << endl;
	}
}
