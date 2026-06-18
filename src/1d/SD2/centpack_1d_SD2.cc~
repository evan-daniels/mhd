////////////////////////////////////////////////////////////////////////////
// 2nd order semi-discrete Central Scheme for 1-d MHD Equations:
// minmod3 reconstruction + 2nd order SSP-RK
////////////////////////////////////////////////////////////////////////////

#include "centpack_1d_SD2.h"

using namespace std;

int CENTPACK::centpack_1d_SD2(int id, int p)
{
	if (id == 0)
		disclaimer();
	
	int tag;		
	MPI::Status status;
		
	double dx, x_left, x_right, t, t_init, t_final, t_out;
	double B1;
	double dt, cfl, dt_out, alpha;
	double dtp;
	double gamma;
	double sum_t = 0.0;
	double dt_cpu = 0.0;
	double t_start = clock();
	
	long J, L;
	long j, l;
	long n = 0;

	int i;

	doublearray1d parameters;
		
	read_in_1d(x_left, x_right, J, L, t_final, dt_out, cfl, alpha, parameters);
	
	doublearray1d x(J+5);
	doublearray1d x_cell(J+4);
	doublearray1d dx_cell(J+4);
	doublearray1d dx_interface(J+3);
	doublearray1d lambda(J+4);
	
	cout.setf(ios::scientific, ios::floatfield);
	
	t = 0.0;
	t_init = 0.0;
	
	doublearray2d un(J+4,L);
	
	mesh(x_left, x_right, x, x_cell, dx_cell, dx_interface, id, p);
	
	initial_conditions(un, parameters, x_cell);
	
	write_mesh(x_cell, id, p);
	writeout(un, t, parameters, n, id, p);
	
	t_out = 0.0;
	n++;
	
	do
	{
		time_step_1d(un, dx_cell, cfl, dtp, t, t_out, dt_out, lambda, parameters);
		
		// MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Reduce(&dtp, &dt, 1, MPI::DOUBLE, MPI::MIN, 0);
		MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);
				
		t += dt;
		t_out += dt;
			
		for (j = 2; j < J+2; j++)
			lambda(j) = dt/dx_cell(j);
				
		// MPI::COMM_WORLD.Barrier();
		
		evolution_1d_SD2(un, lambda, dx_cell, dx_interface, alpha, parameters, id, p);
		
		dt_cpu = (clock() - t_start)/CLOCKS_PER_SEC;
		sum_t = sum_t + dt_cpu;
		
		end_of_step_1d_SD2(un, dt, t, dt_out, t_out, n, sum_t, dt_cpu, parameters, id, p);
		
		t_start = clock();
			
	}while(t < t_final);
	
	writeout(un, t, parameters, n, id, p);
	run_info_1d(dt, sum_t, J, cfl, id, p);
	
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
