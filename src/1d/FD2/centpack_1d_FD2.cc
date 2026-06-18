////////////////////////////////////////////////////////////////////////////
// 2nd order fully-discrete Central Scheme for 1-d hyperbolic conservation laws:
// minmod3 reconstruction + midpoint rule for time integration of fluxes
////////////////////////////////////////////////////////////////////////////

#include "centpack_1d_FD2.h"

using namespace std;

int CENTPACK::centpack_1d_FD2(int id, int p)
{
	if (id == 0)
		disclaimer();
	
	int tag;		
	MPI::Status status;
	
	double dx, x_left, x_right, t, t_init, t_final, t_out;
	double dt, cfl, dt_out, alpha;
	double dtp;
	double sum_t = 0.0;
	double dt_cpu = 0.0;
	double t_start = clock();
	
	long J, L;
	long j, l;
	long n = 0;
	bool odd = false;
	
	int i;
	
	doublearray1d parameters;
	
	read_in_1d(x_left, x_right, J, L, t_final, dt_out, cfl, alpha, parameters);
	
	doublearray1d xe(J+5);
	doublearray1d xo(J+5);
	doublearray1d xe_cell(J+4);
	doublearray1d xo_cell(J+4);
	doublearray1d dx_cell(J+4);
	doublearray1d dxe_cell(J+4);
	doublearray1d dxo_cell(J+4);
	doublearray1d dxe_interface(J+3);
	doublearray1d dxo_interface(J+3);
	doublearray1d lambda(J+4);	
	//cout.setf(ios::scientific, ios::floatfield);
	
	t = 0.0;
	t_init = 0.0;
	
	doublearray2d un(J+4,L), uj_half(J+4,L);

	mesh(x_left, x_right, xe, xo, xe_cell, xo_cell, dxe_cell, dxo_cell, dxe_interface, dxo_interface, id, p);
	
	initial_conditions(un, parameters, xe_cell);

	write_mesh(xe_cell, xo_cell, id, p);
	writeout(un, t, parameters, n, odd, id, p);
	
	t_out = 0.0;
	n++;
	
	do
	{
		if (odd)
			dx_cell = dxo_cell;
		else
			dx_cell = dxe_cell;
		
		time_step_1d(un, dx_cell, cfl, dtp, t, t_out, dt_out, lambda, parameters);
		
		// MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Reduce(&dtp, &dt, 1, MPI::DOUBLE, MPI::MIN, 0);
		MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);
				
		t += dt;
		t_out += dt;
			
		for (j = 2; j < J+2; j++)
			lambda(j) = dt/dx_cell(j);
				
		// MPI::COMM_WORLD.Barrier();
		
		reconstruction_1d_FD2(un, uj_half, alpha, odd);
		evolution_1d_FD2(un, uj_half, lambda, alpha, parameters, odd, id, p);
		
		dt_cpu = (clock() - t_start)/CLOCKS_PER_SEC;
		sum_t += dt_cpu;
		
		end_of_step_1d_FD2(un, dt, t, dt_out, t_out, n, sum_t, dt_cpu, parameters, odd, id, p);
		
		odd = !(odd);
		
		t_start = clock();
		
	}while(t < t_final);
	
	writeout(un, t, parameters, n, odd, id, p);
	run_info_1d(dt, sum_t, J, cfl, id, p);
	
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
