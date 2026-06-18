////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License in the accompanying
// file LICENSE for more details
//
// centpack_main_2d_FD2.cc -- main function in
//
// Requires: initial_conditions.cc, spectral_radii.cc, evolution.cc, writeout.cc
// and those called by these (see README file for more information.
//
// This program approximates the solution of hyperbolic conservation laws in 
// two space dimensions (denoted by x and z),
//
//   u_t + f(u)_x + g(u)_y = 0                 (1)
//
// Using Jiang and Tadmor's 2nd order central scheme
//
// Input (read from file input):
// 
// (1) x_left -- a double type variable holding the left end point of the 
//               x-interval of the solution
//
// (2) x_right -- a double type variable holding the right end point of the 
//                x-interval of the solution
//
// (3) y_bottom -- a double type variable holding the left end point of the 
//               z-interval of the solution
//
// (4) y_top -- a double type variable holding the right end point of the 
//                z-interval of the solution
//
// (5) J  --  a long type variable holding the number of cells along x-dimension
//
// (6) K  --  a long type variable holding the number of cells along y-dimension
//
// (7) L  --  a long type variable holding the number of components in the 
//            system
//
// (8) gamma - a double type variable holding the ratio of specific heats
//
// (9) cfl  --  a double type variable holding the CFL restriction to determine 
//              time step
// 
// (10) t_final  --  a double type variable holding the time of simulation
//
// (11) dt_out  --  a double type variable holding the desired output interval
//
// (12) alpha  --  a double type variable holding the minmod3 parameter
//
// Output:
//
// This function does not produce any direct output except for monitoring 
// information (see README file).  The solutions variables are output via a 
// writeout function called at the time interval indicated by the user in the 
// input file
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"

using namespace std;

int CENTPACK::centpack_2d_FD2(int id, int p)
{
	if (id == 0)
		disclaimer();
	
	int tag;		
	MPI::Status status;
	
	double x_left, x_right, y_bottom, y_top, t_init, t_final;
	double dx, dy, dt, dtp, cfl, dt_out;
	double gamma, alpha;
	double sum_t = 0.0;
	double dt_cpu = 0.0;
	double t_start = clock();
	double t = 0.0;
	double t_out = 0.0;
	double pi = acos(-1.0);

	long J, K, L;
	long j, k, l;
	long n = 0;
	bool odd = false;
	
	doublearray1d parameters;
	
	read_in_2d(x_left, x_right, y_bottom, y_top, J, K, L, t_final, dt_out, cfl, alpha, parameters);
	
	doublearray1d xe(J+5), ye(K+5);
	doublearray1d xo(J+5), yo(K+5);
	doublearray1d xe_cell(J+4), ye_cell(K+4);
	doublearray1d xo_cell(J+4), yo_cell(K+4);
	doublearray1d dx_cell(J+4), dy_cell(K+4);
	doublearray1d dxe_cell(J+4), dye_cell(K+4);
	doublearray1d dxo_cell(J+4), dyo_cell(K+4);
	doublearray1d dxe_interface(J+3), dye_interface(K+3);
	doublearray1d dxo_interface(J+3), dyo_interface(K+3);
	doublearray1d lambda(J+4), mu(K+4);
	
	cout.setf(ios::scientific, ios::floatfield);
	
	t = 0.0;
	t_init = 0.0;
	
	doublearray3d un(J+4,K+4,L), unhalf(J+4,K+4,L);
	
	mesh(x_left, x_right, y_bottom, y_top, xe, xo, xe_cell, xo_cell, dxe_cell, dxo_cell, dxe_interface, dxo_interface, ye, yo, ye_cell, yo_cell, dye_cell, dyo_cell, dye_interface, dyo_interface, id, p);
	
	initial_conditions(un, parameters, dxe_cell, dye_cell, dxe_interface, dye_interface, xe_cell, ye_cell, xe, ye);
	
	write_mesh(xe_cell, xo_cell, ye_cell, yo_cell, id, p);
	writeout(un, t, parameters, n, odd, id, p);

	t_out = 0.0;
	n++;
	
	do
	{
		if (odd)
		{
			dx_cell = dxo_cell;
			dy_cell = dyo_cell;
		}
		
		else
		{
			dx_cell = dxe_cell;
			dy_cell = dye_cell;
		}
		
		time_step_2d(un, dx_cell, dy_cell, cfl, dtp, t, t_out, dt_out, lambda, mu, parameters);
		
		// MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Reduce(&dtp, &dt, 1, MPI::DOUBLE, MPI::MIN, 0);
		MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);
		
		dt = min(dt, dt_out - t_out);
		
		for (j = 2; j < J + 2; j++)
			lambda(j) = dt/dx_cell(j);
	
		for (k = 2; k < K + 2; k++)
			mu(k) = dt/dy_cell(k);
	
		t += dt;
		t_out += dt;
		
		reconstruction_2d_FD2(un, unhalf, alpha, odd);
		evolution_2d_FD2(un, unhalf, alpha, lambda, mu, parameters, odd, id, p);
		
		dt_cpu = (clock() - t_start)/CLOCKS_PER_SEC;
		sum_t = sum_t + dt_cpu;
		
		end_of_step_2d_FD2(un, dt, t, dt_out, t_out, n, sum_t, dt_cpu, parameters, odd, id, p);
		
		odd = !(odd);
		
		t_start = clock();
		
	}while(t < t_final);
	
	writeout(un, t, parameters, n, odd, id, p);
	run_info_2d(dt, sum_t, J, K, cfl, id, p);
	
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
