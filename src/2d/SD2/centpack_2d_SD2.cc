////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2005 Jorge Balbas and Eitan Tadmor
//
// This program is distributed freely for research and instructional use only. 
// You may copy, modify,  and use this software for these non-commercial
// purposes, provided that the copyright notice and associated text is 
// reproduced on all copies
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.
//
// centpack_2d_SD2.cc -- main function of CentPack
//
// Requires: read_in, initial_conditions.cc, evolution.cc, writeout.cc
// and those called by these (see README file for more information)f
//
// This program approximates the solution of hyperbolic conservation laws in 
// two space dimensions (denoted by x and y),
//
//   u_t + f(u)_x + g(u)_y = 0                 (1)
//
// Input (read from file input):
// 
// (1) x_init -- a double type variable holding the left end point of the 
//               x-interval of the solution
//
// (2) x_final -- a double type variable holding the right end point of the 
//                x-interval of the solution
//
// (3) y_init -- a double type variable holding the left end point of the 
//               y-interval of the solution
//
// (4) y_final -- a double type variable holding the right end point of the 
//                y-interval of the solution
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
// (12) alpha  --  minmod limiting parameter
//
// Output:
//
// This function does not produce any direct output except for monitoring 
// information (see README file).  The solutions variables are output via a 
// writeout function called at the time interval indicated by the user in the 
// input file
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

int CENTPACK::centpack_2d_SD2(int id, int p)
{
	if (id == 0)
		disclaimer();
	
	int tag;		
	MPI::Status status;

	double x_left, x_right, y_bottom, y_top, t, t_init, t_final, t_out;
	double dx, dy, dtp, dt, cfl, dt_out, alpha;
	double gamma;
	double sum_t = 0.0;
	double dt_cpu = 0.0;
	double t_start = clock();
	double pi = acos(-1.0);
	long J, K, L;
	long j, k, l;
	long n = 0;
	
	doublearray1d parameters;
	
	read_in_2d(x_left, x_right, y_bottom, y_top, J, K, L, t_final, dt_out, cfl, alpha, parameters);
	
	doublearray1d x(J+5), y(K+5);
	doublearray1d x_cell(J+4), y_cell(K+4);
	doublearray1d dx_cell(J+4), dy_cell(K+4);
	doublearray1d dx_interface(J+3), dy_interface(K+3);
	doublearray1d lambda(J+4), mu(K+4);
	
	cout.setf(ios::scientific, ios::floatfield);
	
	t = 0.0;
	t_init = 0.0;
	
	doublearray3d un(J+4,K+4,L);
	doublearray3d B1_xf(J+4,K+4,1), B2_yf(J+4,K+4,1);
	
	mesh(x_left, x_right, y_bottom, y_top, x, x_cell, dx_cell, dx_interface, y, y_cell, dy_cell, dy_interface, id, p);
	
	write_mesh(x_cell, y_cell, id, p);
	
	initial_conditions(un, parameters, dx_cell, dy_cell, dx_interface, dy_interface, x_cell, y_cell, x, y);

	for (k = 2; k < K+2; k++)
		for (j = 2; j < J+2; j++)
		{
			B1_xf(j,k,0) = 0.5*(un(j,k,4) + un(j+1,k,4));
			B2_yf(j,k,0) = 0.5*(un(j,k,5) + un(j,k+1,5));
		}

	boundary_conditions_ct(B1_xf, id, p);
	boundary_conditions_ct(B2_yf, id, p);
		
	writeout(un, t, parameters, n, id, p);

	t_out = 0.0;
	n++;


	do
	{
		time_step_2d(un, dx_cell, dy_cell, cfl, dtp, t, t_out, dt_out, lambda, mu, parameters);

		if (std::isnan(dtp) || std::isinf(dtp)) {
			MPI::COMM_WORLD.Abort(1);
		}

		MPI::COMM_WORLD.Barrier(); // COMMENT OR UNCOMMENT BASED ON WHAT DR BALBAS RECOMMENDS
		MPI::COMM_WORLD.Reduce(&dtp, &dt, 1, MPI::DOUBLE, MPI::MIN, 0);
		MPI::COMM_WORLD.Bcast(&dt, 1, MPI::DOUBLE, 0);
		
		dt = min(dt, dt_out - t_out);
		
		for (j = 2; j < J + 2; j++)
			lambda(j) = dt/dx_cell(j);
	
		for (k = 2; k < K + 2; k++)
			mu(k) = dt/dy_cell(k);

		t += dt;
		t_out += dt;
		
		evolution_2d_SD2(un, B1_xf, B2_yf, lambda, mu, dx_cell, dx_interface, dy_cell, dy_interface, alpha, parameters, id, p);
		
		double max_divB = 0.0;
		long max_j = -1, max_k = -1;
		for (k = 2; k < K+2; k++)
			for (j = 2; j < J+2; j++)
			{
				double divB = (B1_xf(j,k,0) - B1_xf(j-1,k,0)) / dx_cell(j)
							+ (B2_yf(j,k,0) - B2_yf(j,k-1,0)) / dy_cell(k);
				if (std::fabs(divB) > max_divB) { max_divB = std::fabs(divB); max_j = j; max_k = k; }
			}
		if (id == 0) printf("t=%.4f  max|div B| = %.3e  at (j,k)=(%ld,%ld)\n", t, max_divB, max_j, max_k);

		if (max_j > 0 && max_k > 0)
    	printf("  at max-divB cell: rho=%.3e  Bx_face=%.3e  By_face=%.3e\n",
           un(max_j,max_k,0), B1_xf(max_j,max_k,0), B2_yf(max_j,max_k,0));
		
		if (id == 0) printf("t=%.4f  dt=%.3e  max|div B| = %.3e  at (j,k)=(%ld,%ld)\n", t, dt, max_divB, max_j, max_k);

		// resistivity_step(un, dx_cell, dy_cell, dt, parameters);
		// hall_step(un, dx_cell, dy_cell, dt, parameters);
		// int nsub = 50;
		// for (int isub = 0; isub < nsub; isub++)
		// 	hall_step(un, dx_cell, dy_cell, dt/nsub, parameters);
		
		dt_cpu = (clock() - t_start)/CLOCKS_PER_SEC;
		sum_t = sum_t + dt_cpu;
		
		end_of_step_2d_SD2(un, dt, t, dt_out, t_out, n, sum_t, dt_cpu, parameters, id, p);
				
		t_start = clock();
		
	}while(t < t_final);
	
	writeout(un, t, parameters, n, id, p);
	run_info_2d(dt, sum_t, J, K, cfl, id, p);
	
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
