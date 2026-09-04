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
// time_step_2d.cc -- a core function of CentPack
//
// Function called by centpack_2d_SD2() (main)
//
// Requires: spectral_radii
//
// This function calculates size of the time step to be used in the evolution of
// the solution variables at each time step
//
// Input (passed by reference):
// 
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n over the discretized solution domain
// (2) dx  --  a double type variable holding the space scale dx
//
// (3) dy  --  a double type variable holding the space scale dy
//
// (4) cfl  --  a double type variable holding the CFL restriction for the 
//              central scheme implemented
//
// (5) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) dt  --  a double type variable holding the value of the next time step to 
//             be used in next evolution step
//
// (2) lambda  --  a double type variable holding themesh ratio dt/dx
//
// (3) mu  --  a double type variable holding the mesh ratio dt/dy
//
// (4) t  --  a double type variable holding the simulation time completed so 
//            far
//
// (5) t_out  --  a double type variable holding the simulation time since last 
//                output
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

time_step_em(un, B, Jc, dx_cell, dy_cell, cfl_em, dtp_em, dt_fluid, parameters);

void CENTPACK::time_step_em(const doublearray3d& un, const doublearray3d& B, const doublearray3d& Jc, const doublearray1d& dx_cell, const doublearray1d& dy_cell, const double& cfl, double& dtp, double& t, double& t_out, double& dt_out, const doublearray1d& parameters)
{
	long J = B.getIndex1Size() - 4;
	long K = B.getIndex2Size() - 4;
	//long L = un.getIndex3Size();
	long j, k, l;

	double rx, ry, r_maxx, r_maxy;
	double dxj, dyk;
	double dx_min = dx_cell(2);
	double dy_min = dy_cell(2);
	
	doublearray1d u_vector(L);
	
	r_maxx = 0.0;
	r_maxy = 0.0;
	
	for (j = 2; j < J + 2; j++)
		dx_min = min(dx_min, dx_cell(j));
	
	for (k = 2; k < K + 2; k++)
		dy_min = min(dy_min, dy_cell(k));
	
	for (j = 2; j < J + 2; j++)
	{
		dxj = dx_cell(j);
		for (k = 2; k < K + 2; k++)
		{
			dyk = dy_cell(k);

			for (l = 0; l < L; l++)
				u_vector(l) = un(j,k,l);
		
			spectral_radii_em(u_vector, dxj, dyk, parameters, rx, ry);
		
			if (rx > r_maxx)
				r_maxx = rx;
			if (ry > r_maxy)
				r_maxy = ry;
		}
	}
	
	dtp_em = cfl_em/sqrt((r_maxx/dx_min)*(r_maxx/dx_min) + (r_maxy/dy_min)*(r_maxy/dy_min));
}