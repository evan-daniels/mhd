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
// evolution_2d_SD2.cc -- a core function of CentPack
//
// Function called by centpack_2d_SD2 (main)
//
// Requires: boundary_conditions, norm, reconstruction, and C_flux (see 
// prototypes below)
//
// This function evolves the cell averages of u from t=t^n to t=t^n + dt
// 
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n over the discretized solution domain
//
// (2) dx  --  a double type variable holding the space scale dx
//
// (3) dy  --  a double type variable holding the space scale dy
//
// (4) lambda  --  a double type variable holding the mesh ratio dt/dx
//
// (5) mu  --  a double type variable holding the mesh ratio dt/dy
//
// (6) gamma  --  a double type variable holding the ratio of specific heats
//
// (7) alpha  --  minmod limiting parameter
//
// Output (returned by reference):
//
// (1)  un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n + dt over the discretized solution domain
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

void CENTPACK::evolution_2d_SD2(doublearray3d& un, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha, const doublearray1d& parameters, const doublearray3d& Efield, const int& id, const int& p)
{

	using namespace std;

	long j, k, l, J, K, L;
	
	J = un.getIndex1Size() - 4;
	K = un.getIndex2Size() - 4;
	L = un.getIndex3Size();
	
	doublearray3d u_N(J+4,K+4,L);
	doublearray3d u_S(J+4,K+4,L);
	doublearray3d u_E(J+4,K+4,L);
	doublearray3d u_W(J+4,K+4,L);
	doublearray3d C0(J+4,K+4,L), C1(J+4,K+4,L);
	
	//boundary_conditions(un);


	reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
	C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Efield, C0);
	
	for (l = 0; l < L; l++)
	{
		for (j = 2; j < J + 2; j++)
		{
			for (k = 2; k < K + 2; k++)
				un(j,k,l) += C0(j,k,l);
		}
	}
	
	boundary_conditions(un, parameters, id, p);	
	reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
	C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Efield, C1);
	
	for (l = 0; l < L; l++)
	{
		for (j = 2; j < J + 2; j++)
		{
			for (k = 2; k < K + 2; k++)
				un(j,k,l) += 0.5*(C1(j,k,l) - C0(j,k,l));
		}
	}
	
	boundary_conditions(un, parameters, id, p);
}
