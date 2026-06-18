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
// reconstruction.cc -- a core function of CentPack
//
// Function called by evolution_2d_SD2
//
// Requires: minmod3
//
// This function implements a 3rd order Kurganov-Levy CWENO reconstruction
// of the interface values of the solution u from its cell averages along
// the x- and z-directions
// 
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t=t^n over the discretized solution domain
//
// Output (returned by reference):
//
// (1) Four doublearray3d type variables
//
//     u_N  --  holds the point values of u at the cell interfaces
//              (x_j, y_k+1/2) over the entire solution domain
//
//     u_S  --  holds the point values of u at the cell interfaces
//              (x_j, y_k-1/2) over the entire solution domain
//
//     u_E  --  holds the point values of u at the cell interfaces
//              (x_j+1/2, y_k) over the entire solution domain
//
//     u_W  --  holds the point values of u at the cell interfaces
//              (x_j-1/2, y_k) over the entire solution domain
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::reconstruction_2d_SD2(doublearray3d& un, doublearray3d& u_N, doublearray3d& u_S, doublearray3d& u_E, doublearray3d& u_W, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface,const double& alpha)
{

	long J = un.getIndex1Size() - 4;
	long K = un.getIndex2Size() - 4;
	long L = un.getIndex3Size();
	long j, k, l;
	
	double ux, uy;
	
	for (l = 0; l < L; l++)
	{		
		for (k = 1; k < K+3; k++)
		{
			for (j = 1; j < J+3; j++)
			{
				ux = minmod3(alpha*(un(j+1,k,l) - un(j,k,l))/dx_interface(j), (un(j+1,k,l) - un(j-1,k,l))/(dx_interface(j-1) + dx_interface(j)), alpha*(un(j,k,l) - un(j-1,k,l))/dx_interface(j-1));
				uy=minmod3(alpha*(un(j,k+1,l) - un(j,k,l))/dy_interface(k), (un(j,k+1,l) - un(j,k-1,l))/(dy_interface(k-1) + dy_interface(k)), alpha*(un(j,k,l) - un(j,k-1,l))/dy_interface(k-1));
				
				u_N(j,k,l) = un(j,k,l) + 0.5*dy_cell(k)*uy;
				u_S(j,k,l) = un(j,k,l) - 0.5*dy_cell(k)*uy;
				u_E(j,k,l) = un(j,k,l) + 0.5*dx_cell(j)*ux;
				u_W(j,k,l) = un(j,k,l) - 0.5*dx_cell(j)*ux;
			}
		}
	}
}
