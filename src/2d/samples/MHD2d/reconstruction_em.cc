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
// reconstruction_em.cc
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

void CENTPACK::reconstruction_em(const doublearray3d& B, const doublearray2d& Bx, const doublearray2d& By, doublearray3d& B_N, doublearray3d& B_S, doublearray3d& B_E, doublearray3d& B_W, doublearray3d& Jc_N, doublearray3d& Jc_S, doublearray3d& Jc_E, doublearray3d& Jc_W, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface,const double& alpha)
{
	long J = B.getIndex1Size() - 4;
	long K = B.getIndex2Size() - 4;
	long j, k, l;
	
	double ux, uy;

	// reconstruct current J from its cell averages
	
	for (l = 0; l < 3; l++)
	{
		for (k = 1; k < K+3; k++)
		{
			for (j = 1; j < J+3; j++)
			{
				ux = minmod3(alpha*(Jc(j+1,k,l) - Jc(j,k,l))/dx_interface(j), (Jc(j+1,k,l) - Jc(j-1,k,l))/(dx_interface(j-1) + dx_interface(j)), alpha*(Jc(j,k,l) - Jc(j-1,k,l))/dx_interface(j-1));
				
				uy=minmod3(alpha*(Jc(j,k+1,l) - Jc(j,k,l))/dy_interface(k), (Jc(j,k+1,l) - Jc(j,k-1,l))/(dy_interface(k-1) + dy_interface(k)), alpha*(Jc(j,k,l) - Jc(j,k-1,l))/dy_interface(k-1));
			
				Jc_N(j,k,l) = Jc(j,k,l) + 0.5*dy_cell(k)*uy;
				Jc_S(j,k,l) = Jc(j,k,l) - 0.5*dy_cell(k)*uy;
				Jc_E(j,k,l) = Jc(j,k,l) + 0.5*dx_cell(j)*ux;
				Jc_W(j,k,l) = Jc(j,k,l) - 0.5*dx_cell(j)*ux;
			}
		}
	}

	// Bx reconstruction

	for (j = 1; j < J+3; j++)
	{
		for (k = 1; k < K+3; k++)
		{
			uy = minmod3(alpha*(B(j,k+1,0) - B(j,k,0))/dy_interface(k), (B(j,k+1,0) - B(j,k-1,0))/(dy_interface(k-1) + dy_interface(k)), alpha*(B(j,k,0) - B(j,k-1,0))/dy_interface(k-1));
				
			B_N(j,k,0) = B(j,k,0) + 0.5*dy_cell(k)*uy;
			B_S(j,k,0) = B(j,k,0) - 0.5*dy_cell(k)*uy;
			B_E(j,k,0) = Bx(j+1,k);
			B_W(j,k,0) = Bx(j,k);
		}
	}

	// By reconstruction

	for (k = 1; k < K+3; k++)
	{
		for (j = 1; j < J+3; j++)
		{
			ux = minmod3(alpha*(B(j+1,k,1) - B(j,k,1))/dx_interface(j), (B(j+1,k,1) - B(j-1,k,1))/(dx_interface(j-1) + dx_interface(j)), alpha*(B(j,k,1) - B(j-1,k,1))/dx_interface(j-1));

			B_N(j,k,1) = By(j,k+1);
			B_S(j,k,1) = By(j,k);
			B_E(j,k,1) = B(j,k,1) + 0.5*dx_cell(j)*ux;
			B_W(j,k,1) = B(j,k,1) - 0.5*dx_cell(j)*ux;
		}
	}

	// Bz reconstruction

	for (k = 1; k < K+3; k++)
	{
		for (j = 1; j < J+3; j++)
		{
			ux = minmod3(alpha*(B(j+1,k,2) - B(j,k,2))/dx_interface(j), (B(j+1,k,2) - B(j-1,k,2))/(dx_interface(j-1) + dx_interface(j)), alpha*(B(j,k,2) - B(j-1,k,2))/dx_interface(j-1));

			uy=minmod3(alpha*(B(j,k+1,2) - B(j,k,2))/dy_interface(k), (B(j,k+1,2) - B(j,k-1,2))/(dy_interface(k-1) + dy_interface(k)), alpha*(B(j,k,2) - B(j,k-1,2))/dy_interface(k-1));
				
			B_N(j,k,2) = B(j,k,2) + 0.5*dy_cell(k)*uy;
			B_S(j,k,2) = B(j,k,2) - 0.5*dy_cell(k)*uy;
			B_E(j,k,2) = B(j,k,2) + 0.5*dx_cell(j)*ux;
			B_W(j,k,2) = B(j,k,2) - 0.5*dx_cell(j)*ux;
		}
	}
}
