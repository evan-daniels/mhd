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

void CENTPACK::electric_field(const doublearray2d& rho_N, const doublearray2d& rho_S, const doublearray2d& rho_E, const doublearray2d& rho_W,const doublearray3d& v_N, const doublearray3d& v_S, const doublearray3d& v_E, const doublearray3d& v_W, const doublearray3d& B_N, const doublearray3d& B_S, const doublearray3d& B_E, const doublearray3d& B_W, const doublearray3d& J_N, const doublearray3d& J_S, const doublearray3d& J_E, const doublearray3d& J_W, doublearray3d& E_N, doublearray3d& E_S, doublearray3d& E_E, doublearray3d& E_W, const double& di, const double& eta)
{
	long J = B_N.getIndex1Size() - 4;
	long K = B_N.getIndex2Size() - 4;
	long j, k, l;
	
	double ux, uy;
	
	for (k = 1; k < K+3; k++)
	{
		for (j = 1; j < J+3; j++)
		{
			rho = rho_N(j,k);
			
			vx = v_N(j,k,0);
			vy = v_N(j,k,1);
			vz = v_N(j,k,2);

			Bx = B_N(j,k,0);
			By = B_N(j,k,1);
			Bz = B_N(j,k,2);

			Jx = J_N(j,k,0);
			Jy = J_N(j,k,1);
			Jz = J_N(j,k,2);

			E_N(j,k,0) = -1.0*vy*Bz + vz*By + eta*Jx + (di/rho)*(Jy*Bz - Jz*By);
			E_N(j,k,1) = vx*Bz - vz*Bx + eta*Jy - (di/rho)*(Jx*Bz - Jz*Bx);
			E_N(j,k,2) = -1.0*vx*By + vy*Bx + eta*Jz + (di/rho)*(Jx*By - Jy*Bx);
			
			rho = rho_S(j,k);
			
			vx = v_S(j,k,0);
			vy = v_S(j,k,1);
			vz = v_S(j,k,2);

			Bx = B_S(j,k,0);
			By = B_S(j,k,1);
			Bz = B_S(j,k,2);

			Jx = J_S(j,k,0);
			Jy = J_S(j,k,1);
			Jz = J_S(j,k,2);

			E_S(j,k,0) = -1.0*vy*Bz + vz*By + eta*Jx + (di/rho)*(Jy*Bz - Jz*By);
			E_S(j,k,1) = vx*Bz - vz*Bx + eta*Jy - (di/rho)*(Jx*Bz - Jz*Bx);
			E_S(j,k,2) = -1.0*vx*By + vy*Bx + eta*Jz + (di/rho)*(Jx*By - Jy*Bx);

			rho = rho_E(j,k);
			
			vx = v_E(j,k,0);
			vy = v_E(j,k,1);
			vz = v_E(j,k,2);

			Bx = B_E(j,k,0);
			By = B_E(j,k,1);
			Bz = B_E(j,k,2);

			Jx = J_E(j,k,0);
			Jy = J_E(j,k,1);
			Jz = J_E(j,k,2);

			E_E(j,k,0) = -1.0*vy*Bz + vz*By + eta*Jx + (di/rho)*(Jy*Bz - Jz*By);
			E_E(j,k,1) = vx*Bz - vz*Bx + eta*Jy - (di/rho)*(Jx*Bz - Jz*Bx);
			E_E(j,k,2) = -1.0*vx*By + vy*Bx + eta*Jz + (di/rho)*(Jx*By - Jy*Bx);

			rho = rho_W(j,k);
			
			vx = v_W(j,k,0);
			vy = v_W(j,k,1);
			vz = v_W(j,k,2);

			Bx = B_W(j,k,0);
			By = B_W(j,k,1);
			Bz = B_W(j,k,2);

			Jx = J_W(j,k,0);
			Jy = J_W(j,k,1);
			Jz = J_W(j,k,2);

			E_W(j,k,0) = -1.0*vy*Bz + vz*By + eta*Jx + (di/rho)*(Jy*Bz - Jz*By);
			E_W(j,k,1) = vx*Bz - vz*Bx + eta*Jy - (di/rho)*(Jx*Bz - Jz*Bx);
			E_W(j,k,2) = -1.0*vx*By + vy*Bx + eta*Jz + (di/rho)*(Jx*By - Jy*Bx);
		}
	}
}
