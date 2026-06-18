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
// C_flux.cc -- a core function of CentPack
//
// Function called by evolution
//
// Requires: Hx_flux_2d_SD2 and Hy_flux_2d_SD2 (see prototypes below)
//
// This function calculates the SSP-RK fluxes employed for the evolution of the
// cell averages of u
//
// Input (passed by reference):
//
// (1) Eight doublearray3d type variables:
//
//     u_N  --  holds the point values of u at the cell interfaces
//                (x_j, y_k+1/2) over the entire solution domain
//
//     u_S  --  holds the point values of u at the cell interfaces
//                (x_j, y_k-1/2) over the entire solution domain
//
//     u_E  --  holds the point values of u at the cell interfaces
//                (x_j+1/2, y_k) over the entire solution domain
//
//     u_W  --  holds the point values of u at the cell interfaces
//                (x_j-1/2, y_k) over the entire solution domain
//
// (2) lambda  --  a double type variable holding the mesh ratio dt/dx
//
// (3) mu  --  a double type variable holding the mesh ratio dt/dy
//
// (4) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) C -- a doublearray3d type variable holding the value of the SSP-RK fluxes
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& Efield, doublearray3d& C)
{
	
	long j, k, l, J, K, L;
	
	J = C.getIndex1Size() - 4;
	K = C.getIndex2Size() - 4;
	L = C.getIndex3Size();
	
	doublearray1d u_nkm1(L);
	doublearray1d u_n(L);
	doublearray1d u_s(L);
	doublearray1d u_skp1(L);
	doublearray1d u_ejm1(L);
	doublearray1d u_e(L);
	doublearray1d u_w(L);
	doublearray1d u_wjp1(L);
	doublearray1d e_here_m(3);
	doublearray1d e_here_p(3);
	
	doublearray1d Hx_halfp(L), Hx_halfm(L), Hy_halfp(L), Hy_halfm(L);

	for (k = 2; k < K + 2; k++)
	{
		for (j = 2; j < J + 2; j++)
		{
			for (l = 0; l < L; l++)
			{
				u_nkm1(l) = u_N(j,k-1,l);
				u_n(l) = u_N(j,k,l);
				u_s(l) = u_S(j,k,l);
				u_skp1(l) = u_S(j,k+1,l);
				u_ejm1(l) = u_E(j-1,k,l);
				u_e(l) = u_E(j,k,l);
				u_w(l) = u_W(j,k,l);
				u_wjp1(l) = u_W(j+1,k,l);
			}
		
			// x direction
			e_here_m(0) = Efield(j,k,0);
			e_here_m(1) = Efield(j,k,1);
			e_here_m(2) = Efield(j,k,2);

			e_here_p(0) = Efield(j+1,k,0);
			e_here_p(1) = Efield(j+1,k,1);
			e_here_p(2) = Efield(j+1,k,2);

			Hx_flux_2d_SD2(u_w, u_ejm1, parameters, e_here_m, Hx_halfm);
			Hx_flux_2d_SD2(u_wjp1, u_e, parameters, e_here_p, Hx_halfp);

			// y direction 
			e_here_m(0) = Efield(j,k,0);
			e_here_m(1) = Efield(j,k,1);
			e_here_m(2) = Efield(j,k,2);

			e_here_p(0) = Efield(j,k+1,0);
			e_here_p(1) = Efield(j,k+1,1);
			e_here_p(2) = Efield(j,k+1,2);

			Hy_flux_2d_SD2(u_s, u_nkm1, parameters, e_here_m, Hy_halfm);
			Hy_flux_2d_SD2(u_skp1, u_n, parameters, e_here_p, Hy_halfp);

			// src below
			// Hx_flux_2d_SD2(u_w, u_ejm1, parameters, e_here, Hx_halfm);
			// Hx_flux_2d_SD2(u_wjp1, u_e, parameters, e_here, Hx_halfp);
			// Hy_flux_2d_SD2(u_s, u_nkm1, parameters, e_here, Hy_halfm);
			// Hy_flux_2d_SD2(u_skp1, u_n, parameters, e_here, Hy_halfp);
			// src above
			
			for (l = 0; l < L; l++)
				C(j,k,l) = -lambda(j)*(Hx_halfp(l) - Hx_halfm(l)) - mu(k)*(Hy_halfp(l) - Hy_halfm(l));
		}
	}
}
