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

void CENTPACK::C_flux_em(const doublearray2d& rho_N, const doublearray2d& rho_S, const doublearray2d& rho_E, const doublearray2d& rho_W, const doublearray3d& B_N, const doublearray3d& B_S, const doublearray3d& B_E, const doublearray3d& B_W, const doublearray3d& E_N, const doublearray3d& E_S, const doublearray3d& E_E, const doublearray3d& E_W, const doublearray1d& lambda_em, const doublearray1d& mu_em, const doublearray1d& dx_interface, const doublearray1d& dy_interface, const doublearray1d& parameters, const double di, doublearray2d& Cx_em, doublearray2d& Cy_em, doublearray2d& Cz_em)
{
	
	long j, k, l, J, K, L;
	
	J = C.getIndex1Size() - 4;
	K = C.getIndex2Size() - 4;
	//L = C.getIndex3Size();

	double ax, axl, axr, ay, ayb, ayt;
	double cxm, cxp, cym, cyp;

	double rho_t, rho_b, rho_r, rho_l;

	double Bt, Bb, Br, Bl;
	
	double Ey_ejm1;
	double Ey_wj;
	double Ey_ej;
	double Ey_wjp1;

	double Ex_nkm1;
	double Ex_nk;
	double Ex_sk;
	double Ex_skp1;
	
	double Bz_ejm1;
	double Bz_ej;
	double Bz_wj;
	double Bz_wjp1;

	double Bz_nkm1;
	double Bz_nk;
	double Bz_sk;
	double Bz_skp1;
	
	double Hx_halfm, Hx_halfp;

	doublearray2d Ez_edgex(J+5,K+4);
	doublearray2d Ez_edgey(J+4,K+5);
	doublearray2d EX(J+5,K+5);

	for (j = 2; j < J+3; j++)
	{
		for (k = 1; k < K+2; k++)
		{
			Bl = sqrt(pow(B_E(j-1,k,0),2.0) + pow(B_E(j-1,k,1),2.0) + pow(B_E(j-1,k,2), 2.0));
			Br = sqrt(pow(B_W(j,k,0),2.0) + pow(B_W(j,k,1),2.0) + pow(B_W(j,k,2), 2.0));

			rho_l = rho_E(j-1,k);
			rho_r = rho_W(j,k);

			axl = di*Bl/(rho_l*dx_interface(j));
			axr = di*Br/(rho_r*dx_interface(j));
			ax = max(axl, axr);
			Ez_edgex(j,k) = 0.5*(E_W(j,k,2) + E_E(j-1,k,2)) - 0.5*ax*(B_W(j,k,1) - B_E(j-1,k,1));
		}
	}

	for (k = 2; k < K+3; k++)
	{
		for (j = 2; j < J+2; j++)
		{
			Bb = sqrt(pow(B_N(j,k-1,0),2.0) + pow(B_N(j,k-1,1),2.0) + pow(B_N(j,k-1,2), 2.0));
			Bt = sqrt(pow(B_S(j,k,0),2.0) + pow(B_S(j,k,1),2.0) + pow(B_S(j,k,2), 2.0));

			rho_b = rho_N(j,k-1);
			rho_t = rho_S(j,k);

			ayb = di*Bb/(rho_b*dy_interface(k));
			ayt = di*Bt/(rho_t*dy_interface(k));
			ay = max(ayb, ayt);
			Ez_edgey(j,k) = 0.5*(E_S(j,k,2) - E_N(j,k-1,2)) - 0.5*ay*(B_S(j,k,0) - B_N(j,k-1,0));
		}
	}

	for (j = 2; j < J+3; j++)
	{
		for (k = 2; k < K+3; k++)
			EX(j,k) = 0.5*(Ez_edgex(j,k-1) + Ez_edgex(j,k) + Ez_edgey(j-1,k) + Ez_edgey(j,k));
	}

	for (j = 2; j < J + 3; j++)
	{
		for (k = 2; k < K + 2; k++)
			Cx_em(j,k) = -1.0*mu_em(k)*(EX(j,k+1) - EX(j,k));
	}
	
	for (k = 2; k < K + 3; k++)
	{
		for (j = 2; j < J + 2; j++)
			Cy_em(j,k) = lambda(j)*(EX(j+1,k) - EX(j,k));
	}

	for (j = 2; j < J + 2; j++)
	{
		for (k = 2; k < K + 2; k++)
		{
			Ey_ejm1 = E_E(j-1,k,1);
			Ey_wj = E_W(j,k,1);

			Bz_ejm1 = B_E(j-1,k,2);
			Bz_wj = B_W(j,k,2);

			Bl = sqrt(pow(B_E(j-1,k,0),2.0) + pow(B_E(j-1,k,1),2.0) + pow(B_E(j-1,k,2), 2.0));
			Br = sqrt(pow(B_W(j,k,0),2.0) + pow(B_W(j,k,1),2.0) + pow(B_W(j,k,2), 2.0));

			rho_l = rho_E(j-1,k);
			rho_r = rho_W(j,k);

			axl = di*Bl/(rho_l*dx_interface(j));
			axr = di*Br/(rho_r*dx_interface(j));
			ax = max(axl, axr);

			Hx_halfm = 0.5*(Ey_wj + Ey_ejm1) - 0.5*ax*(Bz_wj - Bz_ejm1);
			
			Ey_ej = E_E(j,k,1);
			Ey_wjp1 = E_W(j+1,k,1);

			Bz_ej = B_E(j,k,2);
			Bz_wjp1 = B_W(j+1,k,2);

			Bl = sqrt(pow(B_E(j,k,0),2.0) + pow(B_E(j,k,1),2.0) + pow(B_E(j,k,2), 2.0));
			Br = sqrt(pow(B_W(j+1,k,0),2.0) + pow(B_W(j+1,k,1),2.0) + pow(B_W(j+1,k,2), 2.0));

			rho_l = rho_E(j,k);
			rho_r = rho_W(j+1,k);

			axl = di*Bl/(rho_l*dx_interface(j+1));
			axr = di*Br/(rho_r*dx_interface(j+1));
			ax = max(axl, axr);

			Hx_halfp = 0.5*(Ey_wjp1 + Ey_ej) - 0.5*ax*(Bz_wjp1 - Bz_ej);

			Ex_nkm1 = E_E(j,k-1,0);
			Ex_sk = E_W(j,k,0);

			Bz_nkm1 = B_N(j,k-1,2);
			Bz_sk = B_S(j,k,2);

			Bb = sqrt(pow(B_N(j,k-1,0),2.0) + pow(B_N(j,k-1,1),2.0) + pow(B_N(j,k-1,2), 2.0));
			Bt = sqrt(pow(B_S(j,k,0),2.0) + pow(B_S(j,k,1),2.0) + pow(B_S(j,k,2), 2.0));

			rho_b = rho_N(j,k-1);
			rho_t = rho_S(j,k);

			ayb = di*Bb/(rho_b*dy_interface(k));
			ayt = di*Bt/(rho_t*dy_interface(k));
			ay = max(ayb, ayt);

			Hy_halfm = -0.5*(Ex_sk + Ex_nkm1) - 0.5*ay(Bz_sk - Bz_nkm1);
			
			Ex_nk = E_N(j,k,0);
			Ex_skp1 = E_S(j,k+1,0);

			Bz_nk = B_N(j,k,2);
			Bz_skp1 = B_S(j,k+1,2);

			Bb = sqrt(pow(B_N(j,k,0),2.0) + pow(B_N(j,k,1),2.0) + pow(B_N(j,k,2), 2.0));
			Bt = sqrt(pow(B_S(j,k+1,0),2.0) + pow(B_S(j,k+1,1),2.0) + pow(B_S(j,k+1,2), 2.0));

			rho_b = rho_N(j,k);
			rho_t = rho_S(j,k+1);

			ayb = di*Bb/(rho_b*dy_interface(k+1));
			ayt = di*Bt/(rho_t*dy_interface(k+1));
			ay = max(ayb, ayt);

			Hy_halfp = -0.5*(Ex_skp1 + Ex_nk) - 0.5*ay*(Bz_skp1 - Bz_ej);
			
			Cz_em(j,k) = -1.0*lambda_em(j)*(Hx_halfp - Hx_halfm) - mu_em(j)*(Hy_halfp - Hy_halfm);
		}
	}
}
