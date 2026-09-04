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

void CENTPACK::evolution_2d_SD2(doublearray3d& un, doublearray2d& Bx, doublearray2d& By, doublearray2d& Ez, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha, const doublearray1d& parameters, const double& dt_fluid, const int& id, const int& p)
{
	using namespace std;

	int r;
	long j, k, l, m, J, K, L;
	double dt_em, t_em;
	
	J = un.getIndex1Size() - 4;
	K = un.getIndex2Size() - 4;
	L = un.getIndex3Size();
	
	doublearray1d lambda_em(J+4);
	doublearray1d mu_em(K+4);

	doublearray2d Bz(J+4,K+4);
	doublearray2d E_X(J+5,K+5);
	doublearray3d v(J+4,K+4,3);
	doublearray3d B(J+4,K+4,3);
	doublearray3d Jc(J+4,K+4,3);

	doublearray2d rho_N(J+4,K+4);
	doublearray2d rho_S(J+4,K+4);
	doublearray2d rho_E(J+4,K+4);
	doublearray2d rho_W(J+4,K+4);
	
	doublearray3d u_N(J+4,K+4,5);
	doublearray3d u_S(J+4,K+4,5);
	doublearray3d u_E(J+4,K+4,5);
	doublearray3d u_W(J+4,K+4,5);
	
	doublearray3d v_N(J+4,K+4,3);
	doublearray3d v_S(J+4,K+4,3);
	doublearray3d v_E(J+4,K+4,3);
	doublearray3d v_W(J+4,K+4,3);
	doublearray3d B_N(J+4,K+4,3);
	doublearray3d B_S(J+4,K+4,3);
	doublearray3d B_E(J+4,K+4,3);
	doublearray3d B_W(J+4,K+4,3);
	doublearray3d J_N(J+4,K+4,3);
	doublearray3d J_S(J+4,K+4,3);
	doublearray3d J_E(J+4,K+4,3);
	doublearray3d J_W(J+4,K+4,3);
	doublearray3d E_N(J+4,K+4,3);
	doublearray3d E_S(J+4,K+4,3);
	doublearray3d E_E(J+4,K+4,3);
	doublearray3d E_W(J+4,K+4,3);
	doublearray3d C0(J+4,K+4,5), C1(J+4,K+4,5);
	
	doublearray2d Cx_em(J+5,K+4);
	doublearray2d Cy_em(J+4,K+5);
	doublearray2d Cz_em(J+4,K+4);

	for (j = 0; j < J+4; j++)
	{
		for (k = 0; k < K+4; k++)
		{
			B(j,k,0) = un(j,k,5);
			B(j,k,1) = un(j,k,6);
			B(j,k,2) = un(j,k,7);

			Jc(j,k,0) = un(j,k,8);
			Jc(j,k,1) = un(j,k,9);
			Jc(j,k,2) = un(j,k,10);
		}
	}
	
	// reconstruction of fluid variables
	reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);

	for (j = 1; j < J+3; j++)
	{
		for (k = 1; k < K+3; k++)
		{
			rho_N(j,k) = u_N(j,k,0);
			rho_S(j,k) = u_S(j,k,0);
			rho_E(j,k) = u_E(j,k,0);
			rho_W(j,k) = u_W(j,k,0);
			
			v_N(j,k,0) = u_N(j,k,1)/u_N(j,k,0);
			v_N(j,k,1) = u_N(j,k,2)/u_N(j,k,0);
			v_N(j,k,2) = u_N(j,k,3)/u_N(j,k,0);

			v_S(j,k,0) = u_S(j,k,1)/u_S(j,k,0);
			v_S(j,k,1) = u_S(j,k,2)/u_S(j,k,0);
			v_S(j,k,2) = u_S(j,k,3)/u_S(j,k,0);

			v_E(j,k,0) = u_E(j,k,1)/u_E(j,k,0);
			v_E(j,k,1) = u_E(j,k,2)/u_E(j,k,0);
			v_E(j,k,2) = u_E(j,k,3)/u_E(j,k,0);

			v_W(j,k,0) = u_W(j,k,1)/u_W(j,k,0);
			v_W(j,k,1) = u_W(j,k,2)/u_W(j,k,0);
			v_W(j,k,2) = u_W(j,k,3)/u_W(j,k,0);
		}
	}

	// reconstruction of EM variables
	reconstruction_em(B, Bx, By, B_N, B_S, B_E, B_W, Jc, J_N, J_S, J_E, J_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
	r = 0;
	
	// Calculation of Electric Field
	electric_field(rho_N, rho_S, rho_E, rho_W, v_N, v_S, v_E, v_W, B_N, B_S, B_E, B_W, J_N, J_S, J_E, J_W, E_N, E_S, E_E, E_W);

	C_flux_2d_SD2(u_N, u_S, u_E, u_W, B_N, B_S, B_E, B_W, E_N, E_S, E_E, E_W, lambda, mu, parameters, C0);
	
	for (l = 0; l < 5; l++)
	{
		for (j = 2; j < J + 2; j++)
		{
			for (k = 2; k < K + 2; k++)
				un(j,k,l) += C0(j,k,l);
		}
	}
	
	boundary_conditions(un, parameters, id, p);

	// fluid variables have been updated over 0.5*dt

	// apply several time substeps (size dt_em) to evolve EM variables over 
	// dt, we evolve B and recover J and E from B at every substep

	m = 0;
	t_em = 0.0;
	r = 0;

	do
	{
		time_step_em(un, B, Jc, dx_cell, dy_cell, cfl_em, dtp_em, dt_fluid, parameters);

		if (std::isnan(dtp_em) || std::isinf(dtp_em))
			MPI_Abort(MPI_COMM_WORLD, 1);

		MPI_Allreduce(&dtp_em, &dt_em, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		
		dt_em = min(dt_em, 0.5*dt_fluid - t_em);
		
		for (j = 2; j < J + 2; j++)
			lambda_em(j) = dt_em/dx_cell(j);
	
		for (k = 2; k < K + 2; k++)
			mu_em(k) = dt_em/dy_cell(k);

		t_em += dt_em;
		m++;
		
		if (r == 1)
		{
			reconstruction_em(B, Bx, By, B_N, B_S, B_E, B_W, J, J_N, J_S, J_E, J_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);

			electric_field(rho_N, rho_S, rho_E, rho_W, v_N, v_S, v_E, v_W, B_N, B_S, B_E, B_W, J_N, J_S, J_E, J_W, E_N, E_S, E_E, E_W);
		}

		else
			r = 1;
		
		C_flux_em(B_N, B_S, B_E, B_W, E_N, E_S, E_E, E_W, lambda_em, mu_em, dx_interface, dy_interface, dt_em, parameters, Cx_em, Cy_em, Cz_em);

		for (j = 2; j < J+3; j++)
		{
			for (k = 2; k < K+2; k++)
				Bx(j,k) += Cx_em(j,k);
		}

		for (k = 2; k < K+3; k++)
		{
			for (j = 2; j < J+2; j++)
				By(j,k) += Cy_em(j,k);
		}

		for (k = 2; k < K+2; k++)
		{
			for (j = 2; j < J+2; j++)
				Bz(j,k) += Cz_em(j,k);
		}

		for (j = 2; j < J+2; j++)
		{
			for (k = 2; k < K+2; k++)
			{
				B(j,k,0) = 0.5*(Bx(j,k) + Bx(j+1,k));
				B(j,k,1) = 0.5*(By(j,k) + By(j,k+1));
				B(j,k,2) = Bz(j,k);
			}
		}

		boundary_conditions_em(B, parameters, id, p);
		current(B, Jc, id, p);

	}while (t_em < 0.5*dt_fluid);

	// Push B and Jc back into un

	for (j = 0; j <  J+4; j++)
	{
		for (j = 0; j <  J+4; j++)
		{
			un(j,k,5) = B(j,k,0);
			un(j,k,6) = B(j,k,1);
			un(j,k,7) = B(j,k,2);

			un(j,k,8) = Jc(j,k,0);
			un(j,k,9) = Jc(j,k,1);
			un(j,k,10) = Jc(j,k,2);
		}
	}

	reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
	
	for (j = 1; j < J+3; j++)
	{
		for (k = 1; k < K+3; k++)
		{
			rho_N(j,k) = u_N(j,k,0);
			rho_S(j,k) = u_S(j,k,0);
			rho_E(j,k) = u_E(j,k,0);
			rho_W(j,k) = u_W(j,k,0);

			v_N(j,k,0) = u_N(j,k,1)/u_N(j,k,0);
			v_N(j,k,1) = u_N(j,k,2)/u_N(j,k,0);
			v_N(j,k,2) = u_N(j,k,3)/u_N(j,k,0);

			v_S(j,k,0) = u_S(j,k,1)/u_S(j,k,0);
			v_S(j,k,1) = u_S(j,k,2)/u_S(j,k,0);
			v_S(j,k,2) = u_S(j,k,3)/u_S(j,k,0);

			v_E(j,k,0) = u_E(j,k,1)/u_E(j,k,0);
			v_E(j,k,1) = u_E(j,k,2)/u_E(j,k,0);
			v_E(j,k,2) = u_E(j,k,3)/u_E(j,k,0);

			v_W(j,k,0) = u_W(j,k,1)/u_W(j,k,0);
			v_W(j,k,1) = u_W(j,k,2)/u_W(j,k,0);
			v_W(j,k,2) = u_W(j,k,3)/u_W(j,k,0);
		}
	}

	// reconstruction of EM variables
	reconstruction_em(B, Bx, By, B_N, B_S, B_E, B_W, Jc, J_N, J_S, J_E, J_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
	r = 0;
	
	// Calculation of Electric Field
	electric_field(rho_N, rho_S, rho_E, rho_W, v_N, v_S, v_E, v_W, B_N, B_S, B_E, B_W, J_N, J_S, J_E, J_W, E_N, E_S, E_E, E_W);

	C_flux_2d_SD2(u_N, u_S, u_E, u_W, B_N, B_S, B_E, B_W, E_N, E_S, E_E, E_W, lambda, mu, parameters, C1);
	
	for (l = 0; l < 5; l++)
	{
		for (j = 2; j < J + 2; j++)
		{
			for (k = 2; k < K + 2; k++)
				un(j,k,l) += 0.5*(C1(j,k,l) - C0(j,k,l));
		}
	}

	boundary_conditions(un, parameters, id, p);

	// fluid variables have been updated over 0.5*dt

	// apply several time substeps (size dt_em) to evolve EM variables over 
	// dt, we evolve B and recover J and E from B at every substep

	m = 0;
	r = 0;

	do
	{
		time_step_em(un, B, Jc, dx_cell, dy_cell, cfl_em, dtp_em, dt_fluid, parameters);

		if (std::isnan(dtp_em) || std::isinf(dtp_em))
			MPI_Abort(MPI_COMM_WORLD, 1);

		MPI_Allreduce(&dtp_em, &dt_em, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		
		dt_em = min(dt_em, dt_fluid - t_em);
		
		for (j = 2; j < J + 2; j++)
			lambda_em(j) = dt_em/dx_cell(j);
	
		for (k = 2; k < K + 2; k++)
			mu_em(k) = dt_em/dy_cell(k);

		t_em += dt_em;
		m++;
		
		if (r == 1)
		{
			reconstruction_em(B, Bx, By, B_N, B_S, B_E, B_W, Jc, J_N, J_S, J_E, J_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);

			electric_field(rho_N, rho_S, rho_E, rho_W, v_N, v_S, v_E, v_W, B_N, B_S, B_E, B_W, J_N, J_S, J_E, J_W, E_N, E_S, E_E, E_W);
		}

		else
			r = 1;
		
		C_flux_em(B_N, B_S, B_E, B_W, E_N, E_S, E_E, E_W, lambda_em, mu_em, dx_interface, dy_interface, parameters, Cx_em, Cy_em, Cz_em);

		for (j = 2; j < J+3; j++)
		{
			for (k = 2; k < K+2; k++)
				Bx(j,k) += Cx_em(j,k);
		}

		for (k = 2; k < K+3; k++)
		{
			for (j = 2; j < J+2; j++)
				By(j,k) += Cy_em(j,k);
		}

		for (k = 2; k < K+2; k++)
		{
			for (j = 2; j < J+2; j++)
				Bz(j,k) += Cz_em(j,k);
		}

		for (j = 2; j < J + 2; j++)
		{
			for (k = 2; k < K +2; k++)
			{
				B(j,k,0) = 0.5*(Bx(j,k) + Bx(j+1,k));
				B(j,k,1) = 0.5*(By(j,k) + By(j,k+1));
				B(j,k,2) = Bz(j,k);
			}
		}

		boundary_conditions_em(B, parameters, id, p);
		current(B, J, id, p);

	}while (t_em < dt_fluid);

	// Push B and Jc back into un

	for (j = 0; j <  J+4; j++)
	{
		for (j = 0; j <  J+4; j++)
		{
			un(j,k,5) = B(j,k,0);
			un(j,k,6) = B(j,k,1);
			un(j,k,7) = B(j,k,2);

			un(j,k,8) = Jc(j,k,0);
			un(j,k,9) = Jc(j,k,1);
			un(j,k,10) = Jc(j,k,2);
		}
	}
}
