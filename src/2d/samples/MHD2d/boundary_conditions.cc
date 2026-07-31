////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2005 Jorge Balbas and Eitan Tadmor
//
// boundary_conditions.cc
//
////////////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const int& id, const int& p)
{
	long J = u.getIndex1Size() - 4;
	long K = u.getIndex2Size() - 4;
	long L = u.getIndex3Size();
	long j, k, l;

	MPI::Status status;

	// y-boundaries: no communication needed

	int problem = (int) parameters(1);

	if (problem == 1)
	{
		// Reconnection Boundary Conditions

		double B0 = 1.0, nb = 0.5;
		double gamma_bc = parameters(0);
		double p_inf = 0.5 * nb;
		double U_inf = p_inf / (gamma_bc - 1.0) + 0.5 * B0 * B0;

		for (j = 2; j < J+2; j++)
		{
			// bottom ghost cells: Bx = -B0
			u(j,0,0) = nb;     u(j,1,0) = nb;
			u(j,0,1) = 0.0;    u(j,1,1) = 0.0;
			u(j,0,2) = 0.0;    u(j,1,2) = 0.0;
			u(j,0,3) = 0.0;    u(j,1,3) = 0.0;
			u(j,0,4) = -B0;    u(j,1,4) = -B0;
			u(j,0,5) = 0.0;    u(j,1,5) = 0.0;
			u(j,0,6) = 0.0;    u(j,1,6) = 0.0;
			u(j,0,7) = U_inf;  u(j,1,7) = U_inf;

			// top ghost cells: Bx = +B0
			u(j,K+2,0) = nb;    u(j,K+3,0) = nb;
			u(j,K+2,1) = 0.0;   u(j,K+3,1) = 0.0;
			u(j,K+2,2) = 0.0;   u(j,K+3,2) = 0.0;
			u(j,K+2,3) = 0.0;   u(j,K+3,3) = 0.0;
			u(j,K+2,4) = B0;    u(j,K+3,4) = B0;
			u(j,K+2,5) = 0.0;   u(j,K+3,5) = 0.0;
			u(j,K+2,6) = 0.0;   u(j,K+3,6) = 0.0;
			u(j,K+2,7) = U_inf; u(j,K+3,7) = U_inf;
		}
	}

	else if (problem == 0)

	{

		// Orszag-Tang Boundary Conditions

		for (l = 0; l < L; l++)
		{
			for (j = 2; j < J+2; j++)
			{
				u(j,0,l) = u(j,K,l);
				u(j,1,l) = u(j,K+1,l);
				u(j,K+2,l) = u(j,2,l);
				u(j,K+3,l) = u(j,3,l);
			}
		}
	}

	// x-boundaries: use Sendrecv to avoid deadlock
	if (id < p-1) {
		MPI::COMM_WORLD.Sendrecv(
			&u(J,0,0),   2*L*(K+4), MPI::DOUBLE, id+1, 2,
			&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, 1,
			status);
	}
	if (id > 0) {
		MPI::COMM_WORLD.Sendrecv(
			&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, 1,
			&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, 2,
			status);
	}
	if (id == 0) {
		MPI::COMM_WORLD.Sendrecv(
			&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, 3,
			&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, 4,
			status);
	}
	if (id == p-1) {
		MPI::COMM_WORLD.Sendrecv(
			&u(J,0,0),   2*L*(K+4), MPI::DOUBLE, 0, 4,
			&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, 0, 3,
			status);
	}
}

void CENTPACK::boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const bool& odd, const int& id, const int& p)
{
	long J = u.getIndex1Size() - 4;
	long K = u.getIndex2Size() - 4;
	long L = u.getIndex3Size();
	long j, k, l;

	MPI::Status status;

	if (odd)
	{
		// y-boundary conditions
		for (l = 0; l < L; l++)
		{
			for (j = 1; j < J+2; j++)
			{
				u(j,0,l) = u(j,K,l);
				u(j,K+2,l) = u(j,2,l);
				u(j,K+3,l) = u(j,3,l);
			}
		}

		// x-boundaries
		if (id < p-1) {
			MPI::COMM_WORLD.Sendrecv(
				&u(J,0,0),   L*(K+4),   MPI::DOUBLE, id+1, 2,
				&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, 1,
				status);
		}
		if (id > 0) {
			MPI::COMM_WORLD.Sendrecv(
				&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, 1,
				&u(0,0,0), L*(K+4),   MPI::DOUBLE, id-1, 2,
				status);
		}
		if (id == 0) {
			MPI::COMM_WORLD.Sendrecv(
				&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, 3,
				&u(0,0,0), L*(K+4),   MPI::DOUBLE, p-1, 4,
				status);
		}
		if (id == p-1) {
			MPI::COMM_WORLD.Sendrecv(
				&u(J,0,0),   L*(K+4),   MPI::DOUBLE, 0, 4,
				&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, 0, 3,
				status);
		}
	}

	else
	{
		// y-boundary conditions
		for (l = 0; l < L; l++)
		{
			for (j = 2; j < J+3; j++)
			{
				u(j,0,l) = u(j,K,l);
				u(j,1,l) = u(j,K+1,l);
				u(j,K+3,l) = u(j,3,l);
			}
		}

		// x-boundaries
		if (id < p-1) {
			MPI::COMM_WORLD.Sendrecv(
				&u(J,0,0),   2*L*(K+4), MPI::DOUBLE, id+1, 2,
				&u(J+3,0,0), L*(K+4),   MPI::DOUBLE, id+1, 1,
				status);
		}
		if (id > 0) {
			MPI::COMM_WORLD.Sendrecv(
				&u(3,0,0), L*(K+4),   MPI::DOUBLE, id-1, 1,
				&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, 2,
				status);
		}
		if (id == 0) {
			MPI::COMM_WORLD.Sendrecv(
				&u(3,0,0), L*(K+4),   MPI::DOUBLE, p-1, 3,
				&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, 4,
				status);
		}
		if (id == p-1) {
			MPI::COMM_WORLD.Sendrecv(
				&u(J,0,0),   2*L*(K+4), MPI::DOUBLE, 0, 4,
				&u(J+3,0,0), L*(K+4),   MPI::DOUBLE, 0, 3,
				status);
		}
	}
}