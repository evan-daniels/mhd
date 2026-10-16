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

// ----------------------------------------------------------------------------
// Cell-centred conserved variables (SD2 solver)
//   x: periodic (MPI exchange)
//   y: problem 0 (Orszag-Tang) -> periodic
//      problem 1 (reconnection) -> perfectly conducting, impermeable walls (GEM)
// ----------------------------------------------------------------------------
void CENTPACK::boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const int& id, const int& p)
{
	long J = u.getIndex1Size() - 4;
	long K = u.getIndex2Size() - 4;
	long L = u.getIndex3Size();
	long j, l;

	MPI::Status status;

	// parameters: 0 gamma, 1 problem, 2 eta, 3 di  (eta = parameters(2), di = parameters(3) in evolution)
	int problem = (int) parameters(1);
	bool walls  = (problem == 1) && (L == 8);   // parity below assumes the 8-component un layout

	// y-boundaries: no communication needed
	for (l = 0; l < L; l++)
	{
		// odd under reflection across the wall: rho*vy (2) and By (5); everything else even
		double s = (l == 2 || l == 5) ? -1.0 : 1.0;
		for (j = 2; j < J+2; j++)
		{
			if (walls)
			{
				// mirror about the wall faces (below cell 2, above cell K+1)
				u(j,1,l)   = s*u(j,2,l);
				u(j,0,l)   = s*u(j,3,l);
				u(j,K+2,l) = s*u(j,K+1,l);
				u(j,K+3,l) = s*u(j,K,l);
			}
			else
			{
				// periodic
				u(j,0,l)   = u(j,K,l);
				u(j,1,l)   = u(j,K+1,l);
				u(j,K+2,l) = u(j,2,l);
				u(j,K+3,l) = u(j,3,l);
			}
		}
	}

	// x-boundaries: periodic, Sendrecv to avoid deadlock
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

// ----------------------------------------------------------------------------
// Staggered (odd/even) variant -- unchanged (believed to be the FD2 path; still periodic in y)
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// Staggered face fields for constrained transport
//   ywall = 0 : periodic in y (default; Orszag-Tang)
//   ywall = 1 : conducting wall, TANGENTIAL face field (B1_xf = Bx, cell-centred in y): even mirror
//   ywall = 2 : conducting wall, NORMAL face field (B2_yf = By on y-faces; face k = north face of cell k)
//               walls are faces 1 and K+1, where By = 0; odd mirror outside
//   ywall = 3 : conducting wall, cell-centred current density Jc (3 components):
//               Jx = dBz/dy and Jz = dBy/dx - dBx/dy are odd, Jy = -dBz/dx is even
// ----------------------------------------------------------------------------
void CENTPACK::boundary_conditions_ct(doublearray3d& f, const int& id, const int& p, const int& ywall)
{
	long J = f.getIndex1Size() - 4;
	long K = f.getIndex2Size() - 4;
	long L = f.getIndex3Size();   // = 1 for B1_xf / B2_yf
	long j, l;

	MPI::Status status;

	// y-boundaries: no communication needed
	for (l = 0; l < L; l++)
	{
		for (j = 2; j < J+2; j++)
		{
			if (ywall == 1)
			{
				f(j,1,l)   = f(j,2,l);
				f(j,0,l)   = f(j,3,l);
				f(j,K+2,l) = f(j,K+1,l);
				f(j,K+3,l) = f(j,K,l);
			}
			else if (ywall == 2)
			{
				f(j,1,l)   = 0.0;
				f(j,K+1,l) = 0.0;
				f(j,0,l)   = -f(j,2,l);
				f(j,K+2,l) = -f(j,K,l);
				f(j,K+3,l) = -f(j,K-1,l);
			}
			else if (ywall == 3)
			{
				double s = (l == 1) ? 1.0 : -1.0;
				f(j,1,l)   = s*f(j,2,l);
				f(j,0,l)   = s*f(j,3,l);
				f(j,K+2,l) = s*f(j,K+1,l);
				f(j,K+3,l) = s*f(j,K,l);
			}
			else
			{
				f(j,0,l)   = f(j,K,l);
				f(j,1,l)   = f(j,K+1,l);
				f(j,K+2,l) = f(j,2,l);
				f(j,K+3,l) = f(j,3,l);
			}
		}
	}

	// x-boundaries: same Sendrecv pattern as un, distinct tags (11-14) to avoid
	// any collision with un's exchange (tags 1-4)
	if (id < p-1) {
		MPI::COMM_WORLD.Sendrecv(
			&f(J,0,0),   2*L*(K+4), MPI::DOUBLE, id+1, 12,
			&f(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, 11,
			status);
	}
	if (id > 0) {
		MPI::COMM_WORLD.Sendrecv(
			&f(2,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, 11,
			&f(0,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, 12,
			status);
	}
	if (id == 0) {
		MPI::COMM_WORLD.Sendrecv(
			&f(2,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, 13,
			&f(0,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, 14,
			status);
	}
	if (id == p-1) {
		MPI::COMM_WORLD.Sendrecv(
			&f(J,0,0),   2*L*(K+4), MPI::DOUBLE, 0, 14,
			&f(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, 0, 13,
			status);
	}
}