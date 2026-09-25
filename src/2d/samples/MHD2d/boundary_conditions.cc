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

	int left_rank  = (id == 0)   ? p - 1 : id - 1;
	int right_rank = (id == p-1) ? 0     : id + 1;

	

	MPI_Status status;

	// y-boundaries: no communication needed
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

	/// x-boundaries require communication between processors

	// Send left two interior columns (j = 2,3)
	// Receive right ghost columns (j = J+2,J+3)

	MPI_Sendrecv(&u(2,0,0), 2*L*(K+4), MPI_DOUBLE, left_rank, 1, &u(J+2,0,0), 2*L*(K+4), MPI_DOUBLE right_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// Send right two interior columns (j = J,J+1)
	// Receive left ghost columns (j = 0,1)

	MPI_Sendrecv(&u(J,0,0), 2*L*(K+4), MPI_DOUBLE, right_rank, 2, &u(0,0,0), 2*L*(K+4), MPI_DOUBLE, left_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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