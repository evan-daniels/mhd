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

#include"centpack_2d_SD2.h"

using namespace std;

void CENTPACK::boundary_conditions_em(doublearray3d& B, const doublearray1d& parameters, const int& id, const int& p)
{
	void CENTPACK::boundary_conditions_em(doublearray3d& B,
                                      const doublearray1d& parameters,
                                      const int& id,
                                      const int& p)
{
    long J = B.getIndex1Size() - 4;
    long K = B.getIndex2Size() - 4;
    long L = B.getIndex3Size();   // should be 3

    long j, l;

	int left_rank  = (id == 0)   ? p - 1 : id - 1;
    int right_rank = (id == p-1) ? 0     : id + 1;

    // Periodic boundaries in y

    for (l = 0; l < L; l++)
    {
        for (j = 2; j < J + 2; j++)
        {
            B(j,0,l)   = B(j,K,l);
            B(j,1,l)   = B(j,K+1,l);

            B(j,K+2,l) = B(j,2,l);
            B(j,K+3,l) = B(j,3,l);
        }
    }

    // Periodic processor boundaries in x

    // Send columns j=2,3 left
    // Receive j=J+2,J+3 from right

    MPI_Sendrecv(&B(2,0,0), 2*L*(K+4), MPI_DOUBLE, left_rank, 1, &B(J+2,0,0), 2*L*(K+4), MPI_DOUBLE, right_rank, 1, MPI_COMM_WORLD MPI_STATUS_IGNORE);

    // Send columns j=J,J+1 right
    // Receive j=0,1 from left

    MPI_Sendrecv(&B(J,0,0), 2*L*(K+4), MPI_DOUBLE, right_rank, 2, &B(0,0,0), 2*L*(K+4), MPI_DOUBLE, left_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}