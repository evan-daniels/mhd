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

void CENTPACK::compute_current(const doublearray3d& B, doublearray3d& Jc, const doublearray1d& dx_interface, const doublearray1d& dy_interface)
{
    long J = B.getIndex1Size() - 4;
    long K = B.getIndex2Size() - 4;

    long j, k;

    for (j = 1; j < J + 3; j++)
    {
        for (k = 1; k < K + 3; k++)
        {
            // Jx = dBz/dy

            Jc(j,k,0) = (B(j,k+1,2) - B(j,k-1,2))/(dy_interface(k) + dy_interface(k-1));

            // Jy = -dBz/dx

            Jc(j,k,1) = -(B(j+1,k,2) - B(j-1,k,2))/(dx_interface(j) + dx_interface(j-1));

            // Jz = dBy/dx - dBx/dy

            Jc(j,k,2) = (B(j+1,k,1) - B(j-1,k,1))/(dx_interface(j) + dx_interface(j-1)) - (B(j,k+1,0) - B(j,k-1,0))/(dy_interface(k) + dy_interface(k-1));
        }
    }
}