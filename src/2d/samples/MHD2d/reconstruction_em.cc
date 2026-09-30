////////////////////////////////////////////////////////////////////////////////
//
// reconstruction_em.cc — CT-aware reconstruction of B and J for the EM solver
//
// B reconstruction follows the CT rule from Dr. Balbás's notes:
//   Bx: minmod in y only; at E/W (x-faces), use staggered CT value directly
//   By: minmod in x only; at N/S (y-faces), use staggered CT value directly
//   Bz: minmod in both x and y (cell-centered, no CT face)
//
// J reconstruction: full minmod in both x and y, all 3 components
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::reconstruction_em(const doublearray3d& un, const doublearray3d& Jc,
    const doublearray3d& B1_xf, const doublearray3d& B2_yf,
    doublearray3d& B_N, doublearray3d& B_S, doublearray3d& B_E, doublearray3d& B_W,
    doublearray3d& Jc_N, doublearray3d& Jc_S, doublearray3d& Jc_E, doublearray3d& Jc_W,
    const doublearray1d& dx_cell, const doublearray1d& dx_interface,
    const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha)
{
    long J = un.getIndex1Size() - 4;
    long K = un.getIndex2Size() - 4;
    long j, k;

    double ux, uy;

    for (k = 1; k < K+3; k++)
    {
        for (j = 1; j < J+3; j++)
        {
            // --- Bx: minmod in y, E/W from staggered CT variable ---
            uy = minmod3(alpha*(un(j,k+1,4) - un(j,k,4))/dy_interface(k),
                         (un(j,k+1,4) - un(j,k-1,4))/(dy_interface(k-1) + dy_interface(k)),
                         alpha*(un(j,k,4) - un(j,k-1,4))/dy_interface(k-1));

            B_N(j,k,0) = un(j,k,4) + 0.5*dy_cell(k)*uy;
            B_S(j,k,0) = un(j,k,4) - 0.5*dy_cell(k)*uy;
            B_E(j,k,0) = B1_xf(j,k,0);      // Bx at East face (j+1/2, k)
            B_W(j,k,0) = B1_xf(j-1,k,0);    // Bx at West face (j-1/2, k)

            // --- By: minmod in x, N/S from staggered CT variable ---
            ux = minmod3(alpha*(un(j+1,k,5) - un(j,k,5))/dx_interface(j),
                         (un(j+1,k,5) - un(j-1,k,5))/(dx_interface(j-1) + dx_interface(j)),
                         alpha*(un(j,k,5) - un(j-1,k,5))/dx_interface(j-1));

            B_E(j,k,1) = un(j,k,5) + 0.5*dx_cell(j)*ux;
            B_W(j,k,1) = un(j,k,5) - 0.5*dx_cell(j)*ux;
            B_N(j,k,1) = B2_yf(j,k,0);      // By at North face (j, k+1/2)
            B_S(j,k,1) = B2_yf(j,k-1,0);    // By at South face (j, k-1/2)

            // --- Bz: minmod in both x and y ---
            ux = minmod3(alpha*(un(j+1,k,6) - un(j,k,6))/dx_interface(j),
                         (un(j+1,k,6) - un(j-1,k,6))/(dx_interface(j-1) + dx_interface(j)),
                         alpha*(un(j,k,6) - un(j-1,k,6))/dx_interface(j-1));

            uy = minmod3(alpha*(un(j,k+1,6) - un(j,k,6))/dy_interface(k),
                         (un(j,k+1,6) - un(j,k-1,6))/(dy_interface(k-1) + dy_interface(k)),
                         alpha*(un(j,k,6) - un(j,k-1,6))/dy_interface(k-1));

            B_E(j,k,2) = un(j,k,6) + 0.5*dx_cell(j)*ux;
            B_W(j,k,2) = un(j,k,6) - 0.5*dx_cell(j)*ux;
            B_N(j,k,2) = un(j,k,6) + 0.5*dy_cell(k)*uy;
            B_S(j,k,2) = un(j,k,6) - 0.5*dy_cell(k)*uy;

            // --- J: full minmod in both x and y, all 3 components ---
            for (int l = 0; l < 3; l++)
            {
                ux = minmod3(alpha*(Jc(j+1,k,l) - Jc(j,k,l))/dx_interface(j),
                             (Jc(j+1,k,l) - Jc(j-1,k,l))/(dx_interface(j-1) + dx_interface(j)),
                             alpha*(Jc(j,k,l) - Jc(j-1,k,l))/dx_interface(j-1));

                uy = minmod3(alpha*(Jc(j,k+1,l) - Jc(j,k,l))/dy_interface(k),
                             (Jc(j,k+1,l) - Jc(j,k-1,l))/(dy_interface(k-1) + dy_interface(k)),
                             alpha*(Jc(j,k,l) - Jc(j,k-1,l))/dy_interface(k-1));

                Jc_N(j,k,l) = Jc(j,k,l) + 0.5*dy_cell(k)*uy;
                Jc_S(j,k,l) = Jc(j,k,l) - 0.5*dy_cell(k)*uy;
                Jc_E(j,k,l) = Jc(j,k,l) + 0.5*dx_cell(j)*ux;
                Jc_W(j,k,l) = Jc(j,k,l) - 0.5*dx_cell(j)*ux;
            }
        }
    }
}