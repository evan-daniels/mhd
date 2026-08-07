#include "centpack_2d_SD2.h"

void CENTPACK::evolution_2d_SD2(doublearray3d& un, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha, const doublearray1d& parameters, const int& id, const int& p)
{
    using namespace std;

    long j, k, l, J, K, L;

    J = un.getIndex1Size() - 4;
    K = un.getIndex2Size() - 4;
    L = un.getIndex3Size();

    doublearray3d u_N(J+4,K+4,L);
    doublearray3d u_S(J+4,K+4,L);
    doublearray3d u_E(J+4,K+4,L);
    doublearray3d u_W(J+4,K+4,L);
    doublearray3d C0(J+4,K+4,L), C1(J+4,K+4,L);
    doublearray3d jcurl(J+4,K+4,3);  // jx, jy, jz

    // --- Helper: compute jcurl from un ---
    auto compute_jcurl = [&](const doublearray3d& u, doublearray3d& jc) {
        for (k = 2; k < K+2; k++)
        {
            for (j = 2; j < J+2; j++)
            {
                double dx = dx_cell(j);
                double dy = dy_cell(k);
                jc(j,k,0) =  (u(j,k+1,6) - u(j,k-1,6))/(2.0*dy);           // jx = dBz/dy
                jc(j,k,1) = -(u(j+1,k,6) - u(j-1,k,6))/(2.0*dx);           // jy = -dBz/dx
                jc(j,k,2) =  (u(j+1,k,5) - u(j-1,k,5))/(2.0*dx)            // jz = dBy/dx - dBx/dy
                           - (u(j,k+1,4) - u(j,k-1,4))/(2.0*dy);
            }
        }
    };

    // SSP-RK stage 1
    compute_jcurl(un, jcurl);
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, jcurl, C0);

    for (l = 0; l < L; l++)
        for (j = 2; j < J+2; j++)
            for (k = 2; k < K+2; k++)
                un(j,k,l) += C0(j,k,l);

    boundary_conditions(un, parameters, id, p);

    // SSP-RK stage 2
    compute_jcurl(un, jcurl);
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, jcurl, C1);

    for (l = 0; l < L; l++)
        for (j = 2; j < J+2; j++)
            for (k = 2; k < K+2; k++)
                un(j,k,l) += 0.5*(C1(j,k,l) - C0(j,k,l));

    boundary_conditions(un, parameters, id, p);
}