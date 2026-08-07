#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& jcurl, doublearray3d& C)
{
    long j, k, l, J, K, L;

    J = C.getIndex1Size() - 4;
    K = C.getIndex2Size() - 4;
    L = C.getIndex3Size();

    doublearray1d u_nkm1(L), u_n(L), u_s(L), u_skp1(L);
    doublearray1d u_ejm1(L), u_e(L), u_w(L), u_wjp1(L);
    doublearray1d j_here(3);

    doublearray1d Hx_halfp(L), Hx_halfm(L), Hy_halfp(L), Hy_halfm(L);

    for (k = 2; k < K+2; k++)
    {
        for (j = 2; j < J+2; j++)
        {
            for (l = 0; l < L; l++)
            {
                u_nkm1(l) = u_N(j,k-1,l);
                u_n(l)    = u_N(j,k,l);
                u_s(l)    = u_S(j,k,l);
                u_skp1(l) = u_S(j,k+1,l);
                u_ejm1(l) = u_E(j-1,k,l);
                u_e(l)    = u_E(j,k,l);
                u_w(l)    = u_W(j,k,l);
                u_wjp1(l) = u_W(j+1,k,l);
            }

            // Pass local j to flux functions
            j_here(0) = jcurl(j,k,0);
            j_here(1) = jcurl(j,k,1);
            j_here(2) = jcurl(j,k,2);

            Hx_flux_2d_SD2(u_w,    u_ejm1, parameters, j_here, Hx_halfm);
            Hx_flux_2d_SD2(u_wjp1, u_e,    parameters, j_here, Hx_halfp);
            Hy_flux_2d_SD2(u_s,    u_nkm1, parameters, j_here, Hy_halfm);
            Hy_flux_2d_SD2(u_skp1, u_n,    parameters, j_here, Hy_halfp);

            for (l = 0; l < L; l++)
                C(j,k,l) = -lambda(j)*(Hx_halfp(l) - Hx_halfm(l)) - mu(k)*(Hy_halfp(l) - Hy_halfm(l));
        }
    }
}