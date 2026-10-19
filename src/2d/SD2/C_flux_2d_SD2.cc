#include "centpack_2d_SD2.h"

using namespace std;

const bool USE_FACE_AVG_J = true;   // TEMP switch: true = face-averaged J, false = cell-centred j_here on both faces

void CENTPACK::C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& jcurl, doublearray3d& C)
{
    long j, k, l, J, K, L;

    J = C.getIndex1Size() - 4;
    K = C.getIndex2Size() - 4;
    L = C.getIndex3Size();

    doublearray1d u_nkm1(L), u_n(L), u_s(L), u_skp1(L);
    doublearray1d u_ejm1(L), u_e(L), u_w(L), u_wjp1(L);
    doublearray1d j_here(3);
    doublearray1d j_w(3), j_e(3), j_s(3), j_n(3);   // face-centered J (west, east, south, north)

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

            // j_here(0) = jcurl(j,k,0);
            // j_here(1) = jcurl(j,k,1);
            // j_here(2) = jcurl(j,k,2);

            // Hx_flux_2d_SD2(u_w,    u_ejm1, parameters, j_here, Hx_halfm);
            // Hx_flux_2d_SD2(u_wjp1, u_e,    parameters, j_here, Hx_halfp);
            // Hy_flux_2d_SD2(u_s,    u_nkm1, parameters, j_here, Hy_halfm);
            // Hy_flux_2d_SD2(u_skp1, u_n,    parameters, j_here, Hy_halfp);

            // face-centered J, so neighbouring cells use the SAME flux at a shared face
            for (int c = 0; c < 3; c++)
            {
                j_w(c) = 0.5*(jcurl(j-1,k,c) + jcurl(j,  k,c));   // west face  (j-1/2)
                j_e(c) = 0.5*(jcurl(j,  k,c) + jcurl(j+1,k,c));   // east face  (j+1/2)
                j_s(c) = 0.5*(jcurl(j,k-1,c) + jcurl(j,k,  c));   // south face (k-1/2)
                j_n(c) = 0.5*(jcurl(j,k,  c) + jcurl(j,k+1,c));   // north face (k+1/2)
            }

            if (USE_FACE_AVG_J)
            {
                Hx_flux_2d_SD2(u_w,    u_ejm1, parameters, j_w, Hx_halfm);
                Hx_flux_2d_SD2(u_wjp1, u_e,    parameters, j_e, Hx_halfp);
                Hy_flux_2d_SD2(u_s,    u_nkm1, parameters, j_s, Hy_halfm);
                Hy_flux_2d_SD2(u_skp1, u_n,    parameters, j_n, Hy_halfp);
            }
            else
            {
                // previous behaviour: cell's own J on both faces
                j_here(0) = jcurl(j,k,0);
                j_here(1) = jcurl(j,k,1);
                j_here(2) = jcurl(j,k,2);

                Hx_flux_2d_SD2(u_w,    u_ejm1, parameters, j_here, Hx_halfm);
                Hx_flux_2d_SD2(u_wjp1, u_e,    parameters, j_here, Hx_halfp);
                Hy_flux_2d_SD2(u_s,    u_nkm1, parameters, j_here, Hy_halfm);
                Hy_flux_2d_SD2(u_skp1, u_n,    parameters, j_here, Hy_halfp);
            }

            for (l = 0; l < L; l++)
                C(j,k,l) = -lambda(j)*(Hx_halfp(l) - Hx_halfm(l)) - mu(k)*(Hy_halfp(l) - Hy_halfm(l));
        }
    }
}