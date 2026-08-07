#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::Hy_flux_2d_SD2(const doublearray1d& u_s, const doublearray1d& u_n, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& Hy)
{
    long l;
    long L = Hy.getIndex1Size();
    double rw, re, rn, rs, b;

    doublearray1d g_n(L), g_s(L);

    spectral_radii(u_n, parameters, re, rn);
    spectral_radii(u_s, parameters, rw, rs);

    if (rn > rs)
        b = rn;
    else
        b = rs;

    flux_y(u_n, parameters, j_here, g_n);
    flux_y(u_s, parameters, j_here, g_s);

    for (l = 0; l < L; l++)
        Hy(l) = 0.5*((g_s(l) + g_n(l)) - b*(u_s(l) - u_n(l)));
}