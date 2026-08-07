#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::Hx_flux_2d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& Hx)
{
    long l;
    long L = Hx.getIndex1Size();
    double rw, re, rn, rs, a;

    doublearray1d f_e(L), f_w(L);

    spectral_radii(u_e, parameters, re, rn);
    spectral_radii(u_w, parameters, rw, rs);

    if (re > rw)
        a = re;
    else
        a = rw;

    flux_x(u_e, parameters, j_here, f_e);
    flux_x(u_w, parameters, j_here, f_w);

    for (l = 0; l < L; l++)
        Hx(l) = 0.5*((f_w(l) + f_e(l)) - a*(u_w(l) - u_e(l)));
}