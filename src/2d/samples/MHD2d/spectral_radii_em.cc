////////////////////////////////////////////////////////////////////////////////
//
// spectral_radii_em.cc — grid-scale whistler speed for the EM sub-step CFL
//
// Whistler dispersion: omega ~ di * vA * k^2  ->  phase speed ~ di * |B| * k / rho
// At the grid scale k ~ 1/dx, so:
//
//   rx = di * |B| / (rho * dx)
//   ry = di * |B| / (rho * dy)
//
// Input u uses the 8-component layout: B at u(4), u(5), u(6)
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::spectral_radii_em(const doublearray1d& u, const double& dx, const double& dy,
                                 const double& di, double& rx, double& ry)
{
    double rho     = std::max(u(0), 0.1);   // same density floor as electric_field / C_flux_em
    double Bx      = u(4);
    double By      = u(5);
    double Bz      = u(6);
    double Bmag    = sqrt(Bx*Bx + By*By + Bz*Bz);

    rx = di * Bmag / (rho * dx);
    ry = di * Bmag / (rho * dy);
}