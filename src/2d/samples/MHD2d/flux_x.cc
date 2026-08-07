#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& f)
{
    double p, p_star;
    double gamma   = parameters(0);
    int    problem  = (int) parameters(1);
    double di      = (problem == 1 && parameters.getIndex1Size() > 3) ? parameters(3) : 0.0;

    double rho = u(0);
    double vx  = u(1)/rho;
    double vy  = u(2)/rho;
    double vz  = u(3)/rho;
    double Bx  = u(4);
    double By  = u(5);
    double Bz  = u(6);

    p      = (gamma-1.0)*(u(7) - 0.5*(u(1)*vx + u(2)*vy + u(3)*vz)
           - 0.5*(Bx*Bx + By*By + Bz*Bz));
    p_star = p + 0.5*(Bx*Bx + By*By + Bz*Bz);

    // Hall electric field: E_Hall = di * (j x B) / rho
    double jx = j_here(0);
    double jy = j_here(1);
    double jz = j_here(2);

    double rho_eff = std::max(rho, 0.1);
    double Ex_hall = di * (jy*Bz - jz*By) / rho_eff;
    double Ey_hall = di * (jz*Bx - jx*Bz) / rho_eff;
    double Ez_hall = di * (jx*By - jy*Bx) / rho_eff;

    // Ideal MHD x-fluxes
    f(0) = u(1);
    f(1) = u(1)*vx + p_star - Bx*Bx;
    f(2) = u(1)*vy - Bx*By;
    f(3) = u(1)*vz - Bx*Bz;
    f(4) = 0.0;
    f(5) = By*vx - Bx*vy;
    f(6) = Bz*vx - Bx*vz;
    f(7) = (u(7) + p_star)*vx - Bx*(Bx*vx + By*vy + Bz*vz);

    // Hall corrections to induction equation in x-flux
    f(5) -= Ez_hall;   // By: -dEz/dx
    f(6) += Ey_hall;   // Bz:  dEy/dx
}