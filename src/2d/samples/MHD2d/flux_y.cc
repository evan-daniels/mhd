#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::flux_y(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& g)
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

    // Hall electric field
    double jx = j_here(0);
    double jy = j_here(1);
    double jz = j_here(2);

    double rho_eff = std::max(rho, 0.1);
    double Ex_hall = di * (jy*Bz - jz*By) / rho_eff;
    double Ey_hall = di * (jz*Bx - jx*Bz) / rho_eff;
    double Ez_hall = di * (jx*By - jy*Bx) / rho_eff;

    // Ideal MHD y-fluxes (transport velocity is vy)
    g(0) = u(2);                           // rho*vy
    g(1) = u(1)*vy - Bx*By;               // rho*vx*vy - Bx*By
    g(2) = u(2)*vy + p_star - By*By;      // rho*vy^2 + p* - By^2
    g(3) = u(3)*vy - By*Bz;               // rho*vz*vy - By*Bz
    g(4) = Bx*vy - By*vx;                 // dBx/dt contribution
    g(5) = 0.0;                            // dBy/dt = 0 in y-flux
    g(6) = Bz*vy - By*vz;                 // dBz/dt contribution
    g(7) = (u(7) + p_star)*vy
         - By*(Bx*vx + By*vy + Bz*vz);   // energy flux

    // Hall corrections
    g(4) += Ez_hall;   // Bx: dEz/dy
    g(6) -= Ex_hall;   // Bz: -dEx/dy
}