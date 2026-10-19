#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

const bool USE_NONIDEAL_POYNTING = true;   // TEMP switch: include (E_n x B) in energy flux

void CENTPACK::flux_x(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& f)
{
    double p, p_star;
    double gamma = parameters(0);

    double rho = u(0);
    double vx  = u(1)/rho;
    double vy  = u(2)/rho;
    double vz  = u(3)/rho;
    double Bx  = u(4);
    double By  = u(5);
    double Bz  = u(6);
    
    // Non-ideal Poynting flux (Eq. 1.1d): E_n = (di/rho) J x B + eta J
    const double eta = parameters(2), di = parameters(3);
    const double Jx = j_here(0), Jy = j_here(1), Jz = j_here(2);
    const double rho_eff = std::max(rho, 0.1);          // same floor as electric_field / C_flux_em
    const double Enx = di/rho_eff*(Jy*Bz - Jz*By) + eta*Jx;
    const double Eny = di/rho_eff*(Jz*Bx - Jx*Bz) + eta*Jy;
    const double Enz = di/rho_eff*(Jx*By - Jy*Bx) + eta*Jz;

    p      = (gamma-1.0)*(u(7) - 0.5*(u(1)*vx + u(2)*vy + u(3)*vz)
           - 0.5*(Bx*Bx + By*By + Bz*Bz));
    p_star = p + 0.5*(Bx*Bx + By*By + Bz*Bz);

    // Ideal MHD x-fluxes (B components 4-6 are computed but not applied;
    // B is evolved by the EM solver)
    f(0) = u(1);
    f(1) = u(1)*vx + p_star - Bx*Bx;
    f(2) = u(1)*vy - Bx*By;
    f(3) = u(1)*vz - Bx*Bz;
    f(4) = 0.0;
    f(5) = By*vx - Bx*vy;
    f(6) = Bz*vx - Bx*vz;
    // f(7) = (u(7) + p_star)*vx - Bx*(Bx*vx + By*vy + Bz*vz);
    f(7) = (u(7) + p_star)*vx - Bx*(Bx*vx + By*vy + Bz*vz);
    if (USE_NONIDEAL_POYNTING) {
        f(7) += (Eny*Bz - Enz*By);               // + (E_n x B)_x
    }
}