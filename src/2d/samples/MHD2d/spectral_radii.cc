#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::spectral_radii(const doublearray1d& u, const doublearray1d& parameters, double& rx, double& ry)
{
    double rho, vx, vy, vz, p, A, B, cfx, cfy;
    double gamma = parameters(0);

    rho = u(0);
    vx  = u(1)/rho;
    vy  = u(2)/rho;
    vz  = u(3)/rho;

    p = (gamma-1.0)*(u(7) - 0.5*rho*(vx*vx + vy*vy + vz*vz)
      - 0.5*(u(4)*u(4) + u(5)*u(5) + u(6)*u(6)));
    p = std::max(p, 1e-10);

    A = gamma*p/rho;
    B = (u(4)*u(4) + u(5)*u(5) + u(6)*u(6))/rho;

    cfx = sqrt(0.5*(A + B + sqrt((A+B)*(A+B) - 4.0*A*u(4)*u(4)/rho)));
    cfy = sqrt(0.5*(A + B + sqrt((A+B)*(A+B) - 4.0*A*u(6)*u(6)/rho)));

    rx = fabs(vx) + cfx;
    ry = fabs(vy) + cfy;
}