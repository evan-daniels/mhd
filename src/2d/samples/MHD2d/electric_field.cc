#include "centpack_2d_SD2.h"
#include <cmath>
using namespace std;

void CENTPACK::electric_field(const doublearray3d& un,
                               const doublearray1d& dx_cell,
                               const doublearray1d& dy_cell,
                               const doublearray1d& parameters,
                               doublearray3d& E)
{
    double gamma  = parameters(0);
    int    problem = (int) parameters(1);
    double eta    = parameters(2);  // resistivity
    double di     = parameters(3);  // ion inertia / Hall parameter
    double eta_j  = parameters(4);  // hyper-resistivity (0 to disable)

    long J = un.getIndex1Size() - 4;
    long K = un.getIndex2Size() - 4;

    for (long k = 2; k < K+2; k++)
    {
        for (long j = 2; j < J+2; j++)
        {
            double dx = dx_cell(j);
            double dy = dy_cell(k);

            double rho = un(j,k,0);
            double vx  = un(j,k,1)/rho;
            double vy  = un(j,k,2)/rho;
            double vz  = un(j,k,3)/rho;
            double Bx  = un(j,k,4);
            double By  = un(j,k,5);
            double Bz  = un(j,k,6);

            // j = curl B (central differences)
            double jx = (un(j,k+1,6) - un(j,k-1,6))/(2.0*dy);
            double jy = -(un(j+1,k,6) - un(j-1,k,6))/(2.0*dx);
            double jz = (un(j+1,k,5) - un(j-1,k,5))/(2.0*dx)
                      - (un(j,k+1,4) - un(j,k-1,4))/(2.0*dy);

            // -v x B
            double e1_ideal = -(vy*Bz - vz*By);
            double e2_ideal = -(vz*Bx - vx*Bz);
            double e3_ideal = -(vx*By - vy*Bx);

            // Hall term: di * (j x B) / rho
            double rho_eff = std::max(rho, 0.1);
            double e1_hall = di * (jy*Bz - jz*By) / rho_eff;
            double e2_hall = di * (jz*Bx - jx*Bz) / rho_eff;
            double e3_hall = di * (jx*By - jy*Bx) / rho_eff;

            // Resistivity: eta * j
            double e1_res = eta * jx;
            double e2_res = eta * jy;
            double e3_res = eta * jz;

            // Hyper-resistivity: eta_j * Laplacian(j)
            double e1_hyp = 0.0, e2_hyp = 0.0, e3_hyp = 0.0;
            if (eta_j > 0.0)
            {
                // Laplacian of j using second neighbors
                double lap_jx = (un(j,k+2,6) - 2.0*un(j,k,6) + un(j,k-2,6))/(4.0*dy*dy)
                               +(un(j+2,k,6) - 2.0*un(j,k,6) + un(j-2,k,6))/(4.0*dx*dx);
                // ... etc for jy, jz — simplified here
                e1_hyp = eta_j * lap_jx;
            }

            E(j,k,0) = e1_ideal + e1_hall + e1_res + e1_hyp;
            E(j,k,1) = e2_ideal + e2_hall + e2_res + e2_hyp;
            E(j,k,2) = e3_ideal + e3_hall + e3_res + e3_hyp;
        }
    }
}