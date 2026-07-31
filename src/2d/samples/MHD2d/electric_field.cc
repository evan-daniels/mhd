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
    // double di     = parameters(3);  // ion inertia / Hall parameter

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

            // Resistivity: eta * j
            double e1_res = eta * jx;
            double e2_res = eta * jy;
            double e3_res = eta * jz;

            E(j,k,0) = e1_ideal + e1_res;
            E(j,k,1) = e2_ideal + e2_res;
            E(j,k,2) = e3_ideal + e3_res;
        }
    }
}