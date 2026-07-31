#include "centpack_2d_SD2.h"
#include <cmath>
using namespace std;

void CENTPACK::resistivity_step(doublearray3d& un, const doublearray1d& dx_cell,
                                const doublearray1d& dy_cell, double dt,
                                const doublearray1d& parameters)
{
    int problem = (int) parameters(1);
    if (problem != 1) return;  // only apply for reconnection

    double eta = parameters(2);  // resistivity from input
    long J = un.getIndex1Size() - 4;
    long K = un.getIndex2Size() - 4;

    doublearray3d un_new = un;

    for (long k = 2; k < K+2; k++)
    {
        for (long j = 2; j < J+2; j++)
        {
            double dx = dx_cell(j);
            double dy = dy_cell(k);

            double lapBx = (un(j+1,k,4) - 2.0*un(j,k,4) + un(j-1,k,4))/(dx*dx)
                         + (un(j,k+1,4) - 2.0*un(j,k,4) + un(j,k-1,4))/(dy*dy);

            double lapBy = (un(j+1,k,5) - 2.0*un(j,k,5) + un(j-1,k,5))/(dx*dx)
                         + (un(j,k+1,5) - 2.0*un(j,k,5) + un(j,k-1,5))/(dy*dy);

            double lapBz = (un(j+1,k,6) - 2.0*un(j,k,6) + un(j-1,k,6))/(dx*dx)
                         + (un(j,k+1,6) - 2.0*un(j,k,6) + un(j,k-1,6))/(dy*dy);

            // curl B (current density jz = dBy/dx - dBx/dy)
            double jz = (un(j+1,k,5) - un(j-1,k,5))/(2.0*dx)
                      - (un(j,k+1,4) - un(j,k-1,4))/(2.0*dy);

            un_new(j,k,4) = un(j,k,4) + dt * eta * lapBx;
            un_new(j,k,5) = un(j,k,5) + dt * eta * lapBy;
            un_new(j,k,6) = un(j,k,6) + dt * eta * lapBz;
            un_new(j,k,7) = un(j,k,7) + dt * eta * jz * jz;  // Ohmic heating
        }
    }

    un = un_new;
}

// void CENTPACK::hall_step(doublearray3d& un, const doublearray1d& dx_cell,
//                          const doublearray1d& dy_cell, double dt,
//                          const doublearray1d& parameters)
// {
//     int problem = (int) parameters(1);
//     if (problem != 1) return;

//     double di = parameters(3);  // ion inertia length delta_i/L0

//     long J = un.getIndex1Size() - 4;
//     long K = un.getIndex2Size() - 4;

//     doublearray3d un_new = un;

//     for (long k = 2; k < K+2; k++)
//     {
//         for (long j = 2; j < J+2; j++)
//         {
//             double dx = dx_cell(j);
//             double dy = dy_cell(k);

//             // Ez at j+1,k
//             double Ez_jp1, Ez_jm1, Ez_kp1, Ez_km1;
//             {
//                 long jj=j+1, kk=k;
//                 double rr=un(jj,kk,0);
//                 double jxx=(un(jj,kk+1,6)-un(jj,kk-1,6))/(2.0*dy);
//                 double jyy=-(un(jj+1,kk,6)-un(jj-1,kk,6))/(2.0*dx);
//                 Ez_jp1=di*(jxx*un(jj,kk,5)-jyy*un(jj,kk,4))/rr;
//             }
//             // Ez at j-1,k
//             {
//                 long jj=j-1, kk=k;
//                 double rr=un(jj,kk,0);
//                 double jxx=(un(jj,kk+1,6)-un(jj,kk-1,6))/(2.0*dy);
//                 double jyy=-(un(jj+1,kk,6)-un(jj-1,kk,6))/(2.0*dx);
//                 Ez_jm1=di*(jxx*un(jj,kk,5)-jyy*un(jj,kk,4))/rr;
//             }
//             // Ez at j,k+1
//             {
//                 long jj=j, kk=k+1;
//                 double rr=un(jj,kk,0);
//                 double jxx=(un(jj,kk+1,6)-un(jj,kk-1,6))/(2.0*dy);
//                 double jyy=-(un(jj+1,kk,6)-un(jj-1,kk,6))/(2.0*dx);
//                 Ez_kp1=di*(jxx*un(jj,kk,5)-jyy*un(jj,kk,4))/rr;
//             }
//             // Ez at j,k-1
//             {
//                 long jj=j, kk=k-1;
//                 double rr=un(jj,kk,0);
//                 double jxx=(un(jj,kk+1,6)-un(jj,kk-1,6))/(2.0*dy);
//                 double jyy=-(un(jj+1,kk,6)-un(jj-1,kk,6))/(2.0*dx);
//                 Ez_km1=di*(jxx*un(jj,kk,5)-jyy*un(jj,kk,4))/rr;
//             }

//             double dEz_dx = (Ez_jp1 - Ez_jm1)/(2.0*dx);
//             double dEz_dy = (Ez_kp1 - Ez_km1)/(2.0*dy);

//             un_new(j,k,4) += dt * dEz_dy;   // dBx/dt += dEz/dy
//             un_new(j,k,5) -= dt * dEz_dx;   // dBy/dt -= dEz/dx
//         }
//     }

//     un = un_new;
// }