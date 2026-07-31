#include "centpack_2d_SD2.h"
#include <cmath>
using namespace std;

// Apply Hall operator L to B fields, storing result in out
// L(B) = [-dEz/dy for Bx, dEz/dx for By]
static void apply_hall_operator(const doublearray3d& B,
                                 const doublearray3d& un,
                                 const doublearray1d& dx_cell,
                                 const doublearray1d& dy_cell,
                                 double di,
                                 doublearray3d& out)
{
    long J = B.getIndex1Size() - 4;
    long K = B.getIndex2Size() - 4;

    for (long k = 2; k < K+2; k++)
    {
        for (long j = 2; j < J+2; j++)
        {
            double dx = dx_cell(j);
            double dy = dy_cell(k);

            // Compute Ez at neighbors using current B iterate
            auto Ez = [&](long jj, long kk) {
                double rho = std::max(un(jj,kk,0), 0.1);
                double jx  =  (B(jj,kk+1,6) - B(jj,kk-1,6))/(2.0*dy);
                double jy  = -(B(jj+1,kk,6) - B(jj-1,kk,6))/(2.0*dx);
                return di * (jx*B(jj,kk,5) - jy*B(jj,kk,4)) / rho;
            };

            double dEz_dx = (Ez(j+1,k) - Ez(j-1,k))/(2.0*dx);
            double dEz_dy = (Ez(j,k+1) - Ez(j,k-1))/(2.0*dy);

            out(j,k,4) =  dEz_dy;   // L applied to Bx
            out(j,k,5) = -dEz_dx;   // L applied to By
        }
    }
}

// Compute dot product of two 3D arrays over Bx, By components only
static double dot_product(const doublearray3d& a, const doublearray3d& b,
                           long J, long K)
{
    double sum = 0.0;
    for (long k = 2; k < K+2; k++)
        for (long j = 2; j < J+2; j++)
        {
            sum += a(j,k,4)*b(j,k,4);  // Bx
            sum += a(j,k,5)*b(j,k,5);  // By
            sum += a(j,k,6)*b(j,k,6);  // Bz
        }
    return sum;
}

void CENTPACK::hall_step(doublearray3d& un,
                                   const doublearray1d& dx_cell,
                                   const doublearray1d& dy_cell,
                                   double dt,
                                   const doublearray1d& parameters)
{
    int problem = (int) parameters(1);
    if (problem != 1) return;

    double di = parameters(3);

    long J = un.getIndex1Size() - 4;
    long K = un.getIndex2Size() - 4;

    // b = right hand side = current B (after ideal MHD step)
    doublearray3d b = un;

    // x = initial guess = current B
    doublearray3d x = un;

    // Temporaries for CG
    doublearray3d Ax(J+4, K+4, 8);   // A*x
    doublearray3d Ap(J+4, K+4, 8);   // A*p
    doublearray3d r(J+4, K+4, 8);    // residual
    doublearray3d p(J+4, K+4, 8);    // search direction
    doublearray3d Lx(J+4, K+4, 8);   // L*x

    // Compute A*x = x - dt*L*x
    apply_hall_operator(x, un, dx_cell, dy_cell, di, Lx);
    for (long k = 2; k < K+2; k++)
        for (long j = 2; j < J+2; j++)
        {
            Ax(j,k,4) = x(j,k,4) - dt*Lx(j,k,4);
            Ax(j,k,5) = x(j,k,5) - dt*Lx(j,k,5);
            Ax(j,k,6) = x(j,k,6) - dt*Lx(j,k,6);
        }

    // r = b - A*x
    for (long k = 2; k < K+2; k++)
        for (long j = 2; j < J+2; j++)
        {
            r(j,k,4) = b(j,k,4) - Ax(j,k,4);
            r(j,k,5) = b(j,k,5) - Ax(j,k,5);
            r(j,k,6) = b(j,k,6) - Ax(j,k,6);
            p(j,k,4) = r(j,k,4);
            p(j,k,5) = r(j,k,5);
            p(j,k,6) = r(j,k,6);
        }

    double r_dot_r = dot_product(r, r, J, K);

    int max_iter = 200;
    double tol   = 1e-10;

    for (int iter = 0; iter < max_iter; iter++)
    {
        if (sqrt(r_dot_r) < tol) break;

        // Ap = A*p = p - dt*L*p
        doublearray3d Lp(J+4, K+4, 8);
        apply_hall_operator(p, un, dx_cell, dy_cell, di, Lp);
        for (long k = 2; k < K+2; k++)
            for (long j = 2; j < J+2; j++)
            {
                Ap(j,k,4) = p(j,k,4) - dt*Lp(j,k,4);
                Ap(j,k,5) = p(j,k,5) - dt*Lp(j,k,5);
                Ap(j,k,6) = p(j,k,6) - dt*Lp(j,k,6);
            }

        double pAp   = dot_product(p, Ap, J, K);
        double alpha = r_dot_r / pAp;

        // x = x + alpha*p
        // r = r - alpha*Ap
        for (long k = 2; k < K+2; k++)
            for (long j = 2; j < J+2; j++)
            {
                x(j,k,4) += alpha * p(j,k,4);
                x(j,k,5) += alpha * p(j,k,5);
                x(j,k,6) += alpha * p(j,k,6);
                r(j,k,4) -= alpha * Ap(j,k,4);
                r(j,k,5) -= alpha * Ap(j,k,5);
                r(j,k,6) -= alpha * Ap(j,k,6);
            }

        double r_dot_r_new = dot_product(r, r, J, K);
        double beta = r_dot_r_new / r_dot_r;
        r_dot_r = r_dot_r_new;

        // p = r + beta*p
        for (long k = 2; k < K+2; k++)
            for (long j = 2; j < J+2; j++)
            {
                p(j,k,4) = r(j,k,4) + beta*p(j,k,4);
                p(j,k,5) = r(j,k,5) + beta*p(j,k,5);
                p(j,k,6) = r(j,k,6) + beta*p(j,k,6);
            }
    }

    // Update B in un with solution x
    for (long k = 2; k < K+2; k++)
        for (long j = 2; j < J+2; j++)
        {
            double B2_old = un(j,k,4)*un(j,k,4) + un(j,k,5)*un(j,k,5) + un(j,k,6)*un(j,k,6);
            double B2_new = x(j,k,4)*x(j,k,4) + x(j,k,5)*x(j,k,5) + x(j,k,6)*x(j,k,6);
            un(j,k,4) = x(j,k,4);
            un(j,k,5) = x(j,k,5);
            un(j,k,6) = x(j,k,6);
            un(j,k,7) += 0.5*(B2_new - B2_old);  // update magnetic energy
        }
}