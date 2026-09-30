////////////////////////////////////////////////////////////////////////////////
//
// C_flux_em.cc — EM flux computation for the CT sub-stepper
//
// Computes:
//   (a) Upwinded Ez at x-faces and y-faces (for Bx/By CT update)
//   (b) Corner Ez from 4-point average (drives Faraday's law)
//   (c) Cx_em, Cy_em: Faraday increments for Bx and By face arrays
//   (d) Cz_em: FV increment for cell-centered Bz (using upwinded Ey, Ex)
//
// Upwinding speed is the whistler phase speed: di * |B| / (rho * dx)
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::C_flux_em(
    const doublearray2d& rho_N, const doublearray2d& rho_S,
    const doublearray2d& rho_E, const doublearray2d& rho_W,
    const doublearray3d& B_N, const doublearray3d& B_S,
    const doublearray3d& B_E, const doublearray3d& B_W,
    const doublearray3d& E_N, const doublearray3d& E_S,
    const doublearray3d& E_E, const doublearray3d& E_W,
    const doublearray1d& lambda_em, const doublearray1d& mu_em,
    const doublearray1d& dx_interface, const doublearray1d& dy_interface,
    const double& di,
    doublearray2d& Cx_em, doublearray2d& Cy_em, doublearray2d& Cz_em)
{
    long j, k, J, K;
    J = B_N.getIndex1Size() - 4;
    K = B_N.getIndex2Size() - 4;

    double Bl, Br, rl, rr, al, ar, ax;
    double Bb, Bt, rb, rt, ab, at_v, ay;

    // Internal face/corner arrays
    doublearray2d Ez_xf(J+4, K+4);    // upwinded Ez at x-faces
    doublearray2d Ez_yf(J+4, K+4);    // upwinded Ez at y-faces
    doublearray2d Fx_xf(J+4, K+4);    // Ey x-flux for Bz update
    doublearray2d Gy_yf(J+4, K+4);    // -Ex y-flux for Bz update
    doublearray2d Ez_corner(J+4, K+4); // corner Ez

    // ===== (a) Upwinded Ez and Ey at x-faces =====
    // x-face at (j+1/2, k): Left = E_E(j,k), Right = E_W(j+1,k)
    for (k = 1; k < K+3; k++)
    {
        for (j = 1; j < J+2; j++)
        {
            Bl = sqrt(B_E(j,k,0)*B_E(j,k,0) + B_E(j,k,1)*B_E(j,k,1) + B_E(j,k,2)*B_E(j,k,2));
            Br = sqrt(B_W(j+1,k,0)*B_W(j+1,k,0) + B_W(j+1,k,1)*B_W(j+1,k,1) + B_W(j+1,k,2)*B_W(j+1,k,2));
            rl = std::max(rho_E(j,k), 0.1);
            rr = std::max(rho_W(j+1,k), 0.1);
            al = di * Bl / (rl * dx_interface(j));
            ar = di * Br / (rr * dx_interface(j));
            ax = std::max(al, ar);

            // Ez: dBy/dt = +dEz/dx, so flux of By is -Ez
            // LxF gives: Ez = avg + 0.5*a*(By_R - By_L)
            Ez_xf(j,k) = 0.5*(E_E(j,k,2) + E_W(j+1,k,2))
                        + 0.5*ax*(B_W(j+1,k,1) - B_E(j,k,1));

            // Ey: dBz/dt = -dEy/dx + ..., so x-flux of Bz is F = Ey
            // LxF gives: F = avg - 0.5*a*(Bz_R - Bz_L)
            Fx_xf(j,k) = 0.5*(E_E(j,k,1) + E_W(j+1,k,1))
                        - 0.5*ax*(B_W(j+1,k,2) - B_E(j,k,2));
        }
    }

    // ===== (b) Upwinded Ez and -Ex at y-faces =====
    // y-face at (j, k+1/2): Bottom = E_N(j,k), Top = E_S(j,k+1)
    for (k = 1; k < K+2; k++)
    {
        for (j = 1; j < J+3; j++)
        {
            Bb = sqrt(B_N(j,k,0)*B_N(j,k,0) + B_N(j,k,1)*B_N(j,k,1) + B_N(j,k,2)*B_N(j,k,2));
            Bt = sqrt(B_S(j,k+1,0)*B_S(j,k+1,0) + B_S(j,k+1,1)*B_S(j,k+1,1) + B_S(j,k+1,2)*B_S(j,k+1,2));
            rb = std::max(rho_N(j,k), 0.1);
            rt = std::max(rho_S(j,k+1), 0.1);
            ab = di * Bb / (rb * dy_interface(k));
            at_v = di * Bt / (rt * dy_interface(k));
            ay = std::max(ab, at_v);

            // Ez: dBx/dt = -dEz/dy, so flux of Bx is +Ez
            // LxF gives: Ez = avg - 0.5*a*(Bx_T - Bx_B)
            Ez_yf(j,k) = 0.5*(E_N(j,k,2) + E_S(j,k+1,2))
                        - 0.5*ay*(B_S(j,k+1,0) - B_N(j,k,0));

            // -Ex: dBz/dt = ... + dEx/dy, so y-flux of Bz is G = -Ex
            // LxF gives: G = avg(-Ex) - 0.5*a*(Bz_T - Bz_B)
            Gy_yf(j,k) = -0.5*(E_N(j,k,0) + E_S(j,k+1,0))
                         - 0.5*ay*(B_S(j,k+1,2) - B_N(j,k,2));
        }
    }

    // ===== (c) Corner Ez: 4-point average of surrounding face values =====
    for (k = 1; k < K+2; k++)
        for (j = 1; j < J+2; j++)
            Ez_corner(j,k) = 0.25*(Ez_xf(j,k) + Ez_xf(j,k+1)
                                  + Ez_yf(j,k) + Ez_yf(j+1,k));

    // ===== (d) Cx_em: Bx face increment, dBx/dt = -dEz/dy =====
    for (k = 2; k < K+2; k++)
        for (j = 2; j < J+2; j++)
            Cx_em(j,k) = -mu_em(k) * (Ez_corner(j,k) - Ez_corner(j,k-1));

    // ===== (e) Cy_em: By face increment, dBy/dt = +dEz/dx =====
    for (k = 2; k < K+2; k++)
        for (j = 2; j < J+2; j++)
            Cy_em(j,k) = lambda_em(j) * (Ez_corner(j,k) - Ez_corner(j-1,k));

    // ===== (f) Cz_em: Bz cell-center increment =====
    // dBz/dt + d(Ey)/dx + d(-Ex)/dy = 0
    // Bz_new = Bz_old - lambda*(Fx_right - Fx_left) - mu*(Gy_top - Gy_bottom)
    for (k = 2; k < K+2; k++)
        for (j = 2; j < J+2; j++)
            Cz_em(j,k) = -lambda_em(j) * (Fx_xf(j,k) - Fx_xf(j-1,k))
                         - mu_em(k) * (Gy_yf(j,k) - Gy_yf(j,k-1));
}