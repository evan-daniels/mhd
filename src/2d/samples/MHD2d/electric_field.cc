////////////////////////////////////////////////////////////////////////////////
//
// electric_field.cc — 3-component Ohm's law at all four interface positions
//
// E = -v×B + ηJ + (di/ρ)J×B
//
// Uses frozen fluid (rho, v) from the fluid reconstruction and
// EM-reconstructed (B, J) from reconstruction_em
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::electric_field(
    const doublearray2d& rho_N, const doublearray2d& rho_S,
    const doublearray2d& rho_E, const doublearray2d& rho_W,
    const doublearray3d& v_N, const doublearray3d& v_S,
    const doublearray3d& v_E, const doublearray3d& v_W,
    const doublearray3d& B_N, const doublearray3d& B_S,
    const doublearray3d& B_E, const doublearray3d& B_W,
    const doublearray3d& J_N, const doublearray3d& J_S,
    const doublearray3d& J_E, const doublearray3d& J_W,
    doublearray3d& E_N, doublearray3d& E_S,
    doublearray3d& E_E, doublearray3d& E_W,
    const double& di, const double& eta)
{
    long J = B_N.getIndex1Size() - 4;
    long K = B_N.getIndex2Size() - 4;
    long j, k;

    double rho, vx, vy, vz, Bx, By, Bz, Jx, Jy, Jz, rho_eff;

    for (k = 1; k < K+3; k++)
    {
        for (j = 1; j < J+3; j++)
        {
            // --- North interface (j, k+1/2) ---
            rho = rho_N(j,k);
            rho_eff = std::max(rho, 0.1);
            vx = v_N(j,k,0);  vy = v_N(j,k,1);  vz = v_N(j,k,2);
            Bx = B_N(j,k,0);  By = B_N(j,k,1);  Bz = B_N(j,k,2);
            Jx = J_N(j,k,0);  Jy = J_N(j,k,1);  Jz = J_N(j,k,2);

            E_N(j,k,0) = -(vy*Bz - vz*By) + eta*Jx + di*(Jy*Bz - Jz*By)/rho_eff;
            E_N(j,k,1) =  (vx*Bz - vz*Bx) + eta*Jy - di*(Jx*Bz - Jz*Bx)/rho_eff;
            E_N(j,k,2) = -(vx*By - vy*Bx) + eta*Jz + di*(Jx*By - Jy*Bx)/rho_eff;

            // --- South interface (j, k-1/2) ---
            rho = rho_S(j,k);
            rho_eff = std::max(rho, 0.1);
            vx = v_S(j,k,0);  vy = v_S(j,k,1);  vz = v_S(j,k,2);
            Bx = B_S(j,k,0);  By = B_S(j,k,1);  Bz = B_S(j,k,2);
            Jx = J_S(j,k,0);  Jy = J_S(j,k,1);  Jz = J_S(j,k,2);

            E_S(j,k,0) = -(vy*Bz - vz*By) + eta*Jx + di*(Jy*Bz - Jz*By)/rho_eff;
            E_S(j,k,1) =  (vx*Bz - vz*Bx) + eta*Jy - di*(Jx*Bz - Jz*Bx)/rho_eff;
            E_S(j,k,2) = -(vx*By - vy*Bx) + eta*Jz + di*(Jx*By - Jy*Bx)/rho_eff;

            // --- East interface (j+1/2, k) ---
            rho = rho_E(j,k);
            rho_eff = std::max(rho, 0.1);
            vx = v_E(j,k,0);  vy = v_E(j,k,1);  vz = v_E(j,k,2);
            Bx = B_E(j,k,0);  By = B_E(j,k,1);  Bz = B_E(j,k,2);
            Jx = J_E(j,k,0);  Jy = J_E(j,k,1);  Jz = J_E(j,k,2);

            E_E(j,k,0) = -(vy*Bz - vz*By) + eta*Jx + di*(Jy*Bz - Jz*By)/rho_eff;
            E_E(j,k,1) =  (vx*Bz - vz*Bx) + eta*Jy - di*(Jx*Bz - Jz*Bx)/rho_eff;
            E_E(j,k,2) = -(vx*By - vy*Bx) + eta*Jz + di*(Jx*By - Jy*Bx)/rho_eff;

            // --- West interface (j-1/2, k) ---
            rho = rho_W(j,k);
            rho_eff = std::max(rho, 0.1);
            vx = v_W(j,k,0);  vy = v_W(j,k,1);  vz = v_W(j,k,2);
            Bx = B_W(j,k,0);  By = B_W(j,k,1);  Bz = B_W(j,k,2);
            Jx = J_W(j,k,0);  Jy = J_W(j,k,1);  Jz = J_W(j,k,2);

            E_W(j,k,0) = -(vy*Bz - vz*By) + eta*Jx + di*(Jy*Bz - Jz*By)/rho_eff;
            E_W(j,k,1) =  (vx*Bz - vz*Bx) + eta*Jy - di*(Jx*Bz - Jz*Bx)/rho_eff;
            E_W(j,k,2) = -(vx*By - vy*Bx) + eta*Jz + di*(Jx*By - Jy*Bx)/rho_eff;
        }
    }
}