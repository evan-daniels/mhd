////////////////////////////////////////////////////////////////////////////////
//
// evolution_2d_SD2.cc — interleaved fluid / EM time step
//
// Fluid (rho, rho*v, energy) advances with the SSP-RK2 central scheme.
// B advances with a constrained-transport EM solver that sub-cycles on the
// whistler CFL. The EM sub-cycles are interleaved with the fluid RK stages:
//
//   fluid stage 1 (t^n state)          -> un(fluid) += C0
//   EM sub-cycle  0        -> dt/2      (rho, v frozen at t^n interface values)
//   fluid stage 2 (stage-1 state)      -> un(fluid) += 0.5*(C1 - C0)
//   EM sub-cycle  dt/2     -> dt        (rho, v frozen at stage-1 interface values)
//
// un layout (8 components):
//   0 rho, 1-3 rho*v, 4-6 B (cell centered), 7 total energy
//   Fluid update touches 0-3 and 7 only. B is owned by the EM solver:
//   Bx, By live on faces (B1_xf, B2_yf) and are recovered to cell centers;
//   Bz lives in un(:,:,6) and is updated by a finite-volume EM flux.
//
// Face conventions:
//   B1_xf(j,k) : Bx on the east face of cell j   (x_{j+1/2}, y_k)
//   B2_yf(j,k) : By on the north face of cell k  (x_j, y_{k+1/2})
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

void CENTPACK::evolution_2d_SD2(doublearray3d& un, doublearray3d& B1_xf, doublearray3d& B2_yf,
    const doublearray1d& lambda, const doublearray1d& mu,
    const doublearray1d& dx_cell, const doublearray1d& dx_interface,
    const doublearray1d& dy_cell, const doublearray1d& dy_interface,
    const double& alpha, const doublearray1d& parameters, const double& dt_fluid,
    const int& id, const int& p)
{
    using namespace std;

    long j, k, l, J, K, L;
    J = un.getIndex1Size() - 4;
    K = un.getIndex2Size() - 4;
    L = un.getIndex3Size();

    const double eta    = parameters(2);
    const double di     = parameters(3);
    const double cfl_em = 0.4;

    // Fluid components updated by the central scheme (B = 4,5,6 excluded)
    const int n_fluid = 5;
    const int fluid_comp[n_fluid] = {0, 1, 2, 3, 7};

    // ---------- fluid arrays ----------
    doublearray3d u_N(J+4,K+4,L), u_S(J+4,K+4,L), u_E(J+4,K+4,L), u_W(J+4,K+4,L);
    doublearray3d C0(J+4,K+4,L), C1(J+4,K+4,L);

    // ---------- frozen fluid state at interfaces ----------
    doublearray2d rho_N(J+4,K+4), rho_S(J+4,K+4), rho_E(J+4,K+4), rho_W(J+4,K+4);
    doublearray3d v_N(J+4,K+4,3), v_S(J+4,K+4,3), v_E(J+4,K+4,3), v_W(J+4,K+4,3);

    // ---------- EM arrays ----------
    doublearray3d Jc(J+4,K+4,3);
    doublearray3d B_N(J+4,K+4,3), B_S(J+4,K+4,3), B_E(J+4,K+4,3), B_W(J+4,K+4,3);
    doublearray3d J_N(J+4,K+4,3), J_S(J+4,K+4,3), J_E(J+4,K+4,3), J_W(J+4,K+4,3);
    doublearray3d E_N(J+4,K+4,3), E_S(J+4,K+4,3), E_E(J+4,K+4,3), E_W(J+4,K+4,3);
    doublearray2d Cx_em(J+4,K+4), Cy_em(J+4,K+4), Cz_em(J+4,K+4);
    doublearray1d lambda_em(J+4), mu_em(K+4);
    doublearray1d u_cell(L);

    double t_em = 0.0;

    // ------------------------------------------------------------------------
    // J = curl B on the interior; ghosts filled by periodic-y / MPI-x exchange
    // ------------------------------------------------------------------------
    auto compute_J = [&]()
    {
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
            {
                double dx = dx_cell(j), dy = dy_cell(k);
                Jc(j,k,0) =  (un(j,k+1,6) - un(j,k-1,6))/(2.0*dy);
                Jc(j,k,1) = -(un(j+1,k,6) - un(j-1,k,6))/(2.0*dx);
                Jc(j,k,2) =  (un(j+1,k,5) - un(j-1,k,5))/(2.0*dx)
                           - (un(j,k+1,4) - un(j,k-1,4))/(2.0*dy);
            }
        boundary_conditions_ct(Jc, id, p);
    };

    // ------------------------------------------------------------------------
    // Pull rho and v from the fluid reconstruction; frozen during a sub-cycle
    // ------------------------------------------------------------------------
    auto freeze_fluid = [&]()
    {
        for (k = 1; k < K+3; k++)
            for (j = 1; j < J+3; j++)
            {
                rho_N(j,k) = u_N(j,k,0);  rho_S(j,k) = u_S(j,k,0);
                rho_E(j,k) = u_E(j,k,0);  rho_W(j,k) = u_W(j,k,0);
                for (l = 0; l < 3; l++)
                {
                    v_N(j,k,l) = u_N(j,k,l+1)/u_N(j,k,0);
                    v_S(j,k,l) = u_S(j,k,l+1)/u_S(j,k,0);
                    v_E(j,k,l) = u_E(j,k,l+1)/u_E(j,k,0);
                    v_W(j,k,l) = u_W(j,k,l+1)/u_W(j,k,0);
                }
            }
    };

    // ------------------------------------------------------------------------
    // EM sub-cycle: forward Euler steps from the current t_em up to t_target
    // ------------------------------------------------------------------------
    auto em_subcycle = [&](double t_target)
    {
        int nsub = 0;
        while (t_em < t_target)
        {
            // --- whistler CFL, per cell, global min ---
            double dtp_em = 1.0e30;
            double rx, ry;
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    for (l = 0; l < L; l++) u_cell(l) = un(j,k,l);
                    spectral_radii_em(u_cell, dx_cell(j), dy_cell(k), di, rx, ry);
                    double rate = rx/dx_cell(j) + ry/dy_cell(k);
                    if (rate > 0.0)
                        dtp_em = std::min(dtp_em, cfl_em/rate);
                }

            if (std::isnan(dtp_em) || std::isinf(dtp_em))
                MPI_Abort(MPI_COMM_WORLD, 1);

            double dt_em;
            MPI_Allreduce(&dtp_em, &dt_em, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            dt_em = std::min(dt_em, t_target - t_em);   // land exactly on target

            for (j = 2; j < J+2; j++) lambda_em(j) = dt_em/dx_cell(j);
            for (k = 2; k < K+2; k++) mu_em(k)     = dt_em/dy_cell(k);

            // --- E from current B, J and frozen rho, v ---
            compute_J();
            reconstruction_em(un, Jc, B1_xf, B2_yf, B_N, B_S, B_E, B_W,
                              J_N, J_S, J_E, J_W,
                              dx_cell, dx_interface, dy_cell, dy_interface, alpha);
            electric_field(rho_N, rho_S, rho_E, rho_W, v_N, v_S, v_E, v_W,
                           B_N, B_S, B_E, B_W, J_N, J_S, J_E, J_W,
                           E_N, E_S, E_E, E_W, di, eta);
            C_flux_em(rho_N, rho_S, rho_E, rho_W, B_N, B_S, B_E, B_W,
                      E_N, E_S, E_E, E_W, lambda_em, mu_em,
                      dx_interface, dy_interface, di, Cx_em, Cy_em, Cz_em);

            // --- Faraday update: faces for Bx, By; cell centers for Bz ---
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    B1_xf(j,k,0) += Cx_em(j,k);
                    B2_yf(j,k,0) += Cy_em(j,k);
                    un(j,k,6)    += Cz_em(j,k);
                }

            boundary_conditions_ct(B1_xf, id, p);
            boundary_conditions_ct(B2_yf, id, p);

            // --- recover cell-centered Bx, By from faces ---
            for (k = 2; k < K+2; k++)
                for (j = 2; j < J+2; j++)
                {
                    un(j,k,4) = 0.5*(B1_xf(j,k,0) + B1_xf(j-1,k,0));
                    un(j,k,5) = 0.5*(B2_yf(j,k,0) + B2_yf(j,k-1,0));
                }

            boundary_conditions(un, parameters, id, p);

            nsub++;
            t_em += dt_em;
        }
        
        if (id == 0) printf("  EM sub-steps: %d\n", nsub);

    };
    
    // ========================================================================
    // Fluid RK stage 1 (state at t^n)
    // ========================================================================
    compute_J();
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    freeze_fluid();
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C0);

    for (int n = 0; n < n_fluid; n++)
    {
        l = fluid_comp[n];
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
                un(j,k,l) += C0(j,k,l);
    }
    boundary_conditions(un, parameters, id, p);

    // ========================================================================
    // EM sub-cycle 1: 0 -> dt/2, rho and v frozen at t^n interface values
    // ========================================================================
    em_subcycle(0.5*dt_fluid);

    // ========================================================================
    // Fluid RK stage 2 (stage-1 fluid, B at t^n + dt/2)
    // ========================================================================
    compute_J();
    reconstruction_2d_SD2(un, u_N, u_S, u_E, u_W, dx_cell, dx_interface, dy_cell, dy_interface, alpha);
    freeze_fluid();
    C_flux_2d_SD2(u_N, u_S, u_E, u_W, lambda, mu, parameters, Jc, C1);

    for (int n = 0; n < n_fluid; n++)
    {
        l = fluid_comp[n];
        for (k = 2; k < K+2; k++)
            for (j = 2; j < J+2; j++)
                un(j,k,l) += 0.5*(C1(j,k,l) - C0(j,k,l));
    }
    boundary_conditions(un, parameters, id, p);

    // ========================================================================
    // EM sub-cycle 2: dt/2 -> dt, rho and v frozen at stage-1 interface values
    // ========================================================================
    em_subcycle(dt_fluid);
}