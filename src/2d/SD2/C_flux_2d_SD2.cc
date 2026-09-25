#include "centpack_2d_SD2.h"

using namespace std;

// Function made to compute Ez at a point from a reconstructed state
static double Ez_from_state(const doublearray1d& u, const doublearray1d& j_here,
                              double di, double eta)
{
    double rho = u(0);
    double vx  = u(1) / rho;
    double vy  = u(2) / rho;
    double Bx  = u(4);
    double By  = u(5);
    double jz  = j_here(2);
    double rho_eff = std::max(rho, 0.1);

    double Ez_ideal = -(vx * By - vy * Bx);
    double Ez_hall  = di * (j_here(0) * By - j_here(1) * Bx) / rho_eff;
    double Ez_res   = eta * jz;
    // return Ez_ideal + Ez_res;
    return Ez_ideal + Ez_res + Ez_hall;
}

void CENTPACK::C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& jcurl, const doublearray3d& B1_xf, const doublearray3d& B2_yf, doublearray3d& C, doublearray3d& Ez_corner)
{
    long j, k, l, J, K, L;

    J = C.getIndex1Size() - 4;
    K = C.getIndex2Size() - 4;
    L = C.getIndex3Size();

    double eta = parameters(2);
    double di  = parameters(3);

    doublearray3d Ez_xf(J+4, K+4, 1);
    doublearray3d Ez_yf(J+4, K+4, 1);

    doublearray1d u_nkm1(L), u_n(L), u_s(L), u_skp1(L);
    doublearray1d u_ejm1(L), u_e(L), u_w(L), u_wjp1(L);
    doublearray1d j_here(3);

    doublearray1d Hx_halfp(L), Hx_halfm(L), Hy_halfp(L), Hy_halfm(L);

    for (k = 2; k < K+2; k++)
    {
        for (j = 2; j < J+2; j++)
        {
            for (l = 0; l < L; l++)
            {
                u_nkm1(l) = u_N(j,k-1,l);
                u_n(l)    = u_N(j,k,l);
                u_s(l)    = u_S(j,k,l);
                u_skp1(l) = u_S(j,k+1,l);
                u_ejm1(l) = u_E(j-1,k,l);
                u_e(l)    = u_E(j,k,l);
                u_w(l)    = u_W(j,k,l);
                u_wjp1(l) = u_W(j+1,k,l);
            }

            // Pass local j to flux functions
            j_here(0) = jcurl(j,k,0);
            j_here(1) = jcurl(j,k,1);
            j_here(2) = jcurl(j,k,2);

            Hx_flux_2d_SD2(u_w,    u_ejm1, parameters, j_here, Hx_halfm);
            Hx_flux_2d_SD2(u_wjp1, u_e,    parameters, j_here, Hx_halfp);
            Hy_flux_2d_SD2(u_s,    u_nkm1, parameters, j_here, Hy_halfm);
            Hy_flux_2d_SD2(u_skp1, u_n,    parameters, j_here, Hy_halfp);

            for (l = 0; l < L; l++)
                C(j,k,l) = -lambda(j)*(Hx_halfp(l) - Hx_halfm(l)) - mu(k)*(Hy_halfp(l) - Hy_halfm(l));

            // This is not truly constrained transport

            // B1 and B2, change. instead of averaging from neighbouring cells, we need to be calculating B1_xf and B2_yf. Those should be evolved 
            // using constrained transport. 
            // we recover them 

        }
    }

    // Ez_xf needed for j = 1..J+1, k = 1..K+2 (corner loop reads Ez_xf(j,k) and Ez_xf(j,k+1))
    for (k = 1; k < K+3; k++)
    {
        for (j = 1; j < J+2; j++)
        {
            for (l = 0; l < L; l++)
            {
                u_e(l)    = u_E(j,k,l);
                u_wjp1(l) = u_W(j+1,k,l);
            }

            // CT: Bx at x-interface comes from the staggered variable, not reconstruction
            u_e(4)    = B1_xf(j,k,0);
            u_wjp1(4) = B1_xf(j,k,0);

            j_here(0) = 0.5*(jcurl(j,k,0) + jcurl(j+1,k,0));
            j_here(1) = 0.5*(jcurl(j,k,1) + jcurl(j+1,k,1));
            j_here(2) = 0.5*(jcurl(j,k,2) + jcurl(j+1,k,2));
            // send j+1 to EzR, j to EzL, rather than average
            double EzR = Ez_from_state(u_wjp1, j_here, di, eta); // EzL originally
            double EzL = Ez_from_state(u_e,    j_here, di, eta); // EzR originally

            // LxF upwinding on tangential B jump at interface
            double rho_avg = std::max(0.5*(u_e(0) + u_wjp1(0)), 0.1);
            double Bsq = B1_xf(j,k,0)*B1_xf(j,k,0)
                       + 0.25*(u_e(5)+u_wjp1(5))*(u_e(5)+u_wjp1(5))
                       + 0.25*(u_e(6)+u_wjp1(6))*(u_e(6)+u_wjp1(6));
            double a_x = sqrt(Bsq / rho_avg);

            double By_jump = u_wjp1(5) - u_e(5);
            Ez_xf(j, k, 0) = 0.5*(EzL + EzR) + 0.5*a_x*By_jump;

            // Ez_xf(j, k, 0) = 0.5*(EzL + EzR); // add upwinding. y magnetic field i reconstructed across the x interface east and west of the interface. take that jump, multiply by 1/2 speed of propagation
        }
    }

    // Ez_yf needed for j = 1..J+2, k = 1..K+1 (corner loop reads Ez_yf(j,k) and Ez_yf(j+1,k))
    for (k = 1; k < K+2; k++)
    {
        for (j = 1; j < J+3; j++)
        {
            for (l = 0; l < L; l++)
            {
                u_n(l)    = u_N(j,k,l);
                u_skp1(l) = u_S(j,k+1,l);
            }

            // CT: By at y-interface comes from the staggered variable, not reconstruction
            u_n(5)    = B2_yf(j,k,0);
            u_skp1(5) = B2_yf(j,k,0); // bottom of the k+1 cell

            j_here(0) = 0.5*(jcurl(j,k,0) + jcurl(j,k+1,0));
            j_here(1) = 0.5*(jcurl(j,k,1) + jcurl(j,k+1,1));
            j_here(2) = 0.5*(jcurl(j,k,2) + jcurl(j,k+1,2));

            // ^ Dr. Balbas is using the k+1 value to calculate 134, k value for 135
            // after he calculates the E fields from 134 + 135 (where i have 136 (avg)), that is where he adds the upwinding

            double EzT = Ez_from_state(u_skp1, j_here, di, eta);
            double EzB = Ez_from_state(u_n,    j_here, di, eta); // is the upper end of the south cell

            // LxF upwinding on tangential B jump at interface
            double rho_avg = std::max(0.5*(u_n(0) + u_skp1(0)), 0.1);
            double Bsq = 0.25*(u_n(4)+u_skp1(4))*(u_n(4)+u_skp1(4))
                       + B2_yf(j,k,0)*B2_yf(j,k,0)
                       + 0.25*(u_n(6)+u_skp1(6))*(u_n(6)+u_skp1(6));
            double a_y = sqrt(Bsq / rho_avg);

            double Bx_jump = u_skp1(4) - u_n(4);
            Ez_yf(j, k, 0) = 0.5*(EzB + EzT) - 0.5*a_y*Bx_jump;

            // Ez_yf(j, k, 0) = 0.5*(EzB + EzT); // x magnetic field i reconstructed below and above the interface. take that jump, multiply by 1/2 speed of propagation
        }
    }

    // 
    // Below is CT attempt
    // Ez at corners

    // for (k = 1; k < K+2; k++)
    // {
    //     for (j = 1; j < J+2; j++)
    //     {
    //         double Ez_avg = 0.25*(Ez_xf(j, k,   0)
    //                             + Ez_xf(j, k+1, 0)
    //                             + Ez_yf(j,   k, 0)
    //                             + Ez_yf(j+1,   k, 0));

    //         // Upwind dissipation

    //         double dBx_y = 0.5*((B1_xf(j,k+1,0)   - B1_xf(j,k,  0))
    //                            +(B1_xf(j+1,k+1,0) - B1_xf(j+1,k,0)));
    //         double dBy_x = 0.5*((B2_yf(j+1,k,0)   - B2_yf(j,  k,0))
    //                            +(B2_yf(j+1,k+1,0) - B2_yf(j,k+1,0)));

    //         // double Bx = B1_xf(j, k, 0);
    //         // double By = B2_yf(j, k, 0);
    //         // double B2 = Bx * Bx + By * By;
    //         // double rho_loc = std::max(0.5*(u_E(j, k, 0) + u_W(j+1, k, 0)), 0.1);
    //         // double a_loc = sqrt(B2 / rho_loc);

    //         // Ez_corner(j, k, 0) = Ez_avg + 0.5 * a_loc * (dBx_y - dBy_x);

    //         // Dissipation speed from surrounding FACE states, not the corner's own
    //         // (self-referential) B — using the corner's own B here creates a positive
    //         // feedback loop that can blow up instead of damping.
    //         double B1n = B1_xf(j,   k,   0);
    //         double B1s = B1_xf(j+1, k,   0);
    //         double B2e = B2_yf(j,   k,   0);
    //         double B2w = B2_yf(j,   k+1, 0);
    //         double B2max = std::max(std::max(B1n*B1n, B1s*B1s), std::max(B2e*B2e, B2w*B2w));
    //         double rho_loc = std::max(0.5*(u_E(j, k, 0) + u_W(j+1, k, 0)), 0.1);
    //         double a_loc = sqrt(B2max / rho_loc);

    //         // Hard cap as a defensive backstop against transient spikes
    //         const double a_loc_max = 50.0;   // tune against a healthy run's spectral_radii output
    //         a_loc = std::min(a_loc, a_loc_max);

    //         Ez_corner(j, k, 0) = Ez_avg + 0.5 * a_loc * (dBx_y - dBy_x);

    //         // mixing central scheme with constrained transport (doing both at the same time). 
    //         // for constrained transport, we're going to use the central scheme to update the averages of the fluid variables at the center of the cell
    //         // then, evolve the values of Bx only at the x interfaces of the cell. updates values for the magnetic field should be avergage of Xj+1/2, Xj-1/2. those values
    //         // for the x magnetic field, do not need to be reconstructed in the x direction. whatever value is reconstructed is value for both left and right interface
    //         // '' same for y
    //         // for constrained transport, we recover those from the Ez values in the corners. 

    //         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //         ////////////        I need to evolve magnetic field separately from the fluid variables !!!!!            /////////////
    //         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //         // separate equation for the evolution of Bx and By because we are doing it at the cell interfaces
    //         // regardless of whether we are recovering J and reconstructing E, or vice versa, we need a separate solver for the electromagnetic variables at the interfaces
    //         // from the fluid variables at the cell center.

    //         // fluid and hall physics will be evolved on different time steps. the delta t on the hall solver might potentially have the same dt as the fluid step
    //         // and it can all be done in one step

    // // Above is CT attempt
    //     }
    // }
    // 

    for (k = 1; k < K+2; k++)
    {
        for (j = 1; j < J+2; j++)
        {
            Ez_corner(j, k, 0) = 0.25*(Ez_xf(j, k,   0)
                                        + Ez_xf(j, k+1, 0)
                                        + Ez_yf(j,   k, 0)
                                        + Ez_yf(j+1, k, 0));
        }
    }

}