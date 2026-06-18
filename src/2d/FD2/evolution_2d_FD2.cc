////////////////////////////////////////////////////////////////////////////////
//
// evolution.cc -- a core function of
//
// Central Hyperbolic Solver (C) Jorge Balbas and Eitan Tadmor, Nov. 2004
//
// Function called by solver (main)
//
// Requires: boundary_conditions.cc, norm.cc, reconstruction.cc,
// reconstruction_diag.cc, indicators.cc, indicators_diag.cc and C_flux.cc
// (see prototypes below)
//
// This function evolves the cell averages of u from t=t^n to t=t^n + dt
// 
// Input (passed by reference):
//
// (1)  un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n over the discretized solution domain
//
// (2) lambda  --  a double type variable holding the mesh ratio dt/dx
//
// (3) dx  --  a double type variable holding the space scale dx
//
// (4) dz  --  a double type variable holding the space scale dz
//
// (5) mu  --  a double type variable holding the mesh ratio dt/dz
//
// (6) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1)  un  --  a doublearray3d type variable holding the cell averages of u
//              at time t=t^n + dt over the discretized solution domain
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"

using namespace std;

void CENTPACK::evolution_2d_FD2(doublearray3d& un, const doublearray3d& unhalf, const double& alpha, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const bool& odd, const int& id, const int& p)
{
	long j, k, l, J, K, L;
	
	J = un.getIndex1Size() - 4;
	K = un.getIndex2Size() - 4;
	L = un.getIndex3Size();
	
	doublearray3d u_star(J+4,K+4,L);
	
	predictor_2d_FD2(un, u_star, alpha, lambda, mu, parameters);
	corrector_2d_FD2(un, unhalf, u_star, lambda, mu, parameters, odd);
	
	boundary_conditions(un, parameters, odd, id, p);
}
