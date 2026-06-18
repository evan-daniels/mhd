////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2005 Jorge Balbas and Eitan Tadmor
//
// This program is distributed freely for research and instructional use only. 
// You may copy, modify,  and use this software for these non-commercial
// purposes, provided that the copyright notice and associated text is 
// reproduced on all copies
//
// spctral_radii.cc -- an auxiliary function of
//
// Function called by time_step_2d.cc, Hx_flux.cc and Hz_flux.cc (model 
// specific)
//
// This function calculates the spectral radii of the Jacobian matrices of
// f(u) and g(u) for ideal MHD equations
// 
// Input (passed by reference):
//
// (1) u -- a doublearray1d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy over the entire soltion domain
//          at a cell of the solution domain
//
// (2) parameters -- a doublearray1d type variable holding additional parameters 
//					for the specific example
//
// Output (returned by reference):
//
// (1) rx -- a double type variable holding the spectral radius of the 
//           Jacobian matrix of f(u)
//
// (2) rz -- a double type variable holding the spectral radius of the 
//           Jacobian matrix of g(u)
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::spectral_radii(const doublearray1d& u, const doublearray1d& parameters, double& rx, double& ry)
{
	double rho, vx, vy, vz, p, A, B, cfx, cfy;
	double gamma = parameters(0);
	
	rho = u(0);
	vx = u(1)/rho;
	vy = u(2)/rho;
	vz = u(3)/rho;
	p = (gamma - 1.0)*(u(7) - 0.5*rho*(pow(vx,2.0) + pow(vy,2.0) + pow(vz,2.0))  - 0.5*(pow(u(4),2.0) + pow(u(5),2.0) + pow(u(6),2.0)));
	p = std::max(p, 1e-10);
	A = gamma*p/rho;
	B = (pow(u(4),2.0) + pow(u(5),2.0) + pow(u(6),2.0))/rho;
	cfx = sqrt(0.5*(A + B + sqrt(pow(A+B,2.0) - 4.0*A*pow(u(4),2.0)/rho)));
	cfy = sqrt(0.5*(A + B + sqrt(pow(A+B,2.0) - 4.0*A*pow(u(6),2.0)/rho)));
	rx = fabs(vx) + cfx;
	ry = fabs(vy) + cfy;
}
