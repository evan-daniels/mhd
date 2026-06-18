////////////////////////////////////////////////////////////////////////////
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
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Hx_flux_2d_SD2.cc -- a core function of CentPack
//
// Function called by C_flux
//
// Requires: flux_x and spectral_radii (see prototypes below)
//
// This function calculates the numerical flux of the semi-discre formulation 
// for the evolution of the cell averages of u at each side of the cell 
// interfaces x=x_j-1/2 and x=x_j+1/2
//
// Input (passed by reference):
//
// (1) two doublearray1d type variables
//
//     u_e  --  holds the point values of u_E for the cell centered at 
//              (x_j-1,y_k) or (x_j,y_k)
//
//     u_w  --  holds the point values of u_W for the cell centered at
//              (x_j, y_k) or (x_j+1,y_k)
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) Hx -- numerical flux at cell interface x=x_j-1/2 or x=x_j+1/2
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::Hx_flux_2d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, 
  const doublearray1d& parameters, const doublearray1d& e_here, doublearray1d& Hx)
{
  long l;
  long L = Hx.getIndex1Size();
  double rw, re, rn, rs, a;

  doublearray1d f_e(L), f_w(L);
  
  spectral_radii(u_e, parameters, re, rn);
  spectral_radii(u_w, parameters, rw, rs);

  if (re > rw)
	  a = re;
  else
	  a = rw;
  
  flux_x(u_e, parameters, e_here, f_e);
  flux_x(u_w, parameters, e_here,  f_w);
  
  for (l = 0; l < L; l++)
	  Hx(l) = 0.5*((f_w(l) + f_e(l)) - a*(u_w(l) - u_e(l)));
}
