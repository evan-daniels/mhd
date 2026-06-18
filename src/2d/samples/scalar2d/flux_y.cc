//////////////////////////////////////////////////////////////////////////
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
// PARTICULAR PURPOSE
//
// flux_y.cc -- Ideal y-flux function, an auxiliary function of
//
// Function called by Hy_flux
//
// Input (passed by reference):
//
// (1) u -- a doublearray1d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy at a point of the discretized 
//          solution domain.
//
// (2) parameters -- a doublearray1d type variable holding parameters of the 
//                   physical model such as the ratio of specific heats
//
// Output (returned by reference):
//
// (1) g -- a doublearray1d type variable holding g(u), the flux in
//          the y direction.
//
//////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"

using namespace std;

void CENTPACK::flux_y(const doublearray1d& u, const doublearray1d& parameters, doublearray1d& g)
{
	g(0) = (1.0/3.0)*pow(u(0),3.0);
}
