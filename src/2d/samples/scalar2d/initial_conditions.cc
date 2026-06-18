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
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.
//
// initial_conditions.cc -- an auxiliary function of
//
// Function called by cenatpack_main, problem specific.
//
// This function initializes the solution array with the initial conditions of 
// the problem to be solved
// 
// Input (passed by reference):
//
// (1) parameters -- a doublearray1d type variable holding additional parameters 
//					for the specific example
//
// (2) dx_cell -- a doublearray1d type variable holding the x-scale of each cell
//
// (3) dy_cell -- a doublearray1d type variable holding the y-scale of each cell
//
// (4) x_cell -- a doublearray1d type variable holding the x-coordinate of the 
//				 each cell's center
//
// (5) y_cell -- a doublearray1d type variable holding the y-coordinate of the 
//				 each cell's center
//
// (6) dx_interface -- distance between consecutive cell interfaces in 
//					  x-direction
//
// (7) dy_interface -- distance between consecutive cell interfaces in 
//					  y-direction
//
// (8) x -- a doublearray1d type variable holding the x-coordinate of the 
//		    cell interfaces
//
// (9) y -- a doublearray1d type variable holding the y-coordinate of the 
//		    cell interfaces
//
// Output (returned by reference):
//
// (1) un -- a doublearray3d type variable holding the conserved
//           quantities density, momentum (3 components), magnetic field
//           (3 components), and total energy over the entire soltion domain
//           at t=0
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::initial_conditions(doublearray3d& un, const doublearray1d& parameters, const doublearray1d& dx_cell, const doublearray1d& dy_cell, const doublearray1d& dx_interface, const doublearray1d& dy_interface, const doublearray1d& x_cell, const doublearray1d& y_cell, const doublearray1d& x, const doublearray1d& y)
{
	
	long J = un.getIndex1Size() - 4;
	long K = un.getIndex2Size() - 4;
	long L = un.getIndex3Size();
	long j, k, l;
	
	for (k = 0; k < K + 4; k++)
	{
		for (j = 0; j < J + 4; j++)
		{
			un(j,k,0) = sin(x_cell(j) + 0.5)*cos(2.0*x_cell(j) + y_cell(k));
		}
	}
}
