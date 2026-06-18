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
// initial_conditions.cc
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
	long k1 = K/2 + 2;
	
	double gamma = parameters(0);
	
	double p_one = 1.5;
	double p_two = 0.3;
	double p_three = 0.029;
	double p_four = 0.3;
	
	doublearray1d u_one(L), u_two(L), u_three(L), u_four(L);
	
	u_one(0) = 1.5;
	u_one(1) = 0.0;
	u_one(2) = 0.0;
	u_one(3) = p_one/(gamma - 1.0) + .5*(pow(u_one(1),2.0) + pow(u_one(2),2.0))/u_one(0);
	
	u_two(0) = 0.5323;
	u_two(1) = 1.206*u_two(0);
	u_two(2) = 0.0;
	u_two(3) = p_two/(gamma - 1.0) + .5*(pow(u_two(1),2.0) + pow(u_two(2),2.0))/u_two(0);
	
	u_three(0) = 0.138;
	u_three(1) = 1.206*u_three(0);
	u_three(2) = 1.206*u_three(0);
	u_three(3) = p_three/(gamma - 1.0) + .5*(pow(u_three(1),2.0) + pow(u_three(2),2.0))/u_three(0);
	
	u_four(0) = 0.5323;
	u_four(1) = 0.0;
	u_four(2) = 1.206*u_four(0);
	u_four(3) = p_four/(gamma - 1.0) + .5*(pow(u_four(1),2.0) + pow(u_four(2),2.0))/u_four(0);
	
	for(l = 0; l < L; l++)
	{
	  for (k = k1; k < K+4; k++)
	  {
		  for (j = 0; j < J+4; j++)
		  {
			  if (x_cell(j) < 0.5)
				  un(j,k,l) = u_two(l);
			  else
				  un(j,k,l) = u_one(l);
		  }
	  }
	  
	  for (k = 0; k < k1; k++)
	  {
		  for (j = 0; j < J+4; j++)
		  {
			  if (x_cell(j) < 0.5)
				  un(j,k,l) = u_three(l);
			  else
				  un(j,k,l) = u_four(l);
		  }
	  }
	}
}
