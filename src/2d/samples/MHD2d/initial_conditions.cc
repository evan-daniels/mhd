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
  
  double I1, I2, I3, I4;
  double gamma = parameters(0);
  int problem = (int) parameters(1); // 0->OT, 1->Reconnection

  if (problem == 0)
  {
  
	for (k = 0; k < K+4; k++)
	{
		for (j = 0; j < J+4; j++)
		{
			un(j,k,0) = pow(gamma,2.0);
			un(j,k,1) = un(j,k,0)*(1.0/dy_cell(k))*(cos(y(k+1)) - cos(y(k)));
			un(j,k,2) = 0.0;
			un(j,k,3) = -un(j,k,0)*(1.0/dx_cell(j))*(cos(x(j+1)) - cos(x(j)));
			un(j,k,4) = (1.0/dy_cell(k))*(cos(y(k+1)) - cos(y(k)));
			un(j,k,5) = 0.0;
			un(j,k,6) = (-1.0/(2.0*dx_cell(j)))*(cos(2.0*x(j+1)) - cos(2.0*x(j)));
			
			I1 = -(0.125/dy_cell(k))*(un(j,k,0) + 1.0)*(sin(2.0*y(k+1)) - sin(2.0*y(k)));
			I2 = -(0.125/dx_cell(j))*un(j,k,0)*(sin(2.0*x(j+1)) - sin(2.0*x(j)));
			I3 = -(0.0625/dx_cell(j))*(sin(4.0*x(j+1)) - sin(4.0*x(j)));
			un(j,k,7) = 3.0 + 0.5*un(j,k,0) + I1 + I2 + I3;
		}
	}
}

else if (problem == 1)
{
    double B0     = 1.0;
    double n0     = 1.0;
    double nb = 0.5; // was .2
    double lambda = 0.5;
    double psi0   = 0.1;
    double Lx     = 19.6;
    double Ly     = 19.6;

    double Bz = 1e-4 * cos(2.0*M_PI*xc/Lx) * sin(M_PI*yc/Ly);  // small seed

    for (k = 0; k < K+4; k++)
    {
      for (j = 0; j < J+4; j++)
      {
        double xc = x_cell(j);
        double yc = y_cell(k);

        double rho = n0 / pow(cosh(yc / lambda), 2.0) + nb;
        double p   = 0.5 * rho;
        double vx  = 0.0, vy = 0.0, vz = 0.0;

        double Bx = B0 * tanh(yc / lambda)
                    - psi0 * (M_PI / Ly) * cos(2.0*M_PI*xc / Lx) * sin(M_PI*yc / Ly);
        double By = psi0 * (2.0*M_PI / Lx) * sin(2.0*M_PI*xc / Lx) * cos(M_PI*yc / Ly);
        double Bz = 0.0;

        double U = p / (gamma - 1.0)
                 + 0.5 * rho * (vx*vx + vy*vy + vz*vz)
                 + 0.5 * (Bx*Bx + By*By + Bz*Bz);

        un(j,k,0) = rho;
        un(j,k,1) = rho * vx;
        un(j,k,2) = rho * vy;
        un(j,k,3) = rho * vz;
        un(j,k,4) = Bx;
        un(j,k,5) = By;
        un(j,k,6) = Bz;
        un(j,k,7) = U;
      }
    }
}
}