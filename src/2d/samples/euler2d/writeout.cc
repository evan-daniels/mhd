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
// write_out.cc -- an auxiliary function of
//
// Function called by solver.cc (main)
//
// This function outputs files containing the conserved solution variables (or
// corresponding primitive variables) for the ideal MHD model.
//
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t over the discretized solution domain
//
// (2) parameters -- a doublearray1d type variable holding additional parameters
//					for the specific example
//
// (3) t -- time of solution being output
//
// (4) n -- output number
//
// (5) id -- process id number
//
// (6) p -- number of processors for the simulation
//
// Output (returned by reference):
//
// (1) Nine text files holding the solution variables at time t (n keeps track
//     of how many output calls have been executed over the time of the
//     simulation
//
//     rho_id_n -- density
//
//     u1_id_n -- velocity in x-direction
//
//     u2_id_n -- velocity in y-direction
//
//     p_id_n -- pressure
//
//     t_n -- time
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::writeout(const doublearray3d& un, const double& t, const doublearray1d& parameters, const long& n, const int& id, const int& p)
{

	cout.setf(ios::scientific, ios::floatfield);

	long j, k, J, K;
	char rho_file[30];
	char u1_file[30];
	char u2_file[30];
	char p_file[30];
	char t_file[30];

	double gamma = parameters(0);

	J = un.getIndex1Size() - 4;
	K = un.getIndex2Size() - 4;

	doublearray2d rho(J,K);
	doublearray2d u1(J,K);
	doublearray2d u2(J,K);
	doublearray2d press(J,K);

	for (j = 0; j < J; j++)
	{
	  for (k = 0; k < K; k++)
	  {
		  rho(j,k) = un(j+2,k+2,0);
		  u1(j,k) = un(j+2,k+2,1)/un(j+2,k+2,0);
		  u2(j,k) = un(j+2,k+2,2)/un(j+2,k+2,0);
		  press(j,k) = (gamma - 1.0)*(un(j+2,k+2,3)-.5*(pow(un(j+2,k+2,1),2.0) + pow(un(j+2,k+2,2),2.0))/un(j+2,k+2,0));
	  }
	}

	sprintf(rho_file, "rho_files/rho_p%0.3d_%0.3ld", id, n);
	sprintf(u1_file, "u1_files/u1_p%0.3d_%0.3ld", id, n);
	sprintf(u2_file, "u2_files/u2_p%0.3d_%0.3ld", id, n);
	sprintf(p_file, "p_files/p_p%0.3d_%0.3ld", id, n);

	ofstream OutFile;
	OutFile.open(rho_file, ios::out);
	OutFile<< rho;
	OutFile.close();

	OutFile.open(u1_file, ios::out);
	OutFile<< u1;
	OutFile.close();

	OutFile.open(u2_file, ios::out);
	OutFile<< u2;
	OutFile.close();

	OutFile.open(p_file, ios::out);
	OutFile<< press;
	OutFile.close();

	if (id == 0)
	{
	  sprintf(t_file, "t_files/t_%0.3ld", n);

	  OutFile.open(t_file, ios::out);
	  OutFile<< t;
	  OutFile.close();
	}
}

void CENTPACK::writeout(const doublearray3d& un, const double& t, const doublearray1d& parameters, const long& n, const bool& odd, const int& id, const int& p)
{

	cout.setf(ios::scientific, ios::floatfield);

	long j, k, J, K;
	char rho_file[30];
	char u1_file[30];
	char u2_file[30];
	char p_file[30];
	char t_file[30];
	char odd_file[30];

	double gamma = parameters(0);

	J = un.getIndex1Size() - 4;
	K = un.getIndex2Size() - 4;

	doublearray2d rho(J,K);
	doublearray2d u1(J,K);
	doublearray2d u2(J,K);
	doublearray2d press(J,K);

	for (j = 0; j < J; j++)
	{
	  for (k = 0; k < K; k++)
	  {
		  rho(j,k) = un(j+2,k+2,0);
		  u1(j,k) = un(j+2,k+2,1)/un(j+2,k+2,0);
		  u2(j,k) = un(j+2,k+2,2)/un(j+2,k+2,0);
		  press(j,k) = (gamma - 1.0)*(un(j+2,k+2,3)-.5*(pow(un(j+2,k+2,1),2.0) + pow(un(j+2,k+2,2),2.0))/un(j+2,k+2,0));
	  }
	}

	sprintf(rho_file, "rho_files/rho_p%0.3d_%0.3ld", id, n);
	sprintf(u1_file, "u1_files/u1_p%0.3d_%0.3ld", id, n);
	sprintf(u2_file, "u2_files/u2_p%0.3d_%0.3ld", id, n);
	sprintf(p_file, "p_files/p_p%0.3d_%0.3ld", id, n);

	ofstream OutFile;
	OutFile.open(rho_file, ios::out);
	OutFile<< rho;
	OutFile.close();

	OutFile.open(u1_file, ios::out);
	OutFile<< u1;
	OutFile.close();

	OutFile.open(u2_file, ios::out);
	OutFile<< u2;
	OutFile.close();

	OutFile.open(p_file, ios::out);
	OutFile<< press;
	OutFile.close();

	if (id == 0)
	{
		sprintf(t_file, "t_files/t_%0.3ld", n);
		sprintf(odd_file, "odd_files/odd_%0.3ld", n);

		OutFile.open(t_file, ios::out);
		OutFile<< t;
		OutFile.close();

		OutFile.open(odd_file, ios::out);
		OutFile<<odd;
		OutFile.close();
	}
}
