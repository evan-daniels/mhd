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
// read_in_1d.cc -- main function of CentPack
//
// Called by solver (main)
// 
// This function reads input parameters for the simulation provided by the user // in a separate file named "input"
//
// Input (passed by reference):
// 
// (1) x_init -- a double type variable holding the left end point of the 
//               x-interval of the solution
//
// (2) x_final -- a double type variable holding the right end point of the 
//                x-interval of the solution
// (3) J  --  a long type variable holding the number of cells along x-dimension
//
// (4) L  --  a long type variable holding the number of components in the 
//            system
//
// (5) gamma - a double type variable holding the ratio of specific heats
//
// (6) cfl  --  a double type variable holding the CFL restriction to determine 
//              time step
// 
// (7) t_final  --  a double type variable holding the time of simulation
//
// (8) dt_out  --  a double type variable holding the desired output interval
//
// (9) alpha  --  minmod limiting parameter
//
// Output:
//
// This function does not produce any direct output, it only initializes the
// input variables described above with the values provided by the user
//
////////////////////////////////////////////////////////////////////////////////
#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::read_in_1d(double& x_left, double& x_right, long& J, long& L, double& t_final, double& dt_out, double& cfl, double& alpha, doublearray1d& parameters)
{	
	long k, par_count;
	string s;
	
	ifstream InFile;
	
	InFile.open("input", ios::in);
	
	if (!InFile)
	{
		cerr<<"Error in opening input file";
		exit(1);
	}
	
	getline(InFile, s, '\n');
	InFile>> x_left;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> x_right;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> J;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> L;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> t_final;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> dt_out;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> cfl;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> alpha;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> par_count;
	
	parameters.resize(par_count);
	
	for (k = 0; k < par_count; k++)
	{	
		getline(InFile, s, '\n');
		getline(InFile, s, '\n');
		InFile>> parameters(k);
	}
}
