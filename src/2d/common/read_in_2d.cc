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
// read_in_2d.cc -- main function of CentPack
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
//
// (3) y_init -- a double type variable holding the left end point of the 
//               y-interval of the solution
//
// (4) y_final -- a double type variable holding the right end point of the 
//                y-interval of the solution
//
// (5) J  --  a long type variable holding the number of cells along x-dimension
//
// (6) K  --  a long type variable holding the number of cells along y-dimension
//
// (7) L  --  a long type variable holding the number of components in the 
//            system
//
// (8) gamma - a double type variable holding the ratio of specific heats
//
// (9) cfl  --  a double type variable holding the CFL restriction to determine 
//              time step
// 
// (10) t_final  --  a double type variable holding the time of simulation
//
// (11) dt_out  --  a double type variable holding the desired output interval
//
// (12) alpha  --  minmod limiting parameter
//
// Output:
//
// This function does not produce any direct output, it only initializes the
// input variables described above with the values provided by the user
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::read_in_2d(double& x_init, double& x_final, double& y_init, double& y_final, long& J, long& K, long& L, double& t_final, double& dt_out, double& cfl, double& alpha, doublearray1d& parameters)
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
	InFile>> x_init;
	
	getline(InFile, s, '\n');
	getline(InFile, s, '\n');
	InFile>> x_final;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> y_init;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> y_final;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> J;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> K;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> L;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> t_final;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> dt_out;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> cfl;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> alpha;
	
	getline(InFile, s);
	getline(InFile, s);
	InFile>> par_count;
	
	if (par_count != 0)
	{
		parameters.resize(par_count);
	
		for (k = 0; k < par_count; k++)
		{
			getline(InFile, s);
			getline(InFile, s);
			InFile>> parameters(k);
		}
	}
	
	cout<< "x_init = " << x_init << endl;
	cout<< "x_final = " << x_final << endl;
	cout<< "y_init = " << y_init << endl;
	cout<< "y_final = " << y_final << endl;
	cout<< "J = " << J << endl;
	cout<< "K = " << K << endl;
	cout<< "L = " << L << endl;
	cout<< "t_inal = " << t_final << endl;
	cout<< "dt_out = " << dt_out << endl;
	cout<< "cfl = " << cfl << endl;
	cout<< "alpha = " << alpha << endl;
	cout<< "par_count = " << par_count << endl;
	
	if (par_count != 0)
	{
		for (k = 0; k < par_count; k++)
			cout<< "parameters("<< k <<") = " << parameters(k) << endl;
	}
	
	cout<< endl;
}
