////////////////////////////////////////////////////////////////////////////////
//
// writeout.cc -- an auxiliary function of
//
// CentPack (C) 2005 Jorge Balbas and Eitan Tadmor
//
// Function called by centpack_main_2d_FD2()
//
// This function outputs files containing the conserved solution variables (or
// corresponding primitive variables) for the ideal MHD model.
//
// Input (passed by reference):
//
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t over the discretized solution domain
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// (3) t -- time of solution being output
//
// Output (outfiles):
//
// (1) Nine text files holding the solution variables at time t (n keeps track
//     of how many output calls have been executed over the time of the
//     simulation
//
//     rho_n -- density
//
//     u1_n -- velocity in x-direction
//
//     u2_n -- velocity in y-direction
//
//     u3_n -- velocity in z-direction
//
//     b1_n -- magnetic field z-direction
//
//     b2_n -- magnetic field x-direction
//
//     b3_n -- magnetic field z-direction
//
//     p_n -- pressure
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

	char u_file[30];
	char t_file[30];

	J = un.getIndex1Size() - 4;
	K = un.getIndex2Size() - 4;

	doublearray2d u(J,K);

	for (j = 0; j < J; j++)
	{
		for (k = 0; k < K; k++)
			u(j,k) = un(j+2,k+2,0);
	}

	sprintf(u_file, "u_files/u_p%0.3d_%0.3ld", id, n);

	ofstream OutFile;
	OutFile.open(u_file, ios::out);
	OutFile<<u;
	OutFile.close();

	if (id == 0)
	{
	  sprintf(t_file, "t_files/t_%0.3ld", n);

	  OutFile.open(t_file, ios::out);
	  OutFile<<t;
	  OutFile.close();
	}
}

void CENTPACK::writeout(const doublearray3d& un, const double& t, const doublearray1d& parameters, const long& n, const bool& odd, const int& id, const int& p)
{

	cout.setf(ios::scientific, ios::floatfield);

	long j, k, J, K;

	char u_file[30];
	char t_file[30];
	char odd_file[30];

	J = un.getIndex1Size() - 4;
	K = un.getIndex2Size() - 4;

	doublearray2d u(J,K);

	for (j = 0; j < J; j++)
	{
	  for (k = 0; k < K; k++)
			u(j,k) = un(j+2,k+2,0);
	}

	sprintf(u_file, "u_files/u_p%0.3d_%0.3ld", id, n);

	ofstream OutFile;
	OutFile.open(u_file, ios::out);
	OutFile<< u;
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
