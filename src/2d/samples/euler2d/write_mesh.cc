////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
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
// write_mesh.cc
//
// Function called by centpack_2d_SD2, problem specific.
//
// This function writes the one-dimensional arrays containing the x- and 
// y-coordinates of the cell centers of the computational domain
// 
// Input (passed by reference):
//
// (1)  x_cell -- a doublearray1d type variable holding the x coordinates of
//           	 the computational mesh on each subdomain
//
// (2)  y_cell -- a doublearray1d type variable holding the y coordinates of the 
//				 computational mesh, same y-mesh for all subdomains
//
// Output: outfiles x_p and y
//
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::write_mesh(const doublearray1d& x_cell, const doublearray1d& y_cell, const int& id, const int& p)
{	
	long j, k, J, K;
	char x_file[30];
	char y_file[20];
	ofstream OutFile;
	
	J = x_cell.getIndex1Size() - 4;
	K = y_cell.getIndex1Size() - 4;
	
	doublearray1d x(J), y(K);
	
	for (j = 0; j < J; j++)
		x(j) = x_cell(j+2);
	
	sprintf(x_file, "mesh_files/x_p%0.3d", id);
	
	OutFile.open(x_file, ios::out);
	OutFile<<x;
	OutFile.close();
	
	if (id == 0)
	{
		for (k = 0; k < K; k++)
			y(k) = y_cell(k+2);

		sprintf(y_file, "mesh_files/y");
		OutFile.open(y_file, ios::out);
		OutFile<<y;
		OutFile.close();
	}
}

void CENTPACK::write_mesh(const doublearray1d& xe_cell, const doublearray1d& xo_cell, const doublearray1d& ye_cell, const doublearray1d& yo_cell, const int& id, const int& p)
{	
	long j, k, J, K;
	char xe_file[30];
	char xo_file[30];
	char ye_file[20];
	char yo_file[20];
	ofstream OutFile;
	
	J = xe_cell.getIndex1Size() - 4;
	K = ye_cell.getIndex1Size() - 4;
	
	doublearray1d xe(J), xo(J), ye(K), yo(K);
	
	for (j = 0; j < J; j++)
	{
		xe(j) = xe_cell(j+2);
		xo(j) = xo_cell(j+1);
	}
	
	sprintf(xe_file, "mesh_files/xe_p%0.3d", id);	
	sprintf(xo_file, "mesh_files/xo_p%0.3d", id);
	
	OutFile.open(xe_file, ios::out);
	OutFile<<xe;
	OutFile.close();
	
	OutFile.open(xo_file, ios::out);
	OutFile<<xo;
	OutFile.close();
	
	if (id == 0)
	{
		for (k = 0; k < K; k++)
		{
			ye(k) = ye_cell(k+2);
			yo(k) = yo_cell(k+1);
		}
		
		sprintf(ye_file, "mesh_files/ye");
		OutFile.open(ye_file, ios::out);
		OutFile<<ye;
		OutFile.close();
		
		sprintf(yo_file, "mesh_files/yo");
		OutFile.open(yo_file, ios::out);
		OutFile<<yo;
		OutFile.close();
	}
}
