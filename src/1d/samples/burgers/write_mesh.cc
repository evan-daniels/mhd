#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::write_mesh(const doublearray1d& x_cell, const int& id, const int& p)
{	
	long j, J;
	char x_file[30];
	ofstream OutFile;
	
	
	J = x_cell.getIndex1Size() - 4;
	
	doublearray1d x(J);
	
	for (j = 0; j < J; j++)
		x(j) = x_cell(j+2);
	
	sprintf(x_file, "mesh_files/x_p%0.3d", id);
	
	OutFile.open(x_file, ios::out);
	OutFile<<x;
	OutFile.close();
}

void CENTPACK::write_mesh(const doublearray1d& xe_cell, const doublearray1d& xo_cell, const int& id, const int& p)
{	
	long j, J;
	char xe_file[20], xo_file[20];
	ofstream OutFile;
	
	J = xe_cell.getIndex1Size() - 4;
	
	doublearray1d xe(J), xo(J);
	
	for (j = 0; j < J; j++)
	{
		xe(j) = xe_cell(j+2);
		xo(j) = xo_cell(j+1);
	}
	
	sprintf(xe_file, "mesh_files/xe_p%0.3d", id);	
	sprintf(xo_file, "mesh_files/xo_p%0.3d", id);	
	
	OutFile.open(xe_file, ios::out);
	OutFile<< xe;
	OutFile.close();
	
	OutFile.open(xo_file, ios::out);
	OutFile<< xo;
	OutFile.close();
}
