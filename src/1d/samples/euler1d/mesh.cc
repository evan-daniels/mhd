#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::mesh(const double& x_left, const double& x_right, doublearray1d& x, doublearray1d& x_cell, doublearray1d& dx_cell, doublearray1d& dx_interface, const int& id, const int& p)
{
	long j, l, J, L;
	double xid0, xidf;
	double dx;
	
	J = x_cell.getIndex1Size() - 4;
	
	dx = (x_right - x_left)/(double)(p*J);
  
	xid0 = x_left + (double)((J*dx)*id);
	xidf = xid0 + J*dx;
  
	x(0) = xid0 - 2.0*dx;
	
	// cell interface points
	for (j = 1; j < J+5; j++)
		x(j) = x(j-1) + dx;
	
	// cell width
	for (j = 0; j < J+4; j++)
		dx_cell(j) = x(j+1) - x(j);
	
	// cell centers
	for (j = 0; j < J+4; j++)
		x_cell(j) = 0.5*(x(j) + x(j+1));
	
	// distance from consecutive cell centers
	for (j = 0; j < J+3; j++)
		dx_interface(j) = x_cell(j+1) - x_cell(j);
}

void CENTPACK::mesh(const double& x_left, const double& x_right, doublearray1d& xe, doublearray1d& xo, doublearray1d& xe_cell, doublearray1d& xo_cell, doublearray1d& dxe_cell, doublearray1d& dxo_cell, doublearray1d& dxe_interface, doublearray1d& dxo_interface, const int& id, const int& p)
{
	long j, l, J, L;
	double xid0, xidf;
	double dx;
	
	J = xe_cell.getIndex1Size() - 4;
	
	dx = (x_right - x_left)/(double)(p*J);
  
	xid0 = x_left + (double)((J*dx)*id);
	xidf = xid0 + J*dx;
	  
	xo(0) = xid0 - 1.5*dx;
	xe(0) = xid0 - 2.0*dx;
	
	// cell interface points
	for (j = 1; j < J+5; j++)
	{
		xo(j) = xo(j-1) + dx;
		xe(j) = xe(j-1) + dx;
	}
	
	// cell width
	for (j = 0; j < J+4; j++)
	{
		dxo_cell(j) = xo(j+1) - xo(j);
		dxe_cell(j) = xe(j+1) - xe(j);
	}
	
	// cell centers
	for (j = 0; j < J+4; j++)
	{
		xo_cell(j) = 0.5*(xo(j) + xo(j+1));
		xe_cell(j) = 0.5*(xe(j) + xe(j+1));
	}
	
	// distance from consecutive cell centers
	for (j = 0; j < J+3; j++)
	{
		dxo_interface(j) = xo_cell(j+1) - xo_cell(j);
		dxe_interface(j) = xe_cell(j+1) - xe_cell(j);
	}
}
