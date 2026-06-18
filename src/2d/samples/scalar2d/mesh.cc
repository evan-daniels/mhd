#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;
using namespace CENTPACK;

void CENTPACK::mesh(const double& x_init, const double& x_final, const double& y_init, const double& y_final, doublearray1d& x, doublearray1d& x_cell, doublearray1d& dx_cell, doublearray1d& dx_interface, doublearray1d& y, doublearray1d& y_cell, doublearray1d& dy_cell, doublearray1d& dy_interface, const int& id, const int& p)
{
	long j, k, l, J, K, L;
	double dx, dy;
	double xid0, xidf;
	
	J = x_cell.getIndex1Size() - 4;
	K = y_cell.getIndex1Size() - 4;
	
	dx = (x_final - x_init)/(double)(p*J);
	dy = (y_final - y_init)/K;
	
	xid0 = x_init + (double)((J*dx)*id);
	xidf = xid0 + J*dx;
  
	x(0) = xid0 - 2.0*dx;
	y(0) = y_init - 2.0*dy;
	
	// cell interface points
	for (j = 1; j < J + 5; j++)
		x(j) = x(j-1) + dx;
	
	for (k = 1; k < K + 5; k++)
		y(k) = y(k-1) + dy;
	
	// cell width
	for (j = 0; j < J + 4; j++)
		dx_cell(j) = x(j+1) - x(j);
	
	for (k = 0; k < K + 4; k++)
		dy_cell(k) = y(k+1) - y(k);
	
	// cell centers
	for (j = 0; j < J + 4; j++)
		x_cell(j) = 0.5*(x(j) + x(j+1));
	
	for (k = 0; k < K + 4; k++)
		y_cell(k) = 0.5*(y(k) + y(k+1));
	
	// distance from consecutive cell centers
	for (j = 0; j < J + 3; j++)
		dx_interface(j) = x_cell(j+1) - x_cell(j);

	for (k = 0; k < K + 3; k++)
		dy_interface(k) = y_cell(k+1) - y_cell(k);
}

void CENTPACK::mesh(const double& x_left, const double& x_right, const double& y_bottom, const double& y_top, doublearray1d& xe, doublearray1d& xo, doublearray1d& xe_cell, doublearray1d& xo_cell, doublearray1d& dxe_cell, doublearray1d& dxo_cell, doublearray1d& dxe_interface, doublearray1d& dxo_interface, doublearray1d& ye, doublearray1d& yo, doublearray1d& ye_cell, doublearray1d& yo_cell, doublearray1d& dye_cell, doublearray1d& dyo_cell, doublearray1d& dye_interface, doublearray1d& dyo_interface, const int& id, const int& p)
{
	long j, k, l, J, K, L;
	double dx, dy;
	double xid0, xidf;
	
	J = xe_cell.getIndex1Size() - 4;
	K = ye_cell.getIndex1Size() - 4;
	
	dx = (x_right - x_left)/(double)(p*J);
	dy = (y_top - y_bottom)/K;
	
	xid0 = x_left + (double)((J*dx)*id);
	xidf = xid0 + J*dx;
  
	xo(0) = xid0 - 1.5*dx;
	xe(0) = xid0 - 2.0*dx;
	
	yo(0) = y_bottom - 1.5*dy;
	ye(0) = y_bottom - 2.0*dy;

	
	// cell interface points
	for (j = 1; j < J + 5; j++)
	{
		xo(j) = xo(j-1) + dx;
		xe(j) = xe(j-1) + dx;
	}
	
	for (k = 1; k < K + 5; k++)
	{
		yo(k) = yo(k-1) + dy;
		ye(k) = ye(k-1) + dy;
	}
	
	// cell width
	for (j = 0; j < J + 4; j++)
	{
		dxo_cell(j) = xo(j+1) - xo(j);
		dxe_cell(j) = xe(j+1) - xe(j);

	}
	
	for (k = 0; k < K + 4; k++)
	{
		dyo_cell(k) = yo(k+1) - yo(k);
		dye_cell(k) = ye(k+1) - ye(k);
	}
	
	// cell centers
	for (j = 0; j < J + 4; j++)
	{
		xo_cell(j) = 0.5*(xo(j) + xo(j+1));
		xe_cell(j) = 0.5*(xe(j) + xe(j+1));

	}
	
	for (k = 0; k < K + 4; k++)
	{
		yo_cell(k) = 0.5*(yo(k) + yo(k+1));
		ye_cell(k) = 0.5*(ye(k) + ye(k+1));
	}
	
	// distance from consecutive cell centers
	for (j = 0; j < J + 3; j++)
	{
		dxo_interface(j) = xo_cell(j+1) - xo_cell(j);
		dxe_interface(j) = xe_cell(j+1) - xe_cell(j);
	}
	
	for (k = 0; k < K + 3; k++)
	{
		dyo_interface(k) = yo_cell(k+1) - yo_cell(k);
		dye_interface(k) = ye_cell(k+1) - ye_cell(k);		
	}
}
