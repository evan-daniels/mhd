#ifndef CENTPACK_1d_FD2_H
#define CENTPACK_1d_FD2_H

#include "arrays.h"
#include "mpi.h"
#include<algorithm>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>

namespace CENTPACK
{

////////////////////////////////////////////////////////////////////////////////
// MAIN CENTPACK FUNCTION 
////////////////////////////////////////////////////////////////////////////////
	
	int centpack_1d_FD2(int, int);

	void disclaimer();
	
	void read_in_1d(double& x_init, double& x_final, long& J, long& L, double& t_final, double& dt_out, double& cfl, double& alpha, doublearray1d& parameters);
	
	void time_step_1d(const doublearray2d& un, const doublearray1d& dx_cell, const double& cfl, double& dt, double& t, double& t_out, double& dt_out, doublearray1d& lambda, const doublearray1d& parameters);
	
	////////////////////////////////////////////////////////////////////////////
	// EVOLUTION
	////////////////////////////////////////////////////////////////////////////
	
	void evolution_1d_FD2(doublearray2d& un, const doublearray2d& uj_half, const doublearray1d& lambda, const double& alpha, const doublearray1d& parameters, const bool& odd, const int& id, const int& p);
	
	void reconstruction_1d_FD2(doublearray2d& un, doublearray2d& uj_half, const double& alpha, const bool& odd);
	
	void predictor_1d_FD2(const doublearray2d& un, doublearray2d& un_half, const double& alpha, const doublearray1d& lambda, const doublearray1d& parameters);
	
	void corrector_1d_FD2(doublearray2d& un, const doublearray2d& uj_half, const doublearray2d& un_half, const doublearray1d& lambda, const doublearray1d& parameters, const bool& odd);
	
	void end_of_step_1d_FD2(const doublearray2d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& sum_t, const double& dt_cpu, const doublearray1d& parameters, const bool& odd, const int& id, const int& p);
	
	double sign(const double& x);
	
	double min(const double& x, const double& y);
	
	double minmod(const double& x, const double& y);
	
	double minmod3(const double& x, const double& y, const double& z);
	
	void run_info_1d(double& dt, double& sum_t, long& J, double& cfl, const int& id, const int& p);
	
////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray2d& u, const bool& odd, const int& id, const int& p);

	void flux_x(const doublearray1d& v, const doublearray1d& parameters, doublearray1d& y);

	void mesh(const double& x_left, const double& x_right, doublearray1d& xe, doublearray1d& xo, doublearray1d& xe_cell, doublearray1d& xo_cell, doublearray1d& dxe_cell, doublearray1d& dxo_cell, doublearray1d& dxe_interface, doublearray1d& dxo_interface, const int& id, const int& p);

	void initial_conditions(doublearray2d& un, const doublearray1d& parameters, const doublearray1d& x);

	void spectral_radius(const doublearray1d& u, const doublearray1d& parameters, double& rx);

	void writeout(const doublearray2d& un, const double& t, const doublearray1d& parameters, const long& n, const bool& odd, const int& id, const int& p);
	
	void write_mesh(const doublearray1d& xe_cell, const doublearray1d& xo_cell, const int& id, const int& p);
}

#endif
