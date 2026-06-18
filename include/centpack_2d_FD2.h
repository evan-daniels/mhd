#ifndef CENTPACK_2d_FD2_H
#define CENTPACK_2d_FD2_H

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

namespace CENTPACK
{

////////////////////////////////////////////////////////////////////////////////
// MAIN CENTPACK FUNCTION 
////////////////////////////////////////////////////////////////////////////////
	
	int centpack_2d_FD2(int id, int p);

	void disclaimer();
	
	void read_in_2d(double& x_left, double& x_right, double& y_bottom, double& y_top, long& J, long& K, long& L, double& t_final, double& dt_out, double& cfl, double& alpha, doublearray1d& parameters);
	
	void time_step_2d(const doublearray3d& un, const doublearray1d& dx_cell, const doublearray1d& dy_cell, const double& cfl, double& dtp, double& t, double& t_out, double& dt_out, doublearray1d& lambda, doublearray1d& mu, const doublearray1d& parameters);
	
	////////////////////////////////////////////////////////////////////////////
	// RECONSTRUCTION
	////////////////////////////////////////////////////////////////////////////
	
	void reconstruction_2d_FD2(const doublearray3d& un, doublearray3d& unhalf, const double& alpha, const bool& odd);

	// LIMITERS	-- minmod functions
	
	double sign(const double& x);

	double min(const double& x, const double& y);

	double minmod(const double& x, const double& y);

	double minmod3(const double& x, const double& y, const double& z);

	// end of reconstruction
	
	////////////////////////////////////////////////////////////////////////////
	// EVOLUTION
	////////////////////////////////////////////////////////////////////////////

	void evolution_2d_FD2(doublearray3d& un, const doublearray3d& unhalf, const double& alpha, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const bool& odd, const int& id, const int& p);

	void predictor_2d_FD2(const doublearray3d& un, doublearray3d& u_star, const double& alpha, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters);
	
	void corrector_2d_FD2(doublearray3d& un, const doublearray3d& unhalf, const doublearray3d& u_star, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const bool& odd);
	
	//end of evolution routines
	
	void end_of_step_2d_FD2(const doublearray3d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& sum_t, const double& dt_cpu, const doublearray1d& parameters, const bool& odd, const int& id, const int& p);

	void run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl, const int& id, const int& p);
	
////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const bool& odd, const int& id, const int& p);

	void flux_x(const doublearray1d& u, const doublearray1d& parameters, doublearray1d& f);

	void flux_y(const doublearray1d& u, const doublearray1d& parameters, doublearray1d& g);
	
	void mesh(const double& x_left, const double& x_right, const double& y_bottom, const double& y_top, doublearray1d& xe, doublearray1d& xo, doublearray1d& xe_cell, doublearray1d& xo_cell, doublearray1d& dxe_cell, doublearray1d& dxo_cell, doublearray1d& dxe_interface, doublearray1d& dxo_interface, doublearray1d& ye, doublearray1d& yo, doublearray1d& ye_cell, doublearray1d& yo_cell, doublearray1d& dye_cell, doublearray1d& dyo_cell, doublearray1d& dye_interface, doublearray1d& dyo_interface, const int& id, const int& p);
	
	void initial_conditions(doublearray3d& un, const doublearray1d& parameters, const doublearray1d& dx_cell, const doublearray1d& dy_cell, const doublearray1d& dx_interface, const doublearray1d& dy_interface, const doublearray1d& x_cell, const doublearray1d& y_cell, const doublearray1d& x, const doublearray1d& y);

	void spectral_radii(const doublearray1d& u, const doublearray1d& parameters, double& rx, double& ry);

	void writeout(const doublearray3d& un, const double& t, const doublearray1d& parameters, const long& n, const bool& odd, const int& id, const int& p);
	
	void write_mesh(const doublearray1d& xe_cell, const doublearray1d& xo_cell, const doublearray1d& ye_cell, const doublearray1d& yo_cell, const int& id, const int& p);
}

#endif
