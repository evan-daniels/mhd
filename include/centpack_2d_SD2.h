#ifndef CENTPACK_2d_SD2_H
#define CENTPACK_2d_SD2_H

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
	
	int centpack_2d_SD2(int id, int p);
	
	void disclaimer();
	
	void read_in_2d(double& x_init, double& x_final, double& y_init, double& y_final, long& J, long& K, long& L, double& t_final, double& dt_out, double& cfl, double& alpha, doublearray1d& parameters);
	
	void time_step_2d(const doublearray3d& un, const doublearray1d& dx, const doublearray1d& dy, const double& cfl, double& dt, double& t, double& t_out, double& dt_out, doublearray1d& lambda, doublearray1d& mu, const doublearray1d& parameters);
	
	////////////////////////////////////////////////////////////////////////////
	// EVOLUTION
	////////////////////////////////////////////////////////////////////////////

	void evolution_2d_SD2(doublearray3d& un, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha, const doublearray1d& parameters, const doublearray3d& Efield, const int& id, const int& p);
	
	void reconstruction_2d_SD2(doublearray3d& un, doublearray3d& u_N, doublearray3d& u_S, doublearray3d& u_E, doublearray3d& u_W, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface,const double& alpha);
	
	// LIMITERS	-- minmod functions
	
	double sign(const double& x);

	double min(const double& x, const double& y);

	double minmod(const double& x, const double& y);
	
	double minmod3(const double& x, const double& y, const double& z);

	void C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& Efield, doublearray3d& C);

	void Hx_flux_2d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, const doublearray1d& parameters, const doublearray1d& e_here, doublearray1d& Hx);

	void Hy_flux_2d_SD2(const doublearray1d& u_s, const doublearray1d& u_n, const doublearray1d& parameters, const doublearray1d& e_here, doublearray1d& Hy);
	
	//end of evolution routines
	
	void end_of_step_2d_SD2(const doublearray3d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& sum_t, const double& dt_cpu, const doublearray1d& parameters, const int& id, const int& p);

	void run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl, const int& id, const int& p);
	
////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const int& id, const int& p);

	void flux_x(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& e_here, doublearray1d& f);

	void flux_y(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& e_here, doublearray1d& g);
	
	void mesh(const double& x_init, const double& x_final, const double& y_init, const double& y_final, doublearray1d& x, doublearray1d& x_cell, doublearray1d& dx_cell, doublearray1d& dx_interface, doublearray1d& y, doublearray1d& y_cell, doublearray1d& dy_cell, doublearray1d& dy_interface, const int& id, const int& p);
	
	void initial_conditions(doublearray3d& un, const doublearray1d& parameters, const doublearray1d& dx_cell, const doublearray1d& dy_cell, const doublearray1d& dx_interface, const doublearray1d& dy_interface, const doublearray1d& x_cell, const doublearray1d& y_cell, const doublearray1d& x, const doublearray1d& y);

	void spectral_radii(const doublearray1d& u, const doublearray1d& parameters, double& rx, double& ry);

	void writeout(const doublearray3d& un, const double& t, const doublearray1d& parameters, const long& n, const int& id, const int& p);
	
	void write_mesh(const doublearray1d& x, const doublearray1d& y, const int& id, const int& p);
	
	void resistivity_step(doublearray3d& un, const doublearray1d& dx_cell,
                      const doublearray1d& dy_cell, double dt,
                      const doublearray1d& parameters);

	void hall_step(doublearray3d& un, const doublearray1d& dx_cell,
							const doublearray1d& dy_cell, double dt,
							const doublearray1d& parameters);

	void electric_field(const doublearray3d& un,
                    const doublearray1d& dx_cell,
                    const doublearray1d& dy_cell,
                    const doublearray1d& parameters,
                    doublearray3d& E);
}

#endif
