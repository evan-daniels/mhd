#ifndef CENTPACK_1d_SD2_H
#define CENTPACK_1d_SD2_H

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

	int centpack_1d_SD2(int, int);

	void disclaimer();

	void read_in_1d(double& x_init, double& x_final, long& J, long& L, double& t_final, double& dt_out, double& cfl, double& alpha, doublearray1d& parameters);

	void time_step_1d(const doublearray2d& un, const doublearray1d& dx_cell, const double& cfl, double& dt, double& t, double& t_out, double& dt_out, doublearray1d& lambda, const doublearray1d& parameters);

	////////////////////////////////////////////////////////////////////////////
	// EVOLUTION
	////////////////////////////////////////////////////////////////////////////

	void evolution_1d_SD2(doublearray2d& un, const doublearray1d& lambda, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const double& alpha, const doublearray1d& parameters, const int& id, const int& p);

	// for SD formulation, reconstruction is embeded in evolution

	void reconstruction_1d_SD2(doublearray2d& un, doublearray2d& u_E, doublearray2d& u_W, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const double& alpha);

	void C_flux_1d_SD2(const doublearray2d& u_E, const doublearray2d& u_W, const doublearray1d& lambda, const doublearray1d& parameters, doublearray2d& C);

	void H_flux_1d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, const doublearray1d& parameters, doublearray1d& H);

	void end_of_step_1d_SD2(const doublearray2d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& sum_t, const double& dt_cpu, const doublearray1d& parameters, const int& id, const int& p);

	double sign(const double& x);

	double min(const double& x, const double& y);

	double minmod(const double& x, const double& y);

	double minmod3(const double& x, const double& y, const double& z);

	void run_info_1d(double& dt, double& sum_t, long& J, double& cfl, const int& id, const int& p);

////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray2d& u, const int& id, const int& p);

	void flux_x(const doublearray1d& v, const doublearray1d& parameters, doublearray1d& y);

	void mesh(const double& x_init, const double& x_final, doublearray1d& x, doublearray1d& x_cell, doublearray1d& dx_cell, doublearray1d& dx_interface, const int& id, const int& p);

	void initial_conditions(doublearray2d& un, const doublearray1d& parameters, const doublearray1d& x);

	void spectral_radius(const doublearray1d& u, const doublearray1d& parameters, double& rx);

	void writeout(const doublearray2d& un, const double& t, const doublearray1d& parameters, const long& n, const int& id, const int& p);

	void write_mesh(const doublearray1d& x, const int& id, const int& p);
}

#endif
