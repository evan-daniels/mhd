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

	// void evolution_2d_SD2(doublearray3d& un, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha, const doublearray1d& parameters, const int& id, const int& p);
	// WORKING CT BENCHMARK v
	// void evolution_2d_SD2(doublearray3d& un, doublearray3d& B1_xf, doublearray3d& B2_yf, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha, const doublearray1d& parameters, const int& id, const int& p);
	void evolution_2d_SD2(doublearray3d& un, doublearray3d& B1_xf, doublearray3d& B2_yf, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha, const doublearray1d& parameters, const double& dt_fluid, const int& id, const int& p);

	void reconstruction_2d_SD2(doublearray3d& un, doublearray3d& u_N, doublearray3d& u_S, doublearray3d& u_E, doublearray3d& u_W, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface,const double& alpha);
	
	// LIMITERS -- minmod functions
	
	double sign(const double& x);

	double min(const double& x, const double& y);

	double minmod(const double& x, const double& y);
	
	double minmod3(const double& x, const double& y, const double& z);

	// void C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& jcurl, doublearray3d& C);
	// WORKING CT BENCHMARK v
	// void C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& jcurl, const doublearray3d& B1_xf, const doublearray3d& B2_yf, doublearray3d& C, doublearray3d& Ez_corner);
	void C_flux_2d_SD2(const doublearray3d& u_N, const doublearray3d& u_S, const doublearray3d& u_E, const doublearray3d& u_W, const doublearray1d& lambda, const doublearray1d& mu, const doublearray1d& parameters, const doublearray3d& jcurl, doublearray3d& C);

	void Hx_flux_2d_SD2(const doublearray1d& u_w, const doublearray1d& u_e, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& Hx);
	
	void Hy_flux_2d_SD2(const doublearray1d& u_s, const doublearray1d& u_n, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& Hy);
	
	//end of evolution routines
	
	void end_of_step_2d_SD2(const doublearray3d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& sum_t, const double& dt_cpu, const doublearray1d& parameters, const int& id, const int& p);

	void run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl, const int& id, const int& p);
	
////////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS -- User defined, see examples
////////////////////////////////////////////////////////////////////////////////

	void boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const int& id, const int& p);

	void boundary_conditions_ct(doublearray3d& f, const int& id, const int& p);

	void reconstruction_em(const doublearray3d& un, const doublearray3d& Jc, const doublearray3d& B1_xf, const doublearray3d& B2_yf, doublearray3d& B_N, doublearray3d& B_S, doublearray3d& B_E, doublearray3d& B_W, doublearray3d& Jc_N, doublearray3d& Jc_S, doublearray3d& Jc_E, doublearray3d& Jc_W, const doublearray1d& dx_cell, const doublearray1d& dx_interface, const doublearray1d& dy_cell, const doublearray1d& dy_interface, const double& alpha);

	void electric_field(const doublearray2d& rho_N, const doublearray2d& rho_S, const doublearray2d& rho_E, const doublearray2d& rho_W, const doublearray3d& v_N, const doublearray3d& v_S, const doublearray3d& v_E, const doublearray3d& v_W, const doublearray3d& B_N, const doublearray3d& B_S, const doublearray3d& B_E, const doublearray3d& B_W, const doublearray3d& J_N, const doublearray3d& J_S, const doublearray3d& J_E, const doublearray3d& J_W, doublearray3d& E_N, doublearray3d& E_S, doublearray3d& E_E, doublearray3d& E_W, const double& di, const double& eta);

	void C_flux_em(const doublearray2d& rho_N, const doublearray2d& rho_S, const doublearray2d& rho_E, const doublearray2d& rho_W, const doublearray3d& B_N, const doublearray3d& B_S, const doublearray3d& B_E, const doublearray3d& B_W, const doublearray3d& E_N, const doublearray3d& E_S, const doublearray3d& E_E, const doublearray3d& E_W, const doublearray1d& lambda_em, const doublearray1d& mu_em, const doublearray1d& dx_interface, const doublearray1d& dy_interface, const double& di, doublearray2d& Cx_em, doublearray2d& Cy_em, doublearray2d& Cz_em);

	void spectral_radii_em(const doublearray1d& u, const double& dx, const double& dy, const double& di, double& rx, double& ry);

	void flux_x(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& f);

	void flux_y(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& j_here, doublearray1d& g);
	
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
}

#endif