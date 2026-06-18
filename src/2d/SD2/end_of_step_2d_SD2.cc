////////////////////////////////////////////////////////////////////////////
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
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.
//
// end_of_step_2d_SD2d.cc -- a core function of CentPack
//
// Function called by solver (main) at the end of each time step; it checks
// whether output needs to be produced and prints monitoring information to the
// screen
//
// Input (read from file input):
// 
// (1) un  --  a doublearray3d type variable holding the cell averages of u
//             at time t=t^n over the discretized solution domain
//
// (2) dt  --  a double type variable holding the size of the time step used in 
//             last time iteration
//
// (3) t  --  a double type variable holding the simulation time completed after 
//            last time iteration
//
// (4) dt_out  --  a double type variable holding the valued of the desired time 
//                 interval between succesive outputs (user specified)
//
// (5) t_out  --  a double type variable holding the simulation time completed 
//                since last output was produced
//
// (6) n  -- a long type variable indicating the number of outputs produced 
//           since beginning of simulation
//
// (7) gamma  --  a double type variable holding the ratio of specific heats
//
// (8) dt_cpu  --  a double type variable holding the CPU time of last time 
//                 iteration
//
// (9) sum_t  --  a double type variable holding the total CPU time employed so //                far
//
// Output:
//
// (1) end_of_step prints a line to the screen (it can be redirected to a file),
//     with the following pieces or information regarding the simulation:
//
//     (a) run -- present simulation time
//  
//     (b) dt -- time step used in the last evolution step
//  
//     (c) dt_cpu -- CPU time employed by last evolution step
//
//     (d) t -- total CPU time
//
// (2) The solutions variables are output via the function writeout, called at
//     the time interval indicated by the user in the input file
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::end_of_step_2d_SD2(const doublearray3d& un, const double& dt, const double& t, double& dt_out, double& t_out, long& n, const double& sum_t, const double& dt_cpu, const doublearray1d& parameters, const int& id, const int& p)
{
	
	if (t_out == dt_out)
	{
		writeout(un, t, parameters, n, id, p);
		n++;
		
		t_out = 0.0;
		
		cout<<"output written at t = "<< t <<endl;
	}
	
	cout<<"run="<<t<<",dt="<<dt<<",cpu_t="<<dt_cpu<<",t="<<sum_t<<endl;
}
