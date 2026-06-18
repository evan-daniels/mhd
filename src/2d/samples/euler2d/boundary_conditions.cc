////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2006 Jorge Balbas and Eitan Tadmor
//
// boundary_conditions.cc
//
// Function called by evolution, problem specific.
//
// This function implements periodic boundary conditions at the four boundaries
// of a rectangular domain by filling the first and last two rows/columns of
// the solution array u according to the periodic conditions u(x,z)=u(x+p,z+p).
// 
// Input (passed by reference):
//
// (1)  u -- a doublearray3d type variable holding the conserved
//           quantities density, momentum (3 components), magnetic field
//           (3 components), and total energy over the entire soltion domain.
//
// (2) parameters -- a doublearray1d type variable holding additional parameters 
//					for the specific example
//
// (3) id -- process id number
//
// (4) p -- number of processors for the simulation
//
// Output (returned by reference):
//
// (1) u -- a doublearray3d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy over the entire soltion domain.
//
////////////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const int& id, const int& p)
{
	
	long J = u.getIndex1Size() - 4;
	long K = u.getIndex2Size() - 4;
	long L = u.getIndex3Size();
	long j, k, l;
	
	MPI::Status status;
	int tag;
	int index;
	
	double p_one = 1.5;
	double p_two = 0.3;
	double p_three = 0.029;
	double p_four = 0.3;
	
	double gamma = parameters(0);
	
	doublearray1d u_one(L), u_two(L), u_three(L), u_four(L);
	
	u_one(0) = 1.5;
	u_one(1) = 0.0;
	u_one(2) = 0.0;
	u_one(3) = p_one/(gamma - 1.0) + 0.5*(pow(u_one(1),2.0) + pow(u_one(2),2.0))/u_one(0);
	
	u_two(0) = 0.5323;
	u_two(1) = 1.206*u_two(0);
	u_two(2) = 0.0;
	u_two(3) = p_two/(gamma - 1.0) + 0.5*(pow(u_two(1),2.0) + pow(u_two(2),2.0))/u_two(0);
	
	u_three(0) = 0.138;
	u_three(1) = 1.206*u_three(0);
	u_three(2) = 1.206*u_three(0);
	u_three(3) = p_three/(gamma - 1.0) + 0.5*(pow(u_three(1),2.0) + pow(u_three(2),2.0))/u_three(0);
	
	u_four(0) = 0.5323;
	u_four(1) = 0.0;
	u_four(2) = 1.206*u_four(0);
	u_four(3) = p_four/(gamma - 1.0) + 0.5*(pow(u_four(1),2.0) + pow(u_four(2),2.0))/u_four(0);
	
	// y-boundary conditions
	
	for (l = 0; l < L; l++)
	{
		for (j = 2; j < J+2; j++)
		{
			for (k = 0; k < 2; k++)
				u(j,k,l) = u(j,2,l);
			for (k = K+2; k < K+4; k++)
				u(j,k,l) = u(j,K+1,l);
		}
	}
	
	// x-boundaries () require communication between processors
	// Send messages
	
	// for processes 1, 2, ..., p-1, the left two nonboundary columns of un in 
	// process id are the the right two boundary columns of un in process id-1

	int left  = (id > 0)   ? id-1 : MPI::PROC_NULL;
	int right = (id < p-1) ? id+1 : MPI::PROC_NULL;

	if (id > 0)
	{
		tag = 1;
		MPI::COMM_WORLD.Send(&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, tag);
	}

	// for processes 0, 1, ..., p-2, the right two nonboundary columns of un in
	// process id are the the left two boundary columns of un in process id-1

	
	if (id < p-1)
	{
		tag = 2;
		MPI::COMM_WORLD.Send(&u(J,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, tag);
	}
	
	// Receive messages
	
	if (id < p-1)
	{
		tag = 1;
		MPI::COMM_WORLD.Recv(&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, tag, status);
	}
		
	if (id > 0)
	{
		tag = 2;
		MPI::COMM_WORLD.Recv(&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, tag, status);
	}

	// boundary conditions at left boundary and left corners
	
	if (id == 0)
	{
		for (l = 0; l < L; l++)
		{
			for (k = 2; k < K+2; k++)
			{
				for (j = 0; j < 2; j++)
					u(j,k,l) = u(2,k,l);
			}
			
			u(0,K+2,l) = u_two(l);
			u(0,K+3,l) = u_two(l);
			u(1,K+2,l) = u_two(l);
			u(1,K+3,l) = u_two(l);
			
			u(0,0,l) = u_three(l);
			u(0,1,l) = u_three(l);
			u(1,0,l) = u_three(l);
			u(1,1,l) = u_three(l);
		}
	}
	
	// boundary conditions at right boundary and right corners

	if (id == p - 1)
	{
		
		for (l = 0; l < L; l++)
		{
			for (k = 2; k < K+2; k++)
			{
				for (j = J+2; j < J+4; j++)
					u(j,k,l) = u(J+1,k,l);
			}
			
			u(J+2,K+2,l) = u_one(l);
			u(J+3,K+2,l) = u_one(l);
			u(J+2,K+3,l) = u_one(l);
			u(J+3,K+3,l) = u_one(l);
			
			u(J+2,0,l) = u_four(l);
			u(J+3,0,l) = u_four(l);
			u(J+2,1,l) = u_four(l);
			u(J+3,1,l) = u_four(l);
		}
	}
}

void CENTPACK::boundary_conditions(doublearray3d& u, const doublearray1d& parameters, const bool& odd, const int& id, const int& p)
{

	long J = u.getIndex1Size() - 4;
	long K = u.getIndex2Size() - 4;
	long L = u.getIndex3Size();
	long j, k, l;
	
	MPI::Status status;
	int tag;
	int index;
	
	double p_one = 1.5;
	double p_two = 0.3;
	double p_three = 0.029;
	double p_four = 0.3;
	
	double gamma = parameters(0);
	
	doublearray1d u_one(L), u_two(L), u_three(L), u_four(L);
	
	u_one(0) = 1.5;
	u_one(1) = 0.0;
	u_one(2) = 0.0;
	u_one(3) = p_one/(gamma - 1.0) + 0.5*(pow(u_one(1),2.0) + pow(u_one(2),2.0))/u_one(0);
	
	u_two(0) = 0.5323;
	u_two(1) = 1.206*u_two(0);
	u_two(2) = 0.0;
	u_two(3) = p_two/(gamma - 1.0) + 0.5*(pow(u_two(1),2.0) + pow(u_two(2),2.0))/u_two(0);
	
	u_three(0) = 0.138;
	u_three(1) = 1.206*u_three(0);
	u_three(2) = 1.206*u_three(0);
	u_three(3) = p_three/(gamma - 1.0) + 0.5*(pow(u_three(1),2.0) + pow(u_three(2),2.0))/u_three(0);
	
	u_four(0) = 0.5323;
	u_four(1) = 0.0;
	u_four(2) = 1.206*u_four(0);
	u_four(3) = p_four/(gamma - 1.0) + 0.5*(pow(u_four(1),2.0) + pow(u_four(2),2.0))/u_four(0);
	
	
	if (odd)
	{
		// y-boundary conditions
		
		for (l = 0; l < L; l++)
		{
			for (j = 1; j < J+2; j++)
			{
				u(j,0,l) = u(j,1,l);
				
				for (k = K+2; k < K+4; k++)
					u(j,k,l) = u(j,K+1,l);
			}
		}
		
		int left  = (id > 0)   ? id-1 : MPI::PROC_NULL;
		int right = (id < p-1) ? id+1 : MPI::PROC_NULL;

		if (id > 0)
		{
			tag = 1;
			MPI::COMM_WORLD.Send(&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, tag);
		}
		
		if (id < p-1)
		{
			tag = 2;
			MPI::COMM_WORLD.Send(&u(J,0,0), L*(K+4), MPI::DOUBLE, id+1, tag);
		}
		
		// Receive messages
		
		if (id < p-1)
		{
			tag = 1;
			MPI::COMM_WORLD.Recv(&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, tag, status);
		}
			
		if (id > 0)
		{
			tag = 2;
			MPI::COMM_WORLD.Recv(&u(0,0,0), L*(K+4), MPI::DOUBLE, id-1, tag, status);
		}
		
		// boundary conditions at left boundary and left corners
		
		if (id == 0)
		{
			for (l = 0; l < L; l++)
			{
				for (k = 1; k < K+2; k++)
					u(0,k,l) = u(1,k,l);
				
				u(0,K+2,l) = u_two(l);
				u(0,K+3,l) = u_two(l);
				
				u(0,0,l) = u_three(l);
			}
		}
		
		// boundary conditions at right boundary and right corners
	
		if (id == p - 1)
		{
			
			for (l = 0; l < L; l++)
			{
				for (k = 1; k < K+2; k++)
				{
					for (j = J+2; j < J+4; j++)
						u(j,k,l) = u(J+1,k,l);
				}
				
				u(J+2,K+2,l) = u_one(l);
				u(J+3,K+2,l) = u_one(l);
				u(J+2,K+3,l) = u_one(l);
				u(J+3,K+3,l) = u_one(l);
				
				u(J+2,0,l) = u_four(l);
				u(J+3,0,l) = u_four(l);
			}
		}
	}
	
	else
	{
		// y-boundary conditions
		
		for (l = 0; l < L; l++)
		{
			for (j = 2; j < J+3; j++)
			{
				for (k = 0; k < 2; k++)
					u(j,k,l) = u(j,2,l);
				
				u(j,K+3,l) = u(j,K+2,l);
			}
		}

		int left  = (id > 0)   ? id-1 : MPI::PROC_NULL;
		int right = (id < p-1) ? id+1 : MPI::PROC_NULL;

		if (id > 0)
		{
			tag = 1;
			MPI::COMM_WORLD.Send(&u(3,0,0), L*(K+4), MPI::DOUBLE, id-1, tag);
		}
		
		if (id < p-1)
		{
			tag = 2;
			MPI::COMM_WORLD.Send(&u(J,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, tag);
		}
		
		// Receive messages
		
		if (id < p-1)
		{
			tag = 1;
			MPI::COMM_WORLD.Recv(&u(J+3,0,0), L*(K+4), MPI::DOUBLE, id+1, tag, status);
		}
			
		if (id > 0)
		{
			tag = 2;
			MPI::COMM_WORLD.Recv(&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, id-1, tag, status);
		}
		
		// boundary conditions at left boundary and left corners
		
		if (id == 0)
		{
			for (l = 0; l < L; l++)
			{
				for (k = 2; k < K+3; k++)
				{
					for (j = 0; j < 2; j++)
						u(j,k,l) = u(2,k,l);
				}
				
				u(0,K+3,l) = u_two(l);
				u(1,K+3,l) = u_two(l);
				
				u(0,0,l) = u_three(l);
				u(0,1,l) = u_three(l);
				u(1,0,l) = u_three(l);
				u(1,1,l) = u_three(l);
			}
		}
		
		// boundary conditions at right boundary and right corners
	
		if (id == p - 1)
		{
			
			for (l = 0; l < L; l++)
			{
				for (k = 2; k < K+3; k++)
					u(J+3,k,l) = u(J+2,k,l);
				
				u(J+3,K+3,l) = u_one(l);
				
				u(J+3,0,l) = u_four(l);
				u(J+3,1,l) = u_four(l);
			}
		}
	}
}

/* void CENTPACK::boundary_conditions(doublearray3d& u, const bool& odd)
{
	long J = u.getIndex1Size() - 4;
	long K = u.getIndex2Size() - 4;
	long L = u.getIndex3Size();
	long j, k, l;
	
	long j1 = J/2 + 2;
	long k1 = K/2 + 2;
	
	double p_one = 1.5;
	double p_two = 0.3;
	double p_three = 0.029;
	double p_four = 0.3;
	
	double gamma = 1.4;
	
	doublearray1d u_one(L), u_two(L), u_three(L), u_four(L);
	
	u_one(0) = 1.5;
	u_one(1) = 0.0;
	u_one(2) = 0.0;
	u_one(3) = p_one/(gamma - 1.0) + 0.5*(pow(u_one(1),2.0) + pow(u_one(2),2.0))/u_one(0);
	
	u_two(0) = 0.5323;
	u_two(1) = 1.206*u_two(0);
	u_two(2) = 0.0;
	u_two(3) = p_two/(gamma - 1.0) + 0.5*(pow(u_two(1),2.0) + pow(u_two(2),2.0))/u_two(0);;
	
	u_three(0) = 0.138;
	u_three(1) = 1.206*u_three(0);
	u_three(2) = 1.206*u_three(0);
	u_three(3) = p_three/(gamma - 1.0) + 0.5*(pow(u_three(1),2.0) + pow(u_three(2),2.0))/u_three(0);
	
	u_four(0) = 0.5323;
	u_four(1) = 0.0;
	u_four(2) = 1.206*u_four(0);
	u_four(3) = p_four/(gamma - 1.0) + 0.5*(pow(u_four(1),2.0) + pow(u_four(2),2.0))/u_four(0);
	
	if (odd)
	{
		for (l = 0; l < L; l++)
		{
			for (j = 1; j < J+2; j++)
			{
				u(j,0,l) = u(j,1,l);
				
				for (k = K+2; k < K+4; k++)
					u(j,k,l) = u(j,K+1,l);
			}
			
			for (k = 1; k < K+2; k++)
			{
				u(0,k,l) = u(1,k,l);
				
				for (j = J+2; j < J+4; j++)
					u(j,k,l) = u(J+1,k,l);
			}
			
			u(J+2,K+2,l) = u_one(l);
			u(J+3,K+2,l) = u_one(l);
			u(J+2,K+3,l) = u_one(l);
			u(J+3,K+3,l) = u_one(l);
			
			u(0,K+2,l) = u_two(l);
			u(0,K+3,l) = u_two(l);
			
			u(0,0,l) = u_three(l);
			u(0,1,l) = u_three(l);
			
			u(J+2,0,l) = u_four(l);
			u(J+3,0,l) = u_four(l);
			u(J+2,1,l) = u_four(l);
			u(J+3,1,l) = u_four(l);
		}
	}
	
	else
	{
		for (l = 0; l < L; l++)
		{
			for (j = 2; j < J+3; j++)
			{
				for (k = 0; k < 2; k++)
					u(j,k,l) = u(j,2,l);
				
				u(j,K+3,l) = u(j,K+2,l);
			}
			
			for (k = 2; k < K+3; k++)
			{
				for (j = 0; j < 2; j++)
					u(j,k,l) = u(2,k,l);
				
				u(J+3,k,l) = u(J+2,k,l);
			}
			
			u(J+3,K+2,l) = u_one(l);
			u(J+3,K+3,l) = u_one(l);
			
			u(0,K+2,l) = u_two(l);
			u(0,K+3,l) = u_two(l);
			u(1,K+2,l) = u_two(l);
			u(1,K+3,l) = u_two(l);
			
			u(0,0,l) = u_three(l);
			u(1,0,l) = u_three(l);
			u(0,1,l) = u_three(l);
			u(1,1,l) = u_three(l);
			
			u(J+3,0,l) = u_four(l);
			u(J+3,1,l) = u_four(l);
		}
	}
} */
