////////////////////////////////////////////////////////////////////////////////
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
// boundary_conditions.cc
//
// Function called by evolution, problem specific.
//
// This function implements periodic boundary conditions at the four boundaries
// of a rectangular domain by filling the first and last two rows/columns of
// the solution array u according to the periodic conditions u(x,z) = 
// u(x+P,z+P).
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
	
	// y-boundaries (top and bottom) are the same for all processors and require 
	// no comunication
	
	for (l = 0; l < L; l++)
	{
		for (j = 2; j < J+2; j++)
		{
			u(j,0,l) = u(j,K,l);
			u(j,1,l) = u(j,K+1,l);
			u(j,K+2,l) = u(j,2,l);
			u(j,K+3,l) = u(j,3,l);
		}
	}
	
	// x-boundaries () require communication between processors
	// Send messages
	
	// for processes 1, 2, ..., p-1, the left two nonboundary columns of un in 
	// process id are the the right two boundary columns of un in process id-1
	
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
	
	// periodic boundary conditions: for process 0, the left two nonboundary 
	// columns of un are the the right boundary coulumns of process p-1
	
	if (id == 0)
	{
		tag = 3;
		MPI::COMM_WORLD.Send(&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, tag);
	}
	
	// periodic boundary conditions: for process p-1 the right two nonboundary 
	// columns of un are the the two left boundary coulumns of process 0
	
	if (id == p-1)
	{	
		tag = 4;
		MPI::COMM_WORLD.Send(&u(J,0,0), 2*L*(K+4), MPI::DOUBLE, 0, tag);
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
	
	if (id == p-1)
	{
		tag = 3;
		MPI::COMM_WORLD.Recv(&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, 0, tag, status);
	}
	
	
	if (id == 0)
	{
		tag = 4;
		MPI::COMM_WORLD.Recv(&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, tag, status);
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
	
	if (odd)
	{
		// y-boundary conditions
		
		for (l = 0; l < L; l++)
		{
			for (j = 1; j < J+2; j++)
			{
				u(j,0,l) = u(j,K,l);
				u(j,K+2,l) = u(j,2,l);
				u(j,K+3,l) = u(j,3,l);
			}
		}
		
		// x-boundaries () require communication between processors
		// Send messages
		
		// for processes 1, 2, ..., p-1, the left two nonboundary columns of un in 
		// process id are the the right two boundary columns of un in process id-1
		
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
		
		if (id == 0)
		{
			tag = 3;
			MPI::COMM_WORLD.Send(&u(2,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, tag);
		}
		
		if (id == p-1)
		{
			tag = 4;
			MPI::COMM_WORLD.Send(&u(J,0,0), L*(K+4), MPI::DOUBLE, 0, tag);
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
		
		if (id == p-1)
		{
			tag = 3;
			MPI::COMM_WORLD.Recv(&u(J+2,0,0), 2*L*(K+4), MPI::DOUBLE, 0, tag, status);
		}
			
		if (id == 0)
		{
			tag = 4;
			MPI::COMM_WORLD.Recv(&u(0,0,0), L*(K+4), MPI::DOUBLE, p-1, tag, status);
		}
	}
	
	else
	{
		// y-boundary conditions
		
		for (l = 0; l < L; l++)
		{
			for (j = 2; j < J+3; j++)
			{
				u(j,0,l) = u(j,K,l);
				u(j,1,l) = u(j,K+1,l);
				u(j,K+3,l) = u(j,3,l);
			}
		}
		
		// x-boundaries () require communication between processors
		// Send messages
		
		// for processes 1, 2, ..., p-1, the left two nonboundary columns of un in 
		// process id are the the right two boundary columns of un in process id-1
		
		if (id > 0)
		{
			tag = 1;
			MPI::COMM_WORLD.Send(&u(3,0,0), L*(K+4), MPI::DOUBLE, id-1, tag);
		}
		
		// for processes 0, 1, ..., p-2, the right two nonboundary columns of un in
		// process id are the the left two boundary columns of un in process id-1
		
		if (id < p-1)
		{
			tag = 2;
			MPI::COMM_WORLD.Send(&u(J,0,0), 2*L*(K+4), MPI::DOUBLE, id+1, tag);
		}
		
		if (id == 0)
		{
			tag = 3;
			MPI::COMM_WORLD.Send(&u(3,0,0), L*(K+4), MPI::DOUBLE, p-1, tag);
		}
		
		if (id == p-1)
		{
			tag = 4;
			MPI::COMM_WORLD.Send(&u(J,0,0), 2*L*(K+4), MPI::DOUBLE, 0, tag);
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
		
		if (id == p-1)
		{
			tag = 3;
			MPI::COMM_WORLD.Recv(&u(J+3,0,0), L*(K+4), MPI::DOUBLE, 0, tag, status);
		}
			
		if (id == 0)
		{
			tag = 4;
			MPI::COMM_WORLD.Recv(&u(0,0,0), 2*L*(K+4), MPI::DOUBLE, p-1, tag, status);
		}
	}
}

/* void CENTPACK::boundary_conditions(doublearray3d& u, const bool& odd)
{
  long J=u.getIndex1Size() - 4;
  long K=u.getIndex2Size() - 4;
  long L=u.getIndex3Size();
  long j, k, l;
  
  if (odd)
  {
	  for (l=0; l<L; l++)
	  {
		  for (j=1; j<J+2; j++)
		  {
			  u(j,0,l)=u(j,K,l);
			  u(j,K+2,l)=u(j,2,l);
			  u(j,K+3,l)=u(j,3,l);
		  }
		  
		  for (k=0; k<K+4; k++)
		  {
			  u(0,k,l)=u(J,k,l);
			  u(J+2,k,l)=u(2,k,l);
			  u(J+3,k,l)=u(3,k,l);
		  }
	  }
  }
  
  else
  {
	  for (l=0; l<L; l++)
	  {
		  for (j=2; j<J+3; j++)
		  {
			  u(j,0,l)=u(j,K,l);
			  u(j,1,l)=u(j,K+1,l);
			  u(j,K+3,l)=u(j,3,l);
		  }
		  
		  for (k=0; k<K+4; k++)
		  {
			  u(0,k,l)=u(J,k,l);
			  u(1,k,l)=u(J+1,k,l);
			  u(J+3,k,l)=u(3,k,l);
		  }
	  }
  }
} */
