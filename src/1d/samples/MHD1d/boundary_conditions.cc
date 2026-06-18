#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::boundary_conditions(doublearray2d& u, const int& id, const int& p)
{
	long J = u.getIndex1Size() - 4;
	long L = u.getIndex2Size();
	long j, l;

	MPI::Status status;
	int tag;
	int index;
	
	doublearray1d u_left(L), u_right(L);
	
	// Send messages
	
	if (id > 0)
	{
		tag = 1;
		MPI::COMM_WORLD.Send(&u(2,0), 2*L, MPI::DOUBLE, id-1, tag);
	}
	
	if (id < p-1)
	{
		tag = 2;
		MPI::COMM_WORLD.Send(&u(J,0), 2*L, MPI::DOUBLE, id+1, tag);
	}
	
	// Receive messages
	
	if (id < p-1)
	{
		tag = 1;
		MPI::COMM_WORLD.Recv(&u(J+2,0), 2*L, MPI::DOUBLE, id+1, tag, status);
	}
		
	if (id > 0)
	{
		tag = 2;
		MPI::COMM_WORLD.Recv(&u(0,0), 2*L, MPI::DOUBLE, id-1, tag, status);
	}
	
	if (id == 0)
	{
		u_left(0) = 1.0;
		u_left(1) = 0.0;
		u_left(2) = 0.0;
		u_left(3) = 0.0;
		u_left(4) = 1.0;
		u_left(5) = 0.0;
		u_left(6) = 1.0 + (25.0/32.0);
		
		for (l = 0; l < L; l++)
		{
			u(0,l) = u_left(l);
			u(1,l) = u_left(l);
		}
	}
	
	if (id == p - 1)
	{
		u_right(0) = 0.125; 
		u_right(1) = 0.0;
		u_right(2) = 0.0;
		u_right(3) = 0.0;
		u_right(4) = -1.0;
		u_right(5) = 0.0;
		u_right(6) = 0.1 + (25.0/32.0);
		
		for (l = 0; l < L; l++)
		{
			u(J+2,l) = u_right(l);
			u(J+3,l) = u_right(l);
		}
	}
}

void CENTPACK::boundary_conditions(doublearray2d& u, const bool& odd, const int& id, const int& p)
{
	long J = u.getIndex1Size() - 4;
	long L = u.getIndex2Size();
	long j, l;

	MPI::Status status;
	int tag;
	int index;
	
	doublearray1d u_left(L), u_right(L);
	
	// Send messages
	if (odd)
	{
		if (id > 0)
		{
			tag = 1;
			MPI::COMM_WORLD.Send(&u(2,0), 2*L, MPI::DOUBLE, id-1, tag);
		}
		
		if (id < p-1)
		{
			tag = 2;
			MPI::COMM_WORLD.Send(&u(J,0), L, MPI::DOUBLE, id+1, tag);
		}
		
		// Receive messages
		
		if (id < p-1)
		{
			tag = 1;
			MPI::COMM_WORLD.Recv(&u(J+2,0), 2*L, MPI::DOUBLE, id+1, tag, status);
		}
			
		if (id > 0)
		{
			tag = 2;
			MPI::COMM_WORLD.Recv(&u(0,0), L, MPI::DOUBLE, id-1, tag, status);
		}
	  
		if (id == 0)
		{
			u_left(0) = 1.0;
			u_left(1) = 0.0;
			u_left(2) = 0.0;
			u_left(3) = 0.0;
			u_left(4) = 1.0;
			u_left(5) = 0.0;
			u_left(6) = 1.0 + (25.0/32.0);
			
			for (l = 0; l < L; l++)
				u(0,l) = u_left(l);
		}
		
		if (id == p - 1)
		{
			u_right(0) = 0.125; 
			u_right(1) = 0.0;
			u_right(2) = 0.0;
			u_right(3) = 0.0;
			u_right(4) = -1.0;
			u_right(5) = 0.0;
			u_right(6) = 0.1 + (25.0/32.0);
			
			for (l = 0; l < L; l++)
			{
				u(J+2,l) = u_right(l);
				u(J+3,l) = u_right(l);
			}
		}
	}
	
	else
	{
		if (id > 0)
		{
			tag = 1;
			MPI::COMM_WORLD.Send(&u(3,0), L, MPI::DOUBLE, id-1, tag);
		}
		
		if (id < p-1)
		{
			tag = 2;
			MPI::COMM_WORLD.Send(&u(J,0), 2*L, MPI::DOUBLE, id+1, tag);
		}
		
		// Receive messages
		
		if (id < p-1)
		{
			tag = 1;
			MPI::COMM_WORLD.Recv(&u(J+3,0), L, MPI::DOUBLE, id+1, tag, status);
		}
			
		if (id > 0)
		{
			tag = 2;
			MPI::COMM_WORLD.Recv(&u(0,0), 2*L, MPI::DOUBLE, id-1, tag, status);
		}
	  
		if (id == 0)
		{
			u_left(0) = 1.0;
			u_left(1) = 0.0;
			u_left(2) = 0.0;
			u_left(3) = 0.0;
			u_left(4) = 1.0;
			u_left(5) = 0.0;
			u_left(6) = 1.0 + (25.0/32.0);
			
			for (l = 0; l < L; l++)
			{
				u(0,l) = u_left(l);
				u(1,l) = u_left(l);
			}
		}
		
		if (id == p - 1)
		{
			u_right(0) = 0.125; 
			u_right(1) = 0.0;
			u_right(2) = 0.0;
			u_right(3) = 0.0;
			u_right(4) = -1.0;
			u_right(5) = 0.0;
			u_right(6) = 0.1 + (25.0/32.0);
			
			for (l = 0; l < L; l++)
				u(J+3,l) = u_right(l);
		}
	}
}
