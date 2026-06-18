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
	
	if (id == 0)
	{
		tag = 3;
		MPI::COMM_WORLD.Send(&u(2,0), 2*L, MPI::DOUBLE, p - 1, tag);
	}
	
	if (id == p-1)
	{
		tag = 4;
		MPI::COMM_WORLD.Send(&u(J,0), 2*L, MPI::DOUBLE, 0, tag);
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
	
	if (id == p-1)
	{
		tag = 3;
		MPI::COMM_WORLD.Recv(&u(J+2,0), 2*L, MPI::DOUBLE, 0, tag, status);
	}
	
	
	if (id == 0)
	{
		tag = 4;
		MPI::COMM_WORLD.Recv(&u(0,0), 2*L, MPI::DOUBLE, p-1, tag, status);
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
	
	if (odd)
	{
		// Send messages
	
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
		
		if (id == 0)
		{
			tag = 3;
			MPI::COMM_WORLD.Send(&u(2,0), 2*L, MPI::DOUBLE, p - 1, tag);
		}
		
		if (id == p-1)
		{
			tag = 4;
			MPI::COMM_WORLD.Send(&u(J,0), L, MPI::DOUBLE, 0, tag);
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
		
		if (id == p-1)
		{
			tag = 3;
			MPI::COMM_WORLD.Recv(&u(J+2,0), 2*L, MPI::DOUBLE, 0, tag, status);
		}
		
		if (id == 0)
		{
			tag = 4;
			MPI::COMM_WORLD.Recv(&u(0,0), L, MPI::DOUBLE, p-1, tag, status);
		}
	}
	
	else
	{
		// Send messages
	
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
		
		if (id == 0)
		{
			tag = 3;
			MPI::COMM_WORLD.Send(&u(3,0), L, MPI::DOUBLE, p - 1, tag);
		}
		
		if (id == p-1)
		{
			tag = 4;
			MPI::COMM_WORLD.Send(&u(J,0), 2*L, MPI::DOUBLE, 0, tag);
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
		
		if (id == p-1)
		{
			tag = 3;
			MPI::COMM_WORLD.Recv(&u(J+3,0), L, MPI::DOUBLE, 0, tag, status);
		}
		
		if (id == 0)
		{
			tag = 4;
			MPI::COMM_WORLD.Recv(&u(0,0), 2*L, MPI::DOUBLE, p-1, tag, status);
		}
	}
}
