#include "centpack_1d_SD2.h"

using namespace std;
using namespace CENTPACK;

int main(int argc, char *argv[])
{
	int id;
	int p;
	double wtime;
	
	MPI::Init (argc, argv);
	
	id = MPI::COMM_WORLD.Get_rank ( );
	
	p = MPI::COMM_WORLD.Get_size ( );
	
	//  Record the starting time.
	
	if (id == 0) 
	{
		wtime = MPI::Wtime ( );
	}
	
	centpack_1d_SD2(id, p);
	
	if (id == 0)
	{
		wtime = MPI::Wtime ( ) - wtime;
	
		cout <<"\n";       
		cout <<"Wall clock elapsed seconds = "<< wtime << "\n";
	}
	
	MPI::Finalize ( );
	
	cout << "\n";
	cout << "CENTPACK_MPI:\n";
	
	return 0;
}
