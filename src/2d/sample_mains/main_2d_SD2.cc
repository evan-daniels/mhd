#include "centpack_2d_SD2.h"

using namespace std;
using namespace CENTPACK;

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	int id, p;
	double wtime;

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//  Record the starting time.
	
	if (id == 0)
		wtime = MPI_Wtime();
	
	centpack_2d_SD2(id, p);
	
	MPI_Barrier(MPI_COMM_WORLD);

	if (id == 0)
	{
		wtime = MPI_Wtime() - wtime;
	
		cout <<"\n";       
		cout <<"Wall clock elapsed seconds = "<< wtime << "\n";
		cout << "\n";
		cout << "CENTPACK_MPI:\n";
	}
	
	MPI_Finalize();
	
	return 0;
}
