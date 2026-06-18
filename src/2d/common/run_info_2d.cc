// #include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl, const int& id, const int& p)
{
	ofstream OutFile;
	
	OutFile.open("run_info.txt");
	OutFile<<"dt = "<< dt <<endl;
	OutFile<<"total time = "<< sum_t <<endl;
	OutFile<<"grid size = "<< p*J <<" x "<< K << endl;
	OutFile<<"CFL = "<< cfl <<endl;
	OutFile<<"Number of processors: "<< p <<endl;
	OutFile<<"Processor "<< p - 1 <<" finished"<<endl;
	OutFile.close();
}
