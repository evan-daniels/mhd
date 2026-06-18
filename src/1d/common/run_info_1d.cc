#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;

void CENTPACK::run_info_1d(double& dt, double& sum_t, long& J, double& cfl, const int& id, const int& p)
{
	ofstream OutFile;
	
	OutFile.open("run_info.txt");
	OutFile<<"dt = "<< dt <<endl;
	OutFile<<"total time = "<< sum_t <<endl;
	OutFile<<"grid size = "<< J << endl;
	OutFile<<"CFL = "<< cfl <<endl;
	OutFile<<"Number of processors: "<< p <<endl;
	OutFile<<"Processor "<< p - 1 <<" finished"<<endl;
	OutFile.close();
}
