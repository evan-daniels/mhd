#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;
using std::string;

void CENTPACK::writeout(const doublearray2d& un, const double& t, const doublearray1d& parameters, const long& n, const int& id, const int& p)
{
	
	long j, J;
	char u_file[30];
	char t_file[30];
	
	J = un.getIndex1Size() - 4;
	
	doublearray1d u(J);
	
	for (j = 0; j < J; j++)
	  u(j) = un(j+2,0);
	
	sprintf(u_file, "u_files/u_p%0.3d_%0.3ld", id, n);
	
	ofstream OutFile;
	OutFile.open(u_file, ios::out);
	OutFile<<u;
	OutFile.close();
	
	if (id == 0)
	{
	  sprintf(t_file, "t_files/t_%0.3ld", n);
	
	  OutFile.open(t_file, ios::out);
	  OutFile<<t;
	  OutFile.close();
	}
    
    cout.setf(ios::scientific, ios::floatfield);
}

void CENTPACK::writeout(const doublearray2d& un, const double& t, const doublearray1d& parameters, const long& n, const bool& odd, const int& id, const int& p)
{
	long j, J;
	char u_file[30];
	char t_file[30];
	char odd_file[30];

	
	J = un.getIndex1Size() - 4;
	
	doublearray1d u(J);
	
	for (j = 0; j < J; j++)
	  u(j) = un(j+2,0);
	
	sprintf(u_file, "u_files/u_p%0.3d_%0.3ld", id, n);
	
	ofstream OutFile;
	OutFile.open(u_file, ios::out);
	OutFile<<u;
	OutFile.close();
	
	if (id == 0)
	{
	  sprintf(t_file, "t_files/t_%0.3ld", n);
  	  sprintf(odd_file, "odd_files/odd_%0.3ld", n);
	  
	  OutFile.open(t_file, ios::out);
	  OutFile<< t;
	  OutFile.close();
	  
	  OutFile.open(odd_file, ios::out);
	  OutFile<< odd;
	  OutFile.close();
	}
    
    cout.setf(ios::scientific, ios::floatfield);
}
