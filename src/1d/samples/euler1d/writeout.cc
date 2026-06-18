#include "centpack_1d_FD2.h"
#include "centpack_1d_SD2.h"

using namespace std;
using std::string;

void CENTPACK::writeout(const doublearray2d& un, const double& t, const doublearray1d& parameters, const long& n, const int& id, const int& p)
{

	cout.setf(ios::scientific, ios::floatfield);

	long j, J;
	char rho_file[30];
	char u1_file[30];
	char p_file[30];
	char t_file[30];

	double gamma = parameters(0);

	J = un.getIndex1Size() - 4;

	doublearray1d rho(J);
	doublearray1d u1(J);
	doublearray1d press(J);

	for (j=0; j<J; j++)
	{
	  rho(j) = un(j+2,0);
	  u1(j) = un(j+2,1)/un(j+2,0);
	  press(j) = (gamma - 1.0)*(un(j+2,2) - .5*(pow(un(j+2,1),2.0)/un(j+2,0)));
	}

	sprintf(rho_file, "rho_files/rho_p%0.3d_%0.3ld", id, n);
	sprintf(u1_file, "u1_files/u1_p%0.3d_%0.3ld", id, n);
	sprintf(p_file, "p_files/p_p%0.3d_%0.3ld", id, n);

	ofstream OutFile;
	OutFile.open(rho_file, ios::out);
	OutFile<< rho;
	OutFile.close();

	OutFile.open(u1_file, ios::out);
	OutFile<< u1;
	OutFile.close();

	OutFile.open(p_file, ios::out);
	OutFile<< press;
	OutFile.close();

	if (id == 0)
	{
	  sprintf(t_file, "t_files/t_%0.3ld", n);

	  OutFile.open(t_file, ios::out);
	  OutFile<<t;
	  OutFile.close();
	}
}

void CENTPACK::writeout(const doublearray2d& un, const double& t, const doublearray1d& parameters, const long& n, const bool& odd, const int& id, const int& p)
{

	cout.setf(ios::scientific, ios::floatfield);

	long j, J;
	char rho_file[30];
	char u1_file[30];
	char p_file[30];
	char t_file[30];
	char odd_file[30];

	double gamma = parameters(0);

	J = un.getIndex1Size() - 4;

	doublearray1d rho(J);
	doublearray1d u1(J);
	doublearray1d press(J);

	for (j=0; j<J; j++)
	{
		rho(j) = un(j+2,0);
		u1(j) = un(j+2,1)/un(j+2,0);
		press(j) = (gamma - 1.0)*(un(j+2,2) - .5*(pow(un(j+2,1),2.0)/un(j+2,0)));
	}

	sprintf(rho_file, "rho_files/rho_p%0.3d_%0.3ld", id, n);
	sprintf(u1_file, "u1_files/u1_p%0.3d_%0.3ld", id, n);
	sprintf(p_file, "p_files/p_p%0.3d_%0.3ld", id, n);

	ofstream OutFile;
	OutFile.open(rho_file, ios::out);
	OutFile<< rho;
	OutFile.close();

	OutFile.open(u1_file, ios::out);
	OutFile<< u1;
	OutFile.close();

	OutFile.open(p_file, ios::out);
	OutFile<< press;
	OutFile.close();

	if (id == 0)
	{
		sprintf(t_file, "t_files/t_%0.3ld", n);
		sprintf(odd_file, "odd_files/odd_%0.3ld", n);

		OutFile.open(t_file, ios::out);
		OutFile<<t;
		OutFile.close();

		OutFile.open(odd_file, ios::out);
		OutFile<<odd;
		OutFile.close();
	}
}
