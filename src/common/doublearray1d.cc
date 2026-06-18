////////////////////////////////////////////////////////////////////////////////
//
//  doublearray1d.cc -- implementation of class doublearray1d
//                      build from the array classes provided by C. Anderson
//                      (UCLA)
//
//  data structures used by CentPack (C) 2006 Jorge Balbas and
//  Eitan Tadmor
//
////////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<iomanip>
#include <cstdio>
#include"doublearray1d.h"

using namespace std;

doublearray1d::doublearray1d()
{
	dataPtr = 0;
	internalAlloc = 0;
	index1Size = 0;
	index1Begin = 0;
	index1End = 0;
}

doublearray1d::doublearray1d(long size)
{
	dataPtr = 0;
	internalAlloc = 0; 
	initialize(size);
}

doublearray1d::doublearray1d(const doublearray1d& d)
{
	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;
	dataPtr = new double[index1Size];
	internalAlloc = 1; 
	
	long i;
	for(i = 0; i < index1Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

doublearray1d::~doublearray1d()
{
	if(internalAlloc == 1)
		delete [] dataPtr;
}

void doublearray1d::initialize(long m)
{
	if(internalAlloc == 1) 
    {
		if(index1Size != m)
		{
			delete [] dataPtr;
			dataPtr = new double[m];
		}
    }
	
    else 
    {
		if(dataPtr == 0)
		{
			dataPtr = new double[m];
			internalAlloc = 1;
		}
    }
    
    index1Size = m;
    index1Begin = 0;
    index1End = index1Begin + (index1Size - 1);
}

void doublearray1d::initialize(const doublearray1d& d)
{
    if(internalAlloc == 1) 
    {
		if(index1Size != d.index1Size)
		{
			delete [] dataPtr;
			dataPtr = new double[d.index1Size];
		}
    }
    
	else 
    {
		if(dataPtr == 0)
		{
			dataPtr = new double[d.index1Size];
			internalAlloc = 1;
		}
    }
    
    index1Size = d.index1Size;
    index1Begin = d.index1Begin;
    index1End = d.index1End;
    
    long i;
    
    for(i = 0; i < index1Size; i++)
    {
	    dataPtr[i] = d.dataPtr[i];
    }
}

ostream& operator << (ostream& out_stream, const doublearray1d& d)
{
	double outvalue;

	cout.setf(ios::scientific);
	cout.setf(ios::floatfield);
	cout.precision(16);
	
	long i;
	for(i = d.index1Begin; i <= d.index1End; i++)
	{
		outvalue = d.dataPtr[i];
		
		if(outvalue < 0 )
			out_stream << setprecision(16)<<outvalue <<" ";
		else
			out_stream << " " << setprecision(16)<<outvalue <<" ";
	}
	
	return out_stream;
}
