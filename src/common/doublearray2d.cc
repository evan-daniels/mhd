////////////////////////////////////////////////////////////////////////////////
//
//  doublearray2d.cc -- implementation of class doublearray2d
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
#include<cstdio>
#include"doublearray2d.h"

using namespace std;

doublearray2d::doublearray2d()
{
	dataPtr = 0;
	internalAlloc = 0;
	index1Size = 0;
	index1Begin = 0;
	index1End = 0;
	index2Size = 0;
	index2Begin = 0;
	index2End = 0;
}

doublearray2d::doublearray2d(long size1, long size2)
{
	dataPtr = 0;
	internalAlloc = 0; 
	initialize(size1, size2);
}

doublearray2d::doublearray2d(const doublearray2d& d)
{
	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;
	index2Size = d.index2Size;
	index2Begin = d.index2Begin;
	index2End = d.index2End;
	
	dataPtr = new double[index1Size*index2Size];
	internalAlloc = 1; 
	
	long i;
	for(i = 0; i < index1Size*index2Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

doublearray2d::~doublearray2d()
{
	if(internalAlloc == 1)
		delete [] dataPtr;
}

void doublearray2d::initialize(long size1, long size2)
{
	if(internalAlloc == 1) 
	{
		if((index1Size != size1) || (index2Size != size2))
		{
			delete [] dataPtr;
			dataPtr = new double[size1*size2];
		}
	}
	
	else
	{
		if(dataPtr == 0)
		{
			dataPtr = new double[size1*size2];
			internalAlloc  = 1;
		}
	}
	
	index1Size = size1;
	index1Begin = 0;
	index1End = index1Begin + (index1Size - 1);
	index2Size = size2;
	index2Begin = 0;
	index2End = index2Begin + (index2Size - 1);
	
}

void doublearray2d::initialize(const doublearray2d& d)
{
    if(internalAlloc == 1) 
    {
		if((index1Size != d.index1Size) ||(index2Size != d.index2Size))
		{
			delete [] dataPtr;
			dataPtr = new double[d.index1Size*d.index2Size];
		}
    }
    
	else 
    {
		if(dataPtr == 0)
		{
			dataPtr = new double[d.index1Size*d.index2Size];
			internalAlloc = 1;
		}
    }
    
    index1Size = d.index1Size;
    index1Begin = d.index1Begin;
    index1End = d.index1End;
    index2Size = d.index2Size;
    index2Begin = d.index2Begin;
    index2End = d.index2End;
    
    long i;
    
    for(i = 0; i < index1Size*index2Size; i++)
    {
	    dataPtr[i] = d.dataPtr[i];
    }
}

ostream& operator<< (ostream& out_stream, const doublearray2d& d)
{
	double outvalue;

	cout.setf(ios::scientific, ios::floatfield);
	cout.precision(16);
	
	long i,j,k;
	for(i = d.index1Begin; i <= d.index1End; i++)
	{
		for(j = d.index2Begin; j <= d.index2End; j++)
		{
			k=(j - d.index2Begin) + (i - d.index1Begin)*d.index2Size;
			outvalue=d.dataPtr[k] ;
		
			if(outvalue < 0 )
				out_stream << setprecision(16)<<outvalue <<" ";
			else
				out_stream <<" "<< setprecision(16)<<outvalue <<" ";
		}
	
		out_stream<<endl;
	}
	
	return out_stream;
}
