////////////////////////////////////////////////////////////////////////////////
//
//  doublearray3d.cc -- implementation of class doublearray3d
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
#include"doublearray3d.h"

using namespace std;

doublearray3d::doublearray3d()
{
	dataPtr = 0;
	internalAlloc = 0;
	index1Size = 0;
	index1Begin = 0;
	index1End = 0;
	index2Size = 0;
	index2Begin = 0;
	index2End = 0;
	index3Size = 0;
	index3Begin = 0;
	index3End = 0;
}

doublearray3d::doublearray3d(long size1, long size2, long size3)
{
	dataPtr = 0;
	internalAlloc = 0;
	initialize(size1, size2, size3);
}

doublearray3d::doublearray3d(const doublearray3d& d)
{
	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;
	index2Size = d.index2Size;
	index2Begin = d.index2Begin;
	index2End = d.index2End;
	index3Size = d.index3Size;
	index3Begin = d.index3Begin;
	index3End = d.index3End;
	
	dataPtr = new double[index1Size*index2Size*index3Size];
	internalAlloc = 1;
	
	long i;
	for(i = 0; i < index1Size*index2Size*index3Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

doublearray3d::~doublearray3d()
{
	if(internalAlloc == 1)
		delete [] dataPtr;
}

void doublearray3d::initialize(long size1, long size2, long size3)
{
	if(internalAlloc == 1) 
	{
		if((index1Size != size1) || (index2Size != size2) || (index3Size != size3))
		{
			delete [] dataPtr;
			dataPtr = new double[size1*size2*size3];
		}
	}
	
	else
	{
		if(dataPtr == 0)
		{
			dataPtr = new double[size1*size2*size3];
			internalAlloc  = 1;
		}
	}
	
	index1Size = size1;
	index1Begin = 0;
	index1End = index1Begin + (index1Size - 1);
	index2Size = size2;
	index2Begin = 0;
	index2End = index2Begin + (index2Size - 1);
	index3Size = size3;
	index3Begin = 0;
	index3End = index3Begin + (index3Size - 1);
	
}

void doublearray3d::initialize(const doublearray3d& d)
{
    if(internalAlloc == 1) 
    {
		if((index1Size != d.index1Size) || (index2Size != d.index2Size) || (index3Size != d.index3Size))
		{
			delete [] dataPtr;
			dataPtr = new double[d.index1Size*d.index2Size*d.index3Size];
		}
    }
    
    else
    {
	    if(dataPtr == 0)
	    {
		    dataPtr = new double[d.index1Size*d.index2Size*d.index3Size];
		    internalAlloc = 1;
	    }
    }
    
    index1Size = d.index1Size;
    index1Begin = d.index1Begin;
    index1End = d.index1End;
    index2Size = d.index2Size;
    index2Begin = d.index2Begin;
    index2End = d.index2End;
    index3Size = d.index3Size;
    index3Begin = d.index3Begin;
    index3End = d.index3End;
    
    long i;
    
    for(i = 0; i < index1Size*index2Size*index3Size; i++)
    {
	    dataPtr[i] = d.dataPtr[i];
    }
}

ostream& operator<< (ostream& out_stream, const doublearray3d& d)
{
	double outvalue;

	cout.setf(ios::scientific);
	cout.setf(ios::floatfield);
	cout.precision(24);
	
	long i,j,k,l;
	
	for (l = 0; l < d.index1Size*d.index2Size*d.index3Size; l++)
    {
		outvalue = d.dataPtr[l];
		
		if(outvalue < 0.0)
			out_stream<<" "<<setprecision(24)<<outvalue<<" ";
		else
			out_stream<<" "<<setprecision(24)<<outvalue<<" ";
		out_stream<<std::endl;
    }
	
	return out_stream;
}
