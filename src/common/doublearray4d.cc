////////////////////////////////////////////////////////////////////////////////
//
//  doublearray4d.cc -- body of class doublearray4d
//                     build partially from the array classes provided by C. Anderson(UCLA)
//                     
//
// data structures used by Central Hyperbolic Solver (C) Jorge Balbas and Eitan Tadmor, Nov. 2004
// Qian Xin reconstruct the 4d class from the 3D class for 3D MHD development in centPack
///////////////////////////////////////////////////////////////////////////////


#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdio>
#include"doublearray4d.h"

using namespace std;

doublearray4d::doublearray4d()
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

    index4Size = 0;
    index4Begin = 0;
    index4End = 0;
}

doublearray4d::doublearray4d(long size1, long size2, long size3,long size4)
{
	dataPtr = 0;
	internalAlloc = 0;
	initialize(size1, size2, size3,size4);
}

doublearray4d::doublearray4d(const doublearray4d& d)
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

	index4Size=d.index4Size;
    index4Begin=d.index4Begin;
    index4End=d.index4End;

	dataPtr = new double[bufLength()];
	internalAlloc = 1;
	
	long i;
	for(i = 0; i < bufLength(); i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

doublearray4d::~doublearray4d()
{
	if(internalAlloc == 1)
		delete [] dataPtr;
}

void doublearray4d::initialize(long size1, long size2, long size3,long size4)
{
	index1Size = size1;
	index1Begin = 0;
	index1End = index1Begin + (index1Size - 1);

	index2Size = size2;
	index2Begin = 0;
	index2End = index2Begin + (index2Size - 1);

	index3Size = size3;
	index3Begin = 0;
	index3End = index3Begin + (index3Size - 1);

    index4Size = size4;
    index4Begin = 0;
    index4End = index4Begin + (index4Size - 1);

	if(internalAlloc == 1) 
	{
		if((index1Size != size1) || (index2Size != size2) || (index3Size != size3) || (index4Size!=size4) )
		{
			delete [] dataPtr;
			dataPtr = new double[bufLength()];
		}
	}
	
	else
	{
		if(dataPtr == NULL)
		{
			dataPtr = new double[bufLength()];
			internalAlloc  = 1;
		}
	}	

}

void doublearray4d::initialize(const doublearray4d& d)
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

    index4Size = d.index4Size;
    index4Begin = d.index4Begin;
    index4End = d.index4End;

    if(internalAlloc == 1) 
    {
		if((index1Size != d.index1Size) || (index2Size != d.index2Size) || (index3Size != d.index3Size)||(index4Size != d.index4Size))
		{
			delete [] dataPtr;
			dataPtr = new double[bufLength()];
		}
    }
    
    else
    {
	    if(dataPtr == 0)
	    {
		    dataPtr = new double[bufLength()];
		    internalAlloc = 1;
	    }
    }
    
    long i;
    
    for(i = 0; i < bufLength(); i++)
    {
	    dataPtr[i] = d.dataPtr[i];
    }
}


ostream& operator<< (ostream& out_stream, const doublearray4d& d)
{
	double outvalue;

	cout.setf(ios::scientific);
	cout.setf(ios::floatfield);
	cout.precision(24);
	
	long i;
    
	for(i = d.index1Begin; i <= d.bufLength(); i++)
	{
        outvalue=d.dataPtr[i];				

		if(outvalue < 0.0)
			out_stream<<" "<<setprecision(24)<<outvalue<<" ";
		else
			out_stream<<" "<<setprecision(24)<<outvalue<<" ";		

            out_stream<<std::endl;
	}	

	return out_stream;
}
