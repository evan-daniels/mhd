////////////////////////////////////////////////////////////////////////////////
//
//  doublearray2d.h -- headers of class doublearray2d
//                     build from the array classes provided by C. Anderson
//                     (UCLA)
//
// data structures used by CentraPack (C) 2006 Jorge Balbas and
// Eitan Tadmor
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __doublearray2d__
#define __doublearray2d__

#include<iostream>

using namespace std;

class doublearray2d
{
	friend ostream& operator << (ostream& outstream, const doublearray2d& d);	
	
	private:
	
	double* dataPtr;
	long index1Begin;
	long index1End;
	long index1Size;
	long index2Begin;
	long index2End;
	long index2Size;
	int internalAlloc;
	
	public:
	
	doublearray2d();
	doublearray2d(long size1, long size2);
	doublearray2d(const doublearray2d& d);
	~doublearray2d();
	void initialize(long size1, long size2);
	void initialize(const doublearray2d& d);
	
	inline double& operator()(long i1, long i2)
	{
		return *(dataPtr + (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
	}
	
	inline const double& operator()(long i1, long i2) const
	{
		return *(dataPtr + (i2 - index2Begin) + (i1 - index1Begin)*index2Size);
	}
	
	inline double* getDataPointer()
	{
		return dataPtr;
	}
	
	inline void setIndex1Begin(long i)
	{
		index1Begin = i;
		index1End = index1Begin + (index1Size - 1);
	}
	
	inline long getIndex1Begin() const
	{
		return index1Begin;
	}
	
	inline long getIndex1End() const
	{
		return index1End;
	}
	
	inline long getIndex1Size() const
	{
		return index1Size;
	}
	
	inline void setIndex2Begin(long i)
	{
		index2Begin = i;
		index2End = index2Begin + (index2Size - 1);
	}
	
	inline long getIndex2Begin() const
	{
		return index2Begin;
	}
	
	inline long getIndex2End() const
	{
		return index2End;
	}
	
	inline long getIndex2Size() const
	{
		return index2Size;
	}
	
	inline void operator=(const doublearray2d& d)
	{
		if(index1Size*index2Size == 0)
			initialize(d.index1Size,d.index2Size);
		
		long i;
		
		for(i = 0; i < d.index1Size*d.index2Size; i++)
			dataPtr[i] = d.dataPtr[i];
	}
	
	inline void setToValue(double val)
	{
		long i;
		for(i = 0; i < index1Size*index2Size; i++)
		{
			dataPtr[i] =  val;
		}
	}
	
	
	inline void addValue(double val)
	{
		long i;
		
		for(i = 0; i < index1Size*index2Size; i++)
		{
			dataPtr[i] += val;
		}
	}
};
#endif
