////////////////////////////////////////////////////////////////////////////////
//
//  doublearray3d.h -- headers of class doublearray3d
//                     build from the array classes provided by C. Anderson
//                     (UCLA)
//
//  data structures used by CentPack (C) 2006 Jorge Balbas and
//  Eitan Tadmor
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __doublearray3d__
#define __doublearray3d__

#include<iostream>

using namespace std;

class doublearray3d
{
	friend ostream& operator << (ostream& outstream, const doublearray3d& d);	
	
	private:
	
	double* dataPtr;
	long    index1Begin;
	long    index1End;
	long    index1Size;
	long    index2Begin;
	long    index2End;
	long    index2Size;
	long    index3Begin;
	long    index3End;
	long    index3Size;
	int     internalAlloc;
	
	public:
	
	doublearray3d();
	doublearray3d(long size1, long size2, long size3);
	doublearray3d(const doublearray3d& d);
	~doublearray3d();
	void initialize(long size1, long size2, long size3);
	void initialize(const doublearray3d& d);
	
	inline double& operator()(long i1, long i2, long i3)
	{
		return *(dataPtr + (i3 - index3Begin) 
		+ index3Size*((i2 - index2Begin) + (i1 - index1Begin)*index2Size));
	}
	
	inline const double& operator()(long i1, long i2, long i3) const
	{
		return *(dataPtr + (i3 - index3Begin) 
		+ index3Size*((i2 - index2Begin) + (i1 - index1Begin)*index2Size));
	}
	
	inline double* getDataPointer()
	{
		return dataPtr;
	}
	
	inline void setIndex1Begin(long i)
	{
		index1Begin = i;
		index1End   = index1Begin + (index1Size - 1);
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
		index2End   = index2Begin + (index2Size - 1);
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
	
	inline void setIndex3Begin(long i)
	{
		index3Begin = i;
		index3End   = index3Begin + (index3Size - 1);
	}
	
	inline long getIndex3Begin() const
	{
		return index3Begin;
	}
	
	inline long getIndex3End() const
	{
		return index3End;
	}
	
	inline long getIndex3Size() const
	{
		return index3Size;
	}
	
	inline void operator=(const doublearray3d& d)
	{
		if(index1Size*index2Size*index3Size == 0)
			initialize(d.index1Size, d.index2Size, d.index3Size);
		
		long i;
		
		for(i = 0; i < d.index1Size*d.index2Size*d.index3Size; i++)
			dataPtr[i] = d.dataPtr[i];
	}
	
	inline void setToValue(double val)
	{
		long i;
		for(i = 0; i < index1Size*index2Size*index3Size; i++)
		{
			dataPtr[i] =  val;
		}
	}
	
	inline void addValue(double val)
	{
		long i;
		
		for(i = 0; i < index1Size*index2Size*index3Size; i++)
		{
			dataPtr[i] += val;
		}
	}
};

#endif
