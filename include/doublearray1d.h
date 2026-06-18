////////////////////////////////////////////////////////////////////////////////
//
//  doublearray1d.h -- headers of class doublearray1d
//                     build from the array classes provided by C. Anderson
//                     (UCLA)
//
//  data structures used by CentPack (C) 2006 Jorge Balbas and
//  Eitan Tadmor
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __doublearray1d__
#define __doublearray1d__

#include<iostream>

using namespace std;

class doublearray1d
{
	friend ostream& operator << (ostream& outstream, const doublearray1d& d);	
	
	private:
	
	double* dataPtr;
	long    index1Begin;
	long    index1End;
	long    index1Size;
	int     internalAlloc;
	
	public:
	
	doublearray1d();
	doublearray1d(long size);
	doublearray1d(const doublearray1d& d);
	~doublearray1d();
	void initialize(long m);
	void initialize(const doublearray1d& d);
	
	inline double& operator()(long i1)
	{
		return *(dataPtr +  (i1 - index1Begin));
	}
	
	inline const double& operator()(long i1) const
	{
		return *(dataPtr +  (i1 - index1Begin));
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

	inline void resize(long newSize)
	{
		long i;
		double*  newDataPtr = new double[newSize];
		double*  tmpDataPtr;
	
		if(newSize > index1Size) 
		{
			for(i = 0; i < index1Size; i++)
				newDataPtr[i] = dataPtr[i];
		}
		
		else
		{
			for(i = 0; i < newSize; i++)
				newDataPtr[i] = dataPtr[i];
		}
	
		index1Size = newSize;
		tmpDataPtr = dataPtr;
		dataPtr    = newDataPtr;
	
		if(internalAlloc == 1) delete [] tmpDataPtr;
		internalAlloc = 1;
	
		index1End = index1Begin + (index1Size - 1);
	}
	
	inline void operator=(const doublearray1d& d)
	{
		initialize(d.index1Size);
		
		long i;
		for(i = 0; i < d.index1Size; i++)
		{
			dataPtr[i] = d.dataPtr[i];
		}
	}
	
	inline void setToValue(double val)
	{
		long i;
		for(i = 0; i < index1Size; i++)
		{
			dataPtr[i] =  val;
		}
	}
	
	
	inline void addValue(double val)
	{
		long i;
		
		for(i = 0; i < index1Size; i++)
		{
			dataPtr[i] += val;
		}
	}
};
#endif
