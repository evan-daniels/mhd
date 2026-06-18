////////////////////////////////////////////////////////////////////////////////
//
//  doublearray4d.h -- headers of class doublearray4d
//                     build partially from the array classes provided by C. 
//					   Anderson(UCLA)
//                     
//
// data structures used by CentPack (C) 2006 Jorge Balbas and Eitan Tadmor, 
// Qian Xin reconstruct the 4d class from the 3D class for 3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __doublearray4d__
#define __doublearray4d__

#include<iostream>
#include<cassert>

using namespace std;

class doublearray4d
{
	friend ostream& operator << (ostream& outstream, const doublearray4d& d);	
	
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
    long    index4Begin;
    long    index4End;
    long    index4Size;
	int     internalAlloc;
	
	public:
	
	doublearray4d();
	doublearray4d(long size1, long size2, long size3,long size4);
	doublearray4d(const doublearray4d& d); //cooy constructor
	~doublearray4d();

	void initialize(long size1, long size2, long size3,long size4);
	void initialize(const doublearray4d& d);

#if 0
    inline double& doublearray4d::operator()(long i1, long i2, long i3,long i4)
    {
    	return *(dataPtr + 
                  (i4 - index4Begin)   
                    index4Size(
                     (i3 - index3Begin) 
                        + index3Size*((i2 - index2Begin) + index2Size*(i1 - index1Begin))));
    }

    inline void setIndex1Begin(long i)
    {
       index1Begin = i;
       index1End   = index1Begin + (index1Size - 1);
    }

    inline void setIndex2Begin(long i)
    {
        index2Begin = i;
        index2End   = index2Begin + (index2Size - 1);
    }

#endif	

   inline double& operator()(long i1, long i2, long i3,long i4) const
   {
        return *(dataPtr + 
                  (i4 - index4Begin)+   
                    index4Size*(
                     (i3 - index3Begin) 
                        + index3Size*((i2 - index2Begin) + index2Size*(i1 - index1Begin))));
    }

    inline double* getDataPointer()
    {
		return dataPtr;
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

    inline void setIndex4Begin(long i)
    {
        index4Begin = i;
        index4End   = index4Begin + (index4Size - 1);
    }

    inline long getIndex4Begin() const
    {
        return index4Begin;
    }

    inline long getIndex4End() const
    {
        return index4End;
    }

    inline  long getIndex4Size() const
    {
        return index4Size;
    }

    inline long bufLength() const
    {
        return index1Size*index2Size*index3Size*index4Size;
    }

    inline void operator=(const doublearray4d& d)
    {
        long i;
    
        if(this->bufLength()== 0)
    		this->initialize(d.index1Size, d.index2Size, d.index3Size,d.index4Size);
    
        assert(this->bufLength()==this->bufLength());//if size not equal, warning will raise up;
    
        for(i = 0; i < d.index1Size*d.index2Size*d.index3Size*d.index4Size; i++)
            dataPtr[i] = d.dataPtr[i];   
    
    }
    
    inline void setToValue(double& val)
    {
    	long i;
    	for(i = 0; i < bufLength(); i++)
    	{
    		dataPtr[i] =  val;
    	}
    }
    
    
    inline void addValue(const double& val)
    {
        long i;
        
        for(i = 0; i < bufLength(); i++)
        {
    	    dataPtr[i] += val;
        }
    }

};

#endif
