//
//  pk_Array.h
//  PumpKin
//
//  Version 1.0
//
//  Created by Aram H. Markosyan on 9/21/13.
//  Copyright (c) 2013 Aram H. Markosyan. All rights reserved.
//
//
// This file orginally was written by Jorge Balbas and Eitan Tadmor, Nov. 2004. 
// The present version extends the original functionality nad now it is part of 
// PumpKin (see http://www.pumpkin-tool.org). Please see the Copyright_and_License 
// file for the copyright notice, disclaimer, contact information and the License.
//
// PumpKin is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
// Public License for more details.
//
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, contact Aram H. Markosyan at 
// aram.math@gmail.com or write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef __PumpKin__pk_Array__
#define __PumpKin__pk_Array__

#include<iostream>
#include "stdlib.h"

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
	double&  operator()(long i1);
	const double&  operator()(long i1) const;
	double* getDataPointer();
    
    friend doublearray1d operator+(const doublearray1d& A1,
                                   const doublearray1d& A2);
    
    friend doublearray1d operator-(const doublearray1d& A1,
                                   const doublearray1d& A2);
    
    friend double operator*(const doublearray1d& A1,
                            const doublearray1d& A2);
    
    friend doublearray1d operator*(const double&        a1,
                                   const doublearray1d& A2);
    
    friend doublearray1d operator*(const doublearray1d& A1,
                                   const double&        a2);
    
	
	void setIndex1Begin(long i);
	long getIndex1Begin() const;
	long getIndex1End() const;
	long getIndex1Size() const;
    
	void resize(long newSize);
	void operator=(const doublearray1d& d);
	
	void setToValue(double d);
	void addValue(double d);
};

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
	double&  operator()(long i1, long i2);
	const double& operator()(long i1, long i2) const;
	double* getDataPointer();
    
    friend doublearray2d operator+(const doublearray2d& A1,
                                   const doublearray2d& A2);
    
    friend doublearray2d operator-(const doublearray2d& A1,
                                   const doublearray2d& A2);
    
    friend doublearray2d operator*(const double&        a1,
                                   const doublearray2d& A2);
    
    friend doublearray2d operator*(const doublearray2d& A1,
                                   const double&        a2);
	
	void setIndex1Begin(long i);
	long getIndex1Begin() const;
	long getIndex1End() const;
	long getIndex1Size() const;
	
	void setIndex2Begin(long i);
	long getIndex2Begin() const;
	long getIndex2End() const;
	long getIndex2Size() const;
	
	void operator=(const doublearray2d& d);
	
	void setToValue(double d);
	void addValue(double d);
};

#endif /* defined(__PumpKin__pk_Array__) */
