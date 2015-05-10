//
//  pk_Array.cpp
//  PumpKin
//
//  Version 1.4
//
//  Created by Aram H. Markosyan on 9/21/13.
//  Copyright (c) 2013 - 2015 Aram H. Markosyan. All rights reserved.
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
// This file contains the source codes of the classes doublearray1d and doublearray2d,
// which are basically the abstractions of an array (1D) and matrix (2D) containing
// a values of 'double' type. These classes provide flexibility and extra
// functionality to the work with a data in the form of an array or matrix.
// Functionality includes the features like array summation, subtractions, scalar product
// etc. This functionality is achieved by overloading the basic operators like +,-,* etc.
//


#include "pk_Array.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include <cstdio>

using namespace std;

// overloading the operator +
doublearray1d operator+(const doublearray1d& A1,
                        const doublearray1d& A2)
{
	if (A1.getIndex1Size() != A2.getIndex1Size()) {
		throw "Array size mismatch!";
	}
	doublearray1d temp(A1.getIndex1Size());
	for (int i = 0; i < A1.getIndex1Size(); i++) {
		temp(i) = A1(i) + A2(i);
	}
	return temp;
}

// overloading the operator -
doublearray1d operator-(const doublearray1d& A1,
                        const doublearray1d& A2)
{
	if (A1.getIndex1Size() != A2.getIndex1Size()) {
		throw "Array size mismatch!";
	}
	doublearray1d temp(A1.getIndex1Size());
	for (int i = 0; i < A1.getIndex1Size(); i++) {
		temp(i) = A1(i) - A2(i);
	}
	return temp;
}

// overloading the operator *
double operator*(const doublearray1d& A1,
                 const doublearray1d& A2)
{
	if (A1.getIndex1Size() != A2.getIndex1Size()) {
		throw "Array size mismatch!";
	}
	double temp = 0;
	for (int i = 0; i < A1.getIndex1Size(); i++) {
		temp = temp + A1(i) * A2(i);
	}
	return temp;
}

// overloading the operator *
doublearray1d operator*(const double&        a1,
                        const doublearray1d& A2)
{
	doublearray1d temp(A2.getIndex1Size());
	for (int i = 0; i < A2.getIndex1Size(); i++) {
		temp(i) = a1 * A2(i);
	}
	return temp;
}

// overloading the operator *
doublearray1d operator*(const doublearray1d& A1,
                        const double&        a2)
{
	doublearray1d temp(A1.getIndex1Size());
	for (int i = 0; i < A1.getIndex1Size(); i++) {
		temp(i) = a2 * A1(i);
	}
	return temp;
}

// default constructor
doublearray1d::doublearray1d()
{
	dataPtr = 0;
	internalAlloc = 0;
	index1Size = 0;
	index1Begin = 0;
	index1End = 0;
}

// constructor
doublearray1d::doublearray1d(long size)
{
	dataPtr = 0;
	internalAlloc = 0;
	initialize(size);
}

// constructor
doublearray1d::doublearray1d(const doublearray1d& d)
{
	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;
	dataPtr = new double[index1Size];
	internalAlloc = 1;

	long i;
	for (i = 0; i < index1Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

// destructor
doublearray1d::~doublearray1d()
{
	if (internalAlloc == 1)
		delete [] dataPtr;
}

// initialization of of array with size of m
void doublearray1d::initialize(long m)
{
	if (internalAlloc == 1)
	{
		if (index1Size != m)
		{
			delete [] dataPtr;
			dataPtr = new double[m];
		}
	}

	else
	{
		if (dataPtr == 0)
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
	if (internalAlloc == 1)
	{
		if (index1Size != d.index1Size)
		{
			delete [] dataPtr;
			dataPtr = new double[d.index1Size];
		}
	}

	else
	{
		if (dataPtr == 0)
		{
			dataPtr = new double[d.index1Size];
			internalAlloc = 1;
		}
	}

	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;

	long i;

	for (i = 0; i < index1Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

double& doublearray1d::operator()(long i1)
{
	if (i1 > index1Size - 1)
		throw "Out of bounds in array!";

	return *(dataPtr +  (i1 - index1Begin));
}

const double& doublearray1d::operator()(long i1) const
{
	if (i1 > index1Size - 1)
		throw "Out of bounds in array!";

	return *(dataPtr +  (i1 - index1Begin));
}

double* doublearray1d::getDataPointer()
{
	return dataPtr;
}

void doublearray1d::setIndex1Begin(long i)
{
	index1Begin = i;
	index1End   = index1Begin + (index1Size - 1);
}

long doublearray1d::getIndex1Begin() const
{
	return index1Begin;
}

long doublearray1d::getIndex1End() const
{
	return index1End;
}

long doublearray1d::getIndex1Size() const
{
	return index1Size;
}

// resizing the array
void doublearray1d::resize(long newSize)
{
	long i;
	double*  newDataPtr = new double[newSize];
	double*  tmpDataPtr;

	if (newSize > index1Size)
	{
		for (i = 0; i < index1Size; i++)
			newDataPtr[i] = dataPtr[i];
	}

	else
	{
		for (i = 0; i < newSize; i++)
			newDataPtr[i] = dataPtr[i];
	}

	index1Size = newSize;
	tmpDataPtr = dataPtr;
	dataPtr    = newDataPtr;

	if (internalAlloc == 1) delete [] tmpDataPtr;
	internalAlloc = 1;

	index1End = index1Begin + (index1Size - 1);
}

void doublearray1d::operator=(const doublearray1d& d)
{
	initialize(d.index1Size);

	long i;
	for (i = 0; i < d.index1Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

// setting all elements to the constant value
void doublearray1d::setToValue(double val)
{
	long i;
	for (i = 0; i < index1Size; i++)
	{
		dataPtr[i] =  val;
	}
}

// adding val to the all elements in the array
void doublearray1d::addValue(double val)
{
	long i;

	for (i = 0; i < index1Size; i++)
	{
		dataPtr[i] += val;
	}
}

ostream& operator << (ostream& out_stream, const doublearray1d& d)
{
	double outvalue;

	cout.setf(ios::scientific);
	cout.setf(ios::floatfield);
	cout.precision(16);

	long i;
	for (i = d.index1Begin; i <= d.index1End; i++)
	{
		outvalue = d.dataPtr[i];

		if (outvalue < 0 )
			out_stream << setprecision(16) << outvalue << " ";
		else
			out_stream << " " << setprecision(16) << outvalue << " ";
	}

	return out_stream;
}

// default constructor
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

// constructor
doublearray2d::doublearray2d(long size1, long size2)
{
	dataPtr = 0;
	internalAlloc = 0;
	initialize(size1, size2);
}

// constructor
doublearray2d::doublearray2d(const doublearray2d& d)
{
	index1Size = d.index1Size;
	index1Begin = d.index1Begin;
	index1End = d.index1End;
	index2Size = d.index2Size;
	index2Begin = d.index2Begin;
	index2End = d.index2End;

	dataPtr = new double[index1Size * index2Size];
	internalAlloc = 1;

	long i;
	for (i = 0; i < index1Size * index2Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

// destructor
doublearray2d::~doublearray2d()
{
	if (internalAlloc == 1)
		delete [] dataPtr;
}

void doublearray2d::initialize(long size1, long size2)
{
	if (internalAlloc == 1)
	{
		if ((index1Size != size1) || (index2Size != size2))
		{
			delete [] dataPtr;
			dataPtr = new double[size1 * size2];
		}
	}

	else
	{
		if (dataPtr == 0)
		{
			dataPtr = new double[size1 * size2];
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
	if (internalAlloc == 1)
	{
		if ((index1Size != d.index1Size) || (index2Size != d.index2Size))
		{
			delete [] dataPtr;
			dataPtr = new double[d.index1Size * d.index2Size];
		}
	}

	else
	{
		if (dataPtr == 0)
		{
			dataPtr = new double[d.index1Size * d.index2Size];
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

	for (i = 0; i < index1Size * index2Size; i++)
	{
		dataPtr[i] = d.dataPtr[i];
	}
}

doublearray2d operator+(const doublearray2d& A1,
                        const doublearray2d& A2)
{
	if (A1.getIndex1Size() != A2.getIndex1Size()) {
		throw "Array size mismatch!";
	}
	if (A1.getIndex2Size() != A2.getIndex2Size()) {
		throw "Array size mismatch!";
	}

	doublearray2d temp(A1.getIndex1Size(), A1.getIndex2Size());
	for (int i = 0; i < A1.getIndex1Size(); i++) {
		for (int j = 0; j < A1.getIndex2Size(); j++) {
			temp(i, j) = A1(i, j) + A2(i, j);
		}
	}
	return temp;
}

doublearray2d operator-(const doublearray2d& A1,
                        const doublearray2d& A2)
{
	if (A1.getIndex1Size() != A2.getIndex1Size()) {
		throw "Array size mismatch!";
	}
	if (A1.getIndex2Size() != A2.getIndex2Size()) {
		throw "Array size mismatch!";
	}

	doublearray2d temp(A1.getIndex1Size(), A1.getIndex2Size());
	for (int i = 0; i < A1.getIndex1Size(); i++) {
		for (int j = 0; j < A1.getIndex2Size(); j++) {
			temp(i, j) = A1(i, j) - A2(i, j);
		}
	}
	return temp;
}

doublearray2d operator*(const double&        a1,
                        const doublearray2d& A2)
{
	doublearray2d temp(A2.getIndex1Size(), A2.getIndex2Size());
	for (int i = 0; i < A2.getIndex1Size(); i++) {
		for (int j = 0; j < A2.getIndex2Size(); j++) {
			temp(i, j) = A2(i, j) * a1;
		}
	}
	return temp;
}

doublearray2d operator*(const doublearray2d& A1,
                        const double&        a2)
{
	doublearray2d temp(A1.getIndex1Size(), A1.getIndex2Size());
	for (int i = 0; i < A1.getIndex1Size(); i++) {
		for (int j = 0; j < A1.getIndex2Size(); j++) {
			temp(i, j) = A1(i, j) * a2;
		}
	}
	return temp;
}

double& doublearray2d::operator()(long i1, long i2)
{
	if (i1 > index1Size - 1)
		throw "Out of bounds in array!";
	if (i2 > index2Size - 1)
		throw "Out of bounds in array!";

	return *(dataPtr + (i2 - index2Begin) + (i1 - index1Begin) * index2Size);
}

const double& doublearray2d::operator()(long i1, long i2) const
{
	if (i1 > index1Size - 1)
		throw "Out of bounds in array!";
	if (i2 > index2Size - 1)
		throw "Out of bounds in array!";

	return *(dataPtr + (i2 - index2Begin) + (i1 - index1Begin) * index2Size);
}

double* doublearray2d::getDataPointer()
{
	return dataPtr;
}

void doublearray2d::setIndex1Begin(long i)
{
	index1Begin = i;
	index1End = index1Begin + (index1Size - 1);
}

long doublearray2d::getIndex1Begin() const
{
	return index1Begin;
}

long doublearray2d::getIndex1End() const
{
	return index1End;
}

long doublearray2d::getIndex1Size() const
{
	return index1Size;
}

void doublearray2d::setIndex2Begin(long i)
{
	index2Begin = i;
	index2End = index2Begin + (index2Size - 1);
}

long doublearray2d::getIndex2Begin() const
{
	return index2Begin;
}

long doublearray2d::getIndex2End() const
{
	return index2End;
}

long doublearray2d::getIndex2Size() const
{
	return index2Size;
}

void doublearray2d::operator=(const doublearray2d& d)
{
	if (index1Size * index2Size == 0)
		initialize(d.index1Size, d.index2Size);

	long i;

	for (i = 0; i < d.index1Size * d.index2Size; i++)
		dataPtr[i] = d.dataPtr[i];
}

void doublearray2d::setToValue(double val)
{
	long i;
	for (i = 0; i < index1Size * index2Size; i++)
	{
		dataPtr[i] =  val;
	}
}


void doublearray2d::addValue(double val)
{
	long i;

	for (i = 0; i < index1Size * index2Size; i++)
	{
		dataPtr[i] += val;
	}
}

ostream& operator<< (ostream& out_stream, const doublearray2d& d)
{
	double outvalue;

	cout.setf(ios::scientific, ios::floatfield);
	cout.precision(16);

	long i, j, k;
	for (i = d.index1Begin; i <= d.index1End; i++)
	{
		for (j = d.index2Begin; j <= d.index2End; j++)
		{
			k = (j - d.index2Begin) + (i - d.index1Begin) * d.index2Size;
			outvalue = d.dataPtr[k] ;

			if (outvalue < 0 )
				out_stream << setprecision(16) << outvalue << " ";
			else
				out_stream << " " << setprecision(16) << outvalue << " ";
		}

		out_stream << endl;
	}

	return out_stream;
}

