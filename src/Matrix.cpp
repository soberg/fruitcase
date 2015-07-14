//
//  RealMatrix.cpp
//  testProject
//
//  Created by Sören Lennart Berg on 10/2/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#include "Matrix.h"

RealMatrix::RealMatrix()
{
    m = n = 0;
    entries = 0;
}

RealMatrix::RealMatrix(const RealMatrix& A)
{
    delete[] entries;
    m = A.rows();
    n = A.columns();
    entries = new double[m * n];
    
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            this->operator()(i,j) = A(i,j);
}

RealMatrix::~RealMatrix()
{
    if(entries)
        delete[] entries;
}

RealMatrix::RealMatrix(int _m)
{
    m = n = _m;
    entries = new double[m * n];
    zero();
}

RealMatrix::RealMatrix(int _m, int _n)
{
    m = _m;
    n = _n;
    //mat = allocateRealMatrix(m, n);
    entries = new double[m * n];

    zero();
}

RealMatrix& RealMatrix::operator= (const RealMatrix& other)
{
    delete[] entries;
    
    m = other.rows();
    n = other.columns();
    
    entries = new double[m*n];
    
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            operator()(i,j) = other(i,j);
    
    return *this;
}

void RealMatrix::zero(void)
{
    set_all(0.0);
}

void RealMatrix::one(void)
{
    set_all(1.0);
}

void RealMatrix::set_all(double d)
{
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            operator()(i,j) = d;
}

void RealMatrix::identity(void)
{
    set_all(0.0);
    for(int i=0; i<n && i<m; ++i)
        operator()(i,i) = 1.0;
}

void RealMatrix::identity(int m)
{
    resize(m,m);
    identity();
}

void RealMatrix::resize(int _m, int _n)
{
    if( _m < m)
    {
        double* temp = new double[_m * n];
        
        for(int i=0; i<_m*n; ++i)
            temp[i] = entries[i];
        
        delete[] entries;
        entries = temp;
        
        m = _m;
    }
    else if( _m > m)
    {
        double* temp = new double[_m * n];
        
        for(int i=0; i<m*n; ++i)
            temp[i] = entries[i];
        for(int i=m*n; i<_m*n; ++i)
            temp[i] = 0.0;
        
        delete[] entries;
        entries = temp;
        
        m = _m;
    }
    
    if( _n < n)
    {
        double* temp = new double[m*_n];
        
        for(int i=0; i<m; i++)
            for(int j=0; j<_n; ++j)
                temp[I(i,j,m,_n)] = entries[I(i,j)];
        delete[] entries;
        entries = temp;
        n = _n;
    }
    else if( _n > n)
    {
        double* temp = new double[m*_n];
        
        for(int i=0; i<m; i++)
        {
            for(int j=0; j<n; ++j)
                temp[I(i,j,m,_n)] = entries[I(i,j)];
            for(int j=n; j<_n; ++j)
                temp[I(i,j,m,_n)] = 0.0;
        }
        delete[] entries;
        entries = temp;
        n = _n;
    }
}

int RealMatrix::I (int _m, int _n) const
{
    assert( 0 <= _m && _m <= m && 0 <= _n && _n <= n );
    return _n + n * _m;
}

int RealMatrix::I(int _m, int _n, int a, int b) const
{
    assert( 0 <= _m && _m <= a && 0 <= _n && _n <= b );
    return _n + n * _m;
}


