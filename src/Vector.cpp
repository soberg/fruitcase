//
//  RealVector.cpp
//  testProject
//
//  Created by Sören Lennart Berg on 10/2/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#include "Vector.h"

RealVector::RealVector()
{
    n = 0;
    vec = 0;
}

RealVector::RealVector(const RealVector& v)
{
    n = v.size();
    vec = new double[n];
    for(int i=0; i<n; ++i)
        vec[i] = v(i);
}

RealVector::RealVector(int _n)
{
    n = _n;
    double* t = new double[n];
    vec = new double[n];
    for(int i=0; i<n; ++i)
        vec[i] = 0;
}

RealVector::~RealVector()
{
    if(vec)
    {
        delete[] vec;
        vec = 0;
    }
}

RealVector& RealVector::operator= (const RealVector& other)
{
    if(vec) delete[] vec;
    
    n = other.size();
    vec = new double[n];
    for(int i=0; i<n; ++i)
        operator()(i) = other(i);
    
    return *this;
}

void RealVector::resize(int _n)
{
    if( _n < n )
    {
        double* temp = new double[_n];
        for(int i=0; i<_n; ++i)
            temp[i] = vec[i];
        delete[] vec;
        vec = temp;
    }
    else if( _n > n )
    {
        double* temp = new double[_n];
        for(int i=0; i<n; ++i)
            temp[i] = vec[i];
        for(int i=n; i<_n; ++i)
            temp[i] = 0.0;
        delete[] vec;
        vec = temp;
    }
    n = _n;
}

void RealVector::zero(void)
{
    set_all(0.0);
}

void RealVector::one(void)
{
    set_all(1.0);
}

void RealVector::set_all(double d)
{
    for(int i=0; i<n; ++i)
        vec[i] = d;
}
