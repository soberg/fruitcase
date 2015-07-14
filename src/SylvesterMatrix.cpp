//
//  SylvesterMatrix.cpp
//  densestlatticepackings
//
//  Created by Sören Lennart Berg on 06/03/14.
//  Copyright (c) 2014 Sören Lennart Berg. All rights reserved.
//

#include "SylvesterMatrix.h"

SylvesterMatrix::SylvesterMatrix(const Polynomial3& f, const Polynomial3& g,const Monomial3 variable)
{
    var = variable;
    
    // calculate degree regarding the desired variable
    varDegree[0] = 0;//f.degree(var);
    varDegree[1] = 0;//g.degree(var);
    
    for(int i=0; i < DEFAULT_POLYNOMIAL_LENGTH; ++i)
    {
        // if e.g. var=x then we transform the polynomials
        // in the form p_n x^n + ... + p_1 x + p_0 (p_i is poly. in y,z)
        // and save the p_i's in ascending order
        // in the array `polynomials'
        if( f.polynomialCoefficient(var, i, polynomials[0][i]))
            varDegree[0] = i;
        
        if( g.polynomialCoefficient(var, i, polynomials[1][i]))
            varDegree[1] = i;
    }
}

int SylvesterMatrix::size(void) const
{
    return (varDegree[0]) + (varDegree[1]);
}

bool SylvesterMatrix::isZero() const
{
    for(int i=0; i<2; ++i)
    {
        for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
        {
            if( !polynomials[i][k].isZero() )
                return false;
        }
    }
    return true;
}

bool SylvesterMatrix::isZero(int i, int j) const
{
    if( i < varDegree[1] )
    {
        int k = varDegree[0] - (j-i);
        if( k <= varDegree[0] && k >=0)
            return polynomials[0][k].isZero();
        else
            return true;
    }
    else
    {
        int k = varDegree[1] - (j-(i-varDegree[1]));
        if( k <= varDegree[1] && k >=0)
            return polynomials[1][k].isZero();
        else
            return true;
    }
}

const Polynomial3& SylvesterMatrix::operator()(int i, int j) const
{
    // TODO assert
    if( i < varDegree[1] )
    {
        int k = varDegree[0] - (j-i);
        if( k <= varDegree[0] && k >=0)
            return polynomials[0][k];
        else
            return nullPolynomial;
    }
    else
    {
        int k = varDegree[1] - (j-(i-varDegree[1]));
        if( k <= varDegree[1] && k >=0)
            return polynomials[1][k];
        else
            return nullPolynomial;
    }
}

Polynomial3 SylvesterMatrix::determinant() const
{    
    bool rows[SYLV_MAT_MAX_SIZE];
    
    std::fill_n(rows, SYLV_MAT_MAX_SIZE, true);
    
    return det_minor(0, rows);
}

Polynomial3 SylvesterMatrix::det_minor(int column, bool rows[SYLV_MAT_MAX_SIZE]) const
{
    int s = size();
    Polynomial3 result;
    
    // end of recursion?
    if(s-1 == column)
    {
        for(int i=0; i<s; ++i)
        {
            if( rows[i])
                return (*this)(i,column);
            // TODO error, <- reaching this point should never happen!
        }
    }
    
    // further recursion
    int sign = -1;
    for(int i=0; i<s; ++i)
    {
        if( !rows[i])
            continue;
        
        sign *= -1;
        
        if( isZero(i,column) )
            continue;
        
        rows[i] = false; // cancel row out
//        Polynomial3 temp = det_minor(column+1, rows);
        result += sign*(*this)(i,column)*det_minor(column+1, rows);
        
//        std::cout << result.print() << std::endl;
//        std::cout << temp.print() << std::endl;
//        std::cout << (*this)(i,column).print() << std::endl;
//        Polynomial3 bla =(*this)(i,column)*det_minor(column+1, rows);
//        std::cout << bla.print() << std::endl;
        rows[i] = true; // take row back in
    }
    
    return result;
}

void SylvesterMatrix::print(void) const
{
    int s = size();
    std::cout << "size " << s << std::endl;
    for(int i=0; i< s; ++i)
    {
        for(int j=0; j<s; ++j)
        {
            std::cout << (*this)(i,j).print() << "\t";
        }
        std::cout << std::endl;
    }
}