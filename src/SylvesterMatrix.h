//
//  SylvesterMatrix.h
//  densestlatticepackings
//
//  Created by Sören Lennart Berg on 06/03/14.
//  Copyright (c) 2014 Sören Lennart Berg. All rights reserved.
//

#ifndef __densestlatticepackings__SylvesterMatrix__
#define __densestlatticepackings__SylvesterMatrix__

#include <iostream>
#include <algorithm>
#include "Polynomial3.h"

const int SYLV_MAT_MAX_SIZE = 2*DEFAULT_POLYNOMIAL_LENGTH;

class SylvesterMatrix
{
public:
    SylvesterMatrix(const Polynomial3& f, const Polynomial3& g,
                    const Monomial3 variable);
    
    int size(void) const;
    
    bool isZero() const;
    bool isZero(int i, int j) const;
    const Polynomial3& operator()(int i, int j) const;
    
    Polynomial3 determinant() const;
    
    void print(void) const;
    
private:
    Polynomial3 det_minor(int column, bool rows[SYLV_MAT_MAX_SIZE]) const;
    
    Monomial3 var;
    int varDegree[2];
    
    Polynomial3 polynomials[2][DEFAULT_POLYNOMIAL_LENGTH];
    Polynomial3 nullPolynomial;
    
    
    //Polynomial3 elements_g[DEFAULT_POLYNOMIAL_VAR_DEGREE];
};

#endif /* defined(__densestlatticepackings__SylvesterMatrix__) */
