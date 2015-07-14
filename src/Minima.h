//
//  RealVector.h
//  testProject
//
//  Created by Sören Lennart Berg on 10/2/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#ifndef __MINIMA__
#define __MINIMA__

#include "Polynomial3.h"
#include "SylvesterMatrix.h"
#include "helper.h"

#include <pari/pari.h>

class Minima
{
public:
    Minima();
    
    
    void MinimizeDeterminant(Polynomial3& polynomial);
    void FindAffineVarieties(Polynomial3* polynomials, int numPolynomials);
    
    // conversion of Polynomials <-> PARI types
//    static GEN getPARIPolynomial(Polynomial3 p);
//    static void getPolynomial3fromPARI(Polynomial3& result, GEN polynomial);
//    static void getPolynomial3fromPARI(Polynomial3& result, GEN polynomial, int v[3]);
    
private:
    double eps;
//    GEN userVariables;
//    GEN monomials;  // monimial[i] is vector of the monomials of the i-th variable of degree 0 to VAR_DEGREE
};
#endif /* defined(__MINIMA__) */
