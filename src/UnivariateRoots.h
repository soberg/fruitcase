//
//  UnivariateRoots.h
//  densestlatticepackings
//
//  Created by Sören Lennart Berg on 10/03/14.
//  Copyright (c) 2014 Sören Lennart Berg. All rights reserved.
//

#ifndef __densestlatticepackings__UnivariateRoots__
#define __densestlatticepackings__UnivariateRoots__

#include <iostream>
#include <pari/pari.h>
#include "Minima.h"
#include "LinearSolution3.h"

const long DEFAULT_PARI_PRECISION = DEFAULTPREC; // ~19 decimal digits

class UnivariateRoots
{
public:
    UnivariateRoots();
    void Roots(const Polynomial3& f, Monomial3 var, LinearSolution3& solutions);
    
    GEN getPARIPolynomial(const Polynomial3& p);
    GEN getLinearPARIPolynomial(const Polynomial3& p, Monomial3 var);
    void getPolynomial3fromPARI(Polynomial3& result, GEN polynomial);
    void getPolynomial3fromPARI(Polynomial3& result, GEN polynomial, int v[3]);
private:
    GEN userVariables;
    GEN monomials;  // monimial[i] is vector of the monomials of the i-th variable of degree 0 to VAR_DEGREE
    static double PARIscalarToDouble(GEN g);
};

#endif /* defined(__densestlatticepackings__UnivariateRoots__) */
