//
//  LinearSolution3.h
//  densestlatticepackings
//
//  Created by Sören Lennart Berg on 10/03/14.
//  Copyright (c) 2014 Sören Lennart Berg. All rights reserved.
//

#ifndef __densestlatticepackings__LinearSolution3__
#define __densestlatticepackings__LinearSolution3__

#include <iostream>
#include <vector>
#include <assert.h>
#include "Polynomial3.h"

// three polynomials(coeff[0][]...coeff[2][])
// with linear coefficients for x,y,z and const
class LinearSolutionAtom3
{
public:
    LinearSolutionAtom3();
    LinearSolutionAtom3(double x, double y, double z);
    LinearSolutionAtom3(Monomial3 var, double d);
    
    void setAll(double d);
    void print() const;
//private:
    double coefficients[3][4];
};

class LinearSolution3
{
public:
    LinearSolution3();
    
    void addSolution(double x, double y, double z);
    void addSolution(Monomial3 var, double d);
    
    void print() const;
    
private:
    std::vector<LinearSolutionAtom3> solutions;
};

#endif /* defined(__densestlatticepackings__LinearSolution3__) */
