//
//  Simplex.h
//  testProject
//
//  Created by Sören Lennart Berg on 10/2/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#ifndef __testProject__Simplex__
#define __testProject__Simplex__

#include <iostream>
#include "Matrix.h"
#include "Vector.h"

class LPSolver
{
public:
    LPSolver(const RealMatrix& A, const RealVector& b, const RealVector& c);
    //LPSolver(const RealMatrix& A, const RealVector& b, const RealVector& c, double eps);
    ~LPSolver();
    
   // double& epsilon(void) { return eps; }
    void epsilon(double eps_) { eps = eps_; }
    double epsilon(void) const { return eps; }
    
    void Pivot(int r, int s);
    bool Simplex(int phase);
    double Solve(RealVector &x);
    
private:
    int m, n;
    RealMatrix D;
    double eps;
    //RealVector b,c;
    
    int* B,*N;
};

#endif /* defined(__testProject__Simplex__) */
