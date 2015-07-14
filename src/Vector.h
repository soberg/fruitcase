//
//  RealVector.h
//  testProject
//
//  Created by Sören Lennart Berg on 10/2/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#ifndef __testProject__RealVector__
#define __testProject__RealVector__

class RealVector
{
public:
    RealVector();
    RealVector(const RealVector& v);
    RealVector(int _n);
    ~RealVector();
    
    RealVector& operator= (const RealVector& other);
    
    double& operator()(int i) { return vec[i]; }
    double  operator()(int i) const { return vec[i]; }
    
    void resize(int _n);
    
    void zero(void);
    void one(void);
    void set_all(double d);
    
    int size(void) const { return n; }
    
private:
    int n;
    double* vec;
};
#endif /* defined(__testProject__RealVector__) */
