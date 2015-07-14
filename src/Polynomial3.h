//
//  Polynomial3.h
//  testProject
//
//  Created by Sören Lennart Berg on 10/9/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#ifndef __testProject__Polynomial3__
#define __testProject__Polynomial3__

#include <string>
#include <sstream>      // std::ostringstream
#include <iostream>
#include <assert.h>
#include <limits>

#include "helper.h"

enum Monomial3 { VAR_X=0, VAR_Y=1, VAR_Z=2, VAR_CONST=3 };

const double    DEFAULT_POLYNOMIAL_EPSILON      = 10E-07;
const int       DEFAULT_POLYNOMIAL_VAR_DEGREE   = 5;
const int       DEFAULT_POLYNOMIAL_LENGTH   = DEFAULT_POLYNOMIAL_VAR_DEGREE+1;
const int       DEFAULT_POLYNOMIAL_NUM_COEFF =
DEFAULT_POLYNOMIAL_LENGTH*DEFAULT_POLYNOMIAL_LENGTH*DEFAULT_POLYNOMIAL_LENGTH;

const Monomial3 POL3_VARS[4] = { VAR_X, VAR_Y, VAR_Z, VAR_CONST };

class Polynomial3
{
public:
    Polynomial3();
    Polynomial3(double _eps);
    
    Polynomial3& operator= (const Polynomial3& other);
    
    bool operator== (const Polynomial3& a) const;
    bool operator!= (const Polynomial3& a) const;
    bool operator< (const Polynomial3& a) const;
    bool operator<= (const Polynomial3& a) const;
    bool operator> (const Polynomial3& a) const;
    bool operator>= (const Polynomial3& a) const;
     
    Polynomial3 operator+ (const Polynomial3& a) const;
    Polynomial3 operator- (const Polynomial3& a) const;
    Polynomial3 operator* (const Polynomial3& a) const;
    Polynomial3 operator* (const double d) const;
    //Polynomial3 operator/ (const Polynomial3& a) const;
    
    Polynomial3& operator+=(const Polynomial3& a);
    Polynomial3& operator-=(const Polynomial3& a);
    Polynomial3& operator*=(const Polynomial3& a);
     
    Polynomial3 operator- () const;
    
    Polynomial3 timesMonomial(int ex, int ey, int ez, double c) const;
    
    void setAll(double d);
    
//    Polynomial3 shift(int i, int j, int k) const;
    
    
    bool isZero() const;
    bool isZero(int i, int j, int k) const;
    double SumOfCoefficients() const;
    
    int totalDegree() const;
    int degree(Monomial3 var) const;
//    int leadCoefficient(int var) const;
    bool lexicLeadTerm(int& i, int& j, int& k) const;
    
    double operator()(int i, int j, int k) const;
    double& operator()(int i, int j, int k);
    
    Polynomial3 substitute(double value, Monomial3 var) const;
    Polynomial3 substitute(const Polynomial3& f, Monomial3 var) const;
    
    Polynomial3 power(int e) const;
    
    bool polynomialCoefficient(Monomial3 var, int e, Polynomial3& result) const;
    
    void straighten();
    
    void derive(Polynomial3& result, Monomial3 var) const;
    
    double coefficients [DEFAULT_POLYNOMIAL_LENGTH]
                        [DEFAULT_POLYNOMIAL_LENGTH]
                        [DEFAULT_POLYNOMIAL_LENGTH];
    
    
    double eps;
    
    std::string print() const;
    

    double& var_permute(int i, Monomial3 var, int j, int k);
    double var_permute(int i, Monomial3 var, int j, int k) const;
    
private:
    Polynomial3 timesMonomial(int e1, Monomial3 var, int e2, int e3, double c) const;
};

Polynomial3 operator* (const double d, const Polynomial3& a);

Polynomial3 resultant( const Polynomial3& p, const Polynomial3& q);

Polynomial3 pdiv(const Polynomial3& p, const Polynomial3& q, Polynomial3& r, int var);
Polynomial3 pdiv2(const Polynomial3& p, const Polynomial3& q, Polynomial3& r, int var);
//Polynomial3 pdivYZ(const Polynomial3& p, const Polynomial3& q, Polynomial3& r);
//Polynomial3 pdivZ(const Polynomial3& p, const Polynomial3& q, Polynomial3& r);

bool poldiv(const Polynomial3& f, const Polynomial3& g, Polynomial3& q, Polynomial3& r);
bool poldiv(const Polynomial3& f, const Polynomial3& g, Polynomial3& q);
bool poldiv(const Polynomial3& f, const Polynomial3& g);


void SylvesterMatrixElement(const Polynomial3& f,
                            const Polynomial3& g,
                            const int varDegree_f,
                            const int varDegree_g,
                            const int variable,
                            int& i, int& j, int& k,
                            double& coeff);

Monomial3 remainingVariable(Monomial3 var, int i);

#endif /* defined(__testProject__Polynomial3__) */
