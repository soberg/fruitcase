//
//  Polynomial3.cpp
//  testProject
//
//  Created by Sören Lennart Berg on 10/9/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#include "Polynomial3.h"

Polynomial3::Polynomial3()
{
    setAll(0.0);
    eps = DEFAULT_POLYNOMIAL_EPSILON;
}

Polynomial3::Polynomial3(double _eps)
{
    setAll(0.0);
    eps = _eps;
}

void Polynomial3::setAll(double d)
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                coefficients[i][j][k] = d;
}

Polynomial3& Polynomial3::operator= (const Polynomial3& other)
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                coefficients[i][j][k] = other.coefficients[i][j][k];
    return *this;
}

bool Polynomial3::operator== (const Polynomial3& a) const
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
            {
                int diff = coefficients[i][j][k] - a.coefficients[i][j][k];
                if( diff > eps || diff < -eps )
                    return false;
            }
    return true;
}

bool Polynomial3::operator!= (const Polynomial3& a) const
{
    return !operator==(a);
}

bool Polynomial3::operator< (const Polynomial3& a) const
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                if( coefficients[i][j][k] >= a.coefficients[i][j][k] )
                    return false;
    
    return true;
}

bool Polynomial3::operator<= (const Polynomial3& a) const
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                if( coefficients[i][j][k] > a.coefficients[i][j][k] )
                    return false;
    
    return true;
}

bool Polynomial3::operator> (const Polynomial3& a) const
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                if( coefficients[i][j][k] <= a.coefficients[i][j][k] )
                    return false;
    
    return true;
}

bool Polynomial3::operator>= (const Polynomial3& a) const
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                if( coefficients[i][j][k] < a.coefficients[i][j][k] )
                    return false;
    
    return true;
}

Polynomial3 Polynomial3::operator+ (const Polynomial3& a) const
{
    Polynomial3 result;
    
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                result.coefficients[i][j][k] = coefficients[i][j][k] + a.coefficients[i][j][k];
    
    return result;
}

Polynomial3 Polynomial3::operator- (const Polynomial3& a) const
{
    Polynomial3 result;
    
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                result.coefficients[i][j][k] = coefficients[i][j][k] - a.coefficients[i][j][k];
    
    return result;
}

Polynomial3 Polynomial3::operator* (const Polynomial3& a) const
{
    // TODO UEBERLAUF ABFANGEN!!!
    
    Polynomial3 result;
    
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
            {
                //result.coefficients[i][j][k]
                for(int px = 0; px<=i; ++px)
                    for(int py=0; py<=j; ++py)
                        for(int pz=0; pz<=k; ++pz)
                            result.coefficients[i][j][k] += coefficients[px][py][pz] * a.coefficients[i-px][j-py][k-pz];
            }
    
    return result;
}

Polynomial3 Polynomial3::operator* (const double d) const
{
    Polynomial3 result;
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                result(i,j,k) = coefficients[i][j][k] * d;
    return result;
}

// multiplies the polynomial with (c * x^ex y^ey z^ez)
Polynomial3 Polynomial3::timesMonomial(int ex, int ey, int ez, double c) const
{
    Polynomial3 result;
    for(int i=ex; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=ey; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=ez; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                result(i,j,k) = coefficients[i-ex][j-ey][k-ez] * c;
    return result;
}

Polynomial3 Polynomial3::timesMonomial(int e1, Monomial3 var, int e2, int e3, double c) const
{
    Polynomial3 result;
    
    for(int i=e1; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
    {
        for(int j=e2; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
        {
            for(int k=e3; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
            {
                switch(var)
                {
                    case(VAR_X):
                        result(i,j,k) = coefficients[i-e1][j-e2][k-e3] * c;
                        break;
                    case(VAR_Y):
                        result(j,i,k) = coefficients[j-e2][i-e1][k-e3] * c;
                        break;
                    case(VAR_Z):
                        result(j,k,i) = coefficients[j-e2][k-e3][i-e1] * c;
                        break;
                }
            }
        }
    }
    return result;
}

//Polynomial3 Polynomial3::operator/ (const Polynomial3& a) const
//{
//    
//}

Polynomial3 Polynomial3::operator- () const
{
    Polynomial3 result;
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                result.coefficients[i][j][k] = -coefficients[i][j][k];
    
    return result;
}

Polynomial3& Polynomial3::operator+=(const Polynomial3& a)
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                (*this)(i,j,k) += a(i,j,k);
    return *this;
}

Polynomial3& Polynomial3::operator-=(const Polynomial3& a)
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                (*this)(i,j,k) -= a(i,j,k);
    return *this;
}

Polynomial3& Polynomial3::operator*=(const Polynomial3& a)
{
    *this = (*this)*a;
    return *this;
}

double Polynomial3::operator()(int i, int j, int k) const
{
    return coefficients[i][j][k];
}

double& Polynomial3::operator()(int i, int j, int k)
{
    return coefficients[i][j][k];
}

// return the (polynomial) coefficient regarding var^e
// e.g. if var=x e=2 and the polynomial is 2x^2 z + 3 x^2 y + 2x -1
// the function returns 2z + 3y written to coefficient
//
bool Polynomial3::polynomialCoefficient(Monomial3 var, int e, Polynomial3& result) const
{
    //Polynomial3 result;
    bool is_zero = true;
    
    result.setAll(0.0);
    
    switch(var)
    {
        case VAR_X:
            
            for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
                for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
                {
                    if( !isZero(e,i,j) )
                    {
                        is_zero = false;
                        result(0,i,j) = (*this)(e,i,j);
                    }
                }
            break;
        case VAR_Y:
            for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
                for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
                {
                    if( !isZero(i,e,j) )
                    {
                        is_zero = false;
                        result(i,0,j) = (*this)(i,e,j);
                    }
                }
            break;
        case VAR_Z:
            for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
                for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
                {
                    if( !isZero(i,j,e) )
                    {
                        is_zero = false;
                        result(i,j,0) += (*this)(i,j,e);
                    }
                }
            break;
    }
    return !is_zero;
}

//double& Polynomial3::coefficient(Monomial3 var, int i, int k)
//{
//    
//}

// obsolete? -> timesMonomial
//Polynomial3 Polynomial3::shift(int _i, int _j, int _k) const
//{
//    // todo i,j, k >=0 asserten
//    Polynomial3 result;// = *this;
//    
//    for(int i=DEFAULT_POLYNOMIAL_LENGTH-1; i - _i >=0; --i)
//        for(int j=DEFAULT_POLYNOMIAL_LENGTH-1; j - _j >=0; --j)
//            for(int k=DEFAULT_POLYNOMIAL_LENGTH-1; k - _k >=0; --k)
//                result(i,j,k) = coefficients[i-_i][j-_j][k-_k];
//    
//    return result;
//}

bool Polynomial3::isZero() const
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                if( !isZero(i,j,k) )
                    return false;
    return true;
}

bool Polynomial3::isZero(int i, int j, int k) const
{
    return iszeroeps(coefficients[i][j][k], eps); // coefficients[i][j][k] > -eps && coefficients[i][j][k] < eps;
}

double Polynomial3::SumOfCoefficients() const
{
    double sum = 0.0;
    
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                sum += (*this)(i,j,k);
    return sum;
}

int Polynomial3::totalDegree() const
{
    int totalDeg = -std::numeric_limits<int>::max();
    
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                if( !isZero(i,j,k) && totalDeg < i+j+k)
                {
                    totalDeg = i+j+k;
                }
    
    return totalDeg;
}

int Polynomial3::degree(Monomial3 var) const
{
    switch(var)
    {
        case VAR_X:
            for(int i=DEFAULT_POLYNOMIAL_VAR_DEGREE; i >=0; --i)
                for(int j=DEFAULT_POLYNOMIAL_VAR_DEGREE; j >=0; --j)
                    for(int k=DEFAULT_POLYNOMIAL_VAR_DEGREE; k >=0; --k)
                    {
                        if( !isZero(i,j,k) )
                            return i;
                    }
            break;
        case VAR_Y:
            for(int j=DEFAULT_POLYNOMIAL_VAR_DEGREE; j >=0; --j)
                for(int i=DEFAULT_POLYNOMIAL_VAR_DEGREE; i >=0; --i)
                    for(int k=DEFAULT_POLYNOMIAL_VAR_DEGREE; k >=0; --k)
                    {
                        if( !isZero(i,j,k) )
                            return j;
                    }
            break;
        case VAR_Z:
            for(int k=DEFAULT_POLYNOMIAL_VAR_DEGREE; k >=0; --k)
                for(int j=DEFAULT_POLYNOMIAL_VAR_DEGREE; j >=0; --j)
                    for(int i=DEFAULT_POLYNOMIAL_VAR_DEGREE; i >=0; --i)
                    {
                        if( !isZero(i,j,k) )
                            return k;
                    }
    }
    return -std::numeric_limits<int>::max();
}

//int Polynomial3::leadCoefficient(int var) const
//{
//    assert( var>=0 && var < 3);
//    
//    int leadCoeffVar = DEFAULT_POLYNOMIAL_LENGTH-1;
//    
//    switch(var)
//    {
//        case 0:
//            while( isZero(leadCoeffVar,0,0) && leadCoeffVar > 0)
//                --leadCoeffVar;
//            break;
//        case 1:
//            while( isZero(0,leadCoeffVar,0) && leadCoeffVar > 0)
//                --leadCoeffVar;
//            break;
//        case 2:
//            while( isZero(0,0,leadCoeffVar) && leadCoeffVar > 0)
//                --leadCoeffVar;
//            break;
//    }
//    return leadCoeffVar;
//}

Polynomial3 Polynomial3::substitute(double value, Monomial3 var) const
{
    Polynomial3 result;
    for(int i=0; i < DEFAULT_POLYNOMIAL_LENGTH;++i)
    {
        for(int j=0; j < DEFAULT_POLYNOMIAL_LENGTH; ++j)
        {
            for(int k=0; k < DEFAULT_POLYNOMIAL_LENGTH; ++k)
            {
                result.var_permute(0,var,j,k) += int_pow(value, i)*var_permute(i,var,j,k);
            }
        }
    }
    return result;
}

Polynomial3 Polynomial3::substitute(const Polynomial3& f, Monomial3 var) const
{
    Polynomial3 result;
    Polynomial3 power;      // power of (*this) polynomial
    power(0,0,0) = 1.0;     // start with ^0
    
    for(int i=0; i < DEFAULT_POLYNOMIAL_LENGTH;++i)
    {
        for(int j=0; j < DEFAULT_POLYNOMIAL_LENGTH; ++j)
        {
            for(int k=0; k < DEFAULT_POLYNOMIAL_LENGTH; ++k)
            {
                result += power.timesMonomial(0,var,j,k, var_permute(i,var,j,k));
            }
        }
        power *= f;
    }
    
    return result;
}

Polynomial3 Polynomial3::power(int e) const
{
    Polynomial3 result;
    result(0,0,0) = 1.0;
    
    for(int i=0; i<e; ++i)
        result *= (*this);
    return result;
}

// returns the first non-zero monomial x^i y^j z^k with (i,j,k)
// lexicographically maximal
// returns true if polynomial is non zero, false otherwise
/*!
    Returns the first non-zero monomial x^i y^j z^k with (i,j,k)
    lexicographically maximal
    @param index i for variable x (written to)
    @param index j for variable y (written to)
    @param index k for variable k (written to)
    @return 1 if polynomial is non-zero, 0 otherwise
*/
bool Polynomial3::lexicLeadTerm(int& i, int& j, int& k) const
{
    for(i=DEFAULT_POLYNOMIAL_VAR_DEGREE; i >=0; --i)
        for(j=DEFAULT_POLYNOMIAL_VAR_DEGREE; j >=0; --j)
            for(k=DEFAULT_POLYNOMIAL_VAR_DEGREE; k >=0; --k)
            {
                if( !isZero(i,j,k) )
                    return true;
            }
    i = j = k = 0;
    return false;
}

// set values in (-eps,eps) to zero
void Polynomial3::straighten()
{
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                if( iszeroeps( (*this)(i,j,k), eps))
                    (*this)(i,j,k) = 0.0;
}

void Polynomial3::derive(Polynomial3& result, Monomial3 var) const
{
    // we basically decrease the index of the variables var time the corresponding the exponent
    switch(var)
    {
        case VAR_X:
            
                for(int j=0; j<DEFAULT_POLYNOMIAL_VAR_DEGREE+1; ++j)
                    for(int k=0; k<DEFAULT_POLYNOMIAL_VAR_DEGREE+1; ++k)
                    {
                        for(int i=0; i<DEFAULT_POLYNOMIAL_VAR_DEGREE; ++i)
                            result(i,j,k) = (i+1)* (*this)(i+1,j,k);
                        result(DEFAULT_POLYNOMIAL_VAR_DEGREE, j, k) = 0.0;  
                    }
            break;
        case VAR_Y:
            for(int i=0; i<DEFAULT_POLYNOMIAL_VAR_DEGREE+1; ++i)
                
                    for(int k=0; k<DEFAULT_POLYNOMIAL_VAR_DEGREE+1; ++k)
                    {
                        for(int j=0; j<DEFAULT_POLYNOMIAL_VAR_DEGREE; ++j)
                            result(i,j,k) = (j+1)* (*this)(i,j+1,k);
                        result(i, DEFAULT_POLYNOMIAL_VAR_DEGREE, k) = 0.0;
                    }
            break;
        case VAR_Z:
            for(int i=0; i<DEFAULT_POLYNOMIAL_VAR_DEGREE+1; ++i)
                for(int j=0; j<DEFAULT_POLYNOMIAL_VAR_DEGREE+1; ++j)
                {
                    for(int k=0; k<DEFAULT_POLYNOMIAL_VAR_DEGREE; ++k)
                        result(i,j,k) = (k+1)* (*this)(i,j,k+1);
                    result(i,j,DEFAULT_POLYNOMIAL_VAR_DEGREE) = 0.0;
                }
            break;
        default:
            std::cout << "error in Polynomial3::derive; invalid var number: " << var << std::endl;
            return;
            break;
    }
}

// this functions return the coefficient belonging to
// the monomial with var^i and the remaining variables(having exponent j,k)
// e.g. if var=y and  then the coefficient belonging
// to the monomial x^j y^i z^k is returned etc.
// the purpose is to avoid using switch statements in further code...
double& Polynomial3::var_permute(int i, Monomial3 var, int j, int k)
{
    switch(var)
    {
        case VAR_X:
            return (*this)(i,j,k);
        case VAR_Y:
            return (*this)(j,i,k);
        case VAR_Z:
            return (*this)(j,k,i);
        default:
            std::cout << "var_permute: invalid VAR" << std::endl;
            return (*this)(0,0,0); // TODO better error handling
    }
}

double Polynomial3::var_permute(int i, Monomial3 var, int j, int k) const
{
    switch(var)
    {
        case VAR_X:
            return (*this)(i,j,k);
        case VAR_Y:
            return (*this)(j,i,k);
        case VAR_Z:
            return (*this)(j,k,i);
        default:
            std::cout << "var_permute: invalid VAR" << std::endl;
            return 0.0; // TODO better error handling
    }
}

std::string Polynomial3::print() const
{
    std::string s;
    bool isFirst = true;
    
        for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
            for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
                for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
                    if(coefficients[i][j][k] < -eps || coefficients[i][j][k] > eps)
                    {
                        std::ostringstream strcoeff, stri, strj, strk;
                        strcoeff << coefficients[i][j][k];
                        stri << i;
                        strj << j;
                        strk << k;
                        //s += "  " + strcoeff.str() + " x^" + stri.str() + " y^" + strj.str() + " z^" + strk.str();
                        if(!isFirst) s += " + ";
                        s += strcoeff.str();
                        if(i>0)
                            s += " x^" + stri.str();
                        if(j>0)
                            s += " y^" + strj.str();
                        if(k>0)
                            s += " z^" + strk.str();
                        //s += "  ";
                        isFirst = false;
                    }
        if(s.length() == 0)
            s = "0";
    
    return s;
}

//std::string Polynomial3::tostring() const
//{
//    std::string s;
//    
//    //    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
//    //        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
//    //            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
//    //                if(coefficients[i][j][k] < -eps || coefficients[i][j][k] > eps)
//    //                   s += "  " + std::to_string(coefficients[i][j][k]) + " x^" + std::to_string(i) + " y^" + std::to_string(j) + " z^" + std::to_string(k);
//    //    if(s.length() == 0)
//    //        s = "0";
//    
//    return s;
//}

Polynomial3 operator* (const double d, const Polynomial3& a)
{
    return a*d;
}

Polynomial3 resultant( const Polynomial3& p, const Polynomial3& q)
{
    // todo: implement!
    Polynomial3 r;
    return r;
}

/*!
    Given two polynomials f and g this function calculates
    * q,r such that f=q*g+r
    @param polynomial f (const)
    @param polynomial g (const)
    @param polynomial q (written to)
    @param polynomial r (written to)
    @return 0 r!=0  and 1 otherwise(g divides f).
*/
bool poldiv(const Polynomial3& f, const Polynomial3& g, Polynomial3& q, Polynomial3& r)
{
    // Algorithm is as follows:
    // 1. Set q=0,r=0,p=f
    // 2. Repeat (until p=0)
    //    IF ( lead term of g divides lead term of p)
    //    THEN u := (lead term of p) / (lead term of g)
    //         q := q + u
    //         p := p - u*g
    //    ELSE r := r + (lead term of p)
    //         p := p - (lead term of p)
    
    Polynomial3 p=f;
    q.setAll(0.0);
    r.setAll(0.0);
    
    bool f_divides_g = true;
    
    int gx, gy, gz;
    g.lexicLeadTerm(gx,gy,gz);
    
    int px,py,pz;
    p.lexicLeadTerm(px,py,pz);
    
    for(int i=0; i <= DEFAULT_POLYNOMIAL_NUM_COEFF; ++i) // basically this is while(true)
                                                         // for loop to avoid infinity loop!
    {
        // check if the lead term of g divides the lead term of p
        if( px>=gx && py>=gy && pz>=gz)
        {
            int ux=px-gx, uy=py-gy, uz=pz-gz;
            q(ux, uy, uz) += p(px,py,pz) / g(gx,gy,gz); // q <- q+u
            p -= g.timesMonomial(ux,uy,uz, q(ux, uy, uz) ); // p <- p - u*g
            
            p(px,py,pz) = 0.0; // numeric. stability, (should be zero anyway)
        }
        else
        {
            f_divides_g = false;
            r(px,py,pz) += p(px,py,pz); // r <- r + (lead coeff. of p)
            p(px,py,pz) = 0.0; // p <- p - (lead coeff. of p)
        }
        
        if( !p.lexicLeadTerm(px,py,pz) ) // update indices for the lead coefficient
										 // if p==0 we are done!
        {
            break;	// we're done!
        }
    }
    return f_divides_g;
}

bool poldiv(const Polynomial3& f, const Polynomial3& g, Polynomial3& q)
{
    Polynomial3 r;
    return poldiv(f,g,q,r);
}

bool poldiv(const Polynomial3& f, const Polynomial3& g)
{
    Polynomial3 r,q;
    return poldiv(f,g,q,r);
}

void SylvesterMatrixElement(const Polynomial3& f,
                            const Polynomial3& g,
                            const int varDegree_f,
                            const int varDegree_g,
                            const int variable,
                            const int i, const int j,
                            int& ex, int& ey, int& ez,
                            double& coeff)
{
    
}

