//
//  minima.c
//  
//
//  Created by Sören Lennart Berg on 9/6/13.
//
//

#include "Minima.h"

Minima::Minima()
{
//    // intialize PARI
//    char ch_x[] = "x", ch_y[] = "y", ch_z[] = "z";
//	userVariables = cgetg(4, t_VECSMALL);
//	userVariables[1] = fetch_user_var(ch_x);
//	userVariables[2] = fetch_user_var(ch_y);
//	userVariables[3] = fetch_user_var(ch_z);
//    
//    monomials = cgetg(4, t_VECSMALL);
//    for(int i=0; i<3; ++i)
//    {
//        gel(monomials,i+1) = pol_x_powers(DEFAULT_POLYNOMIAL_LENGTH, userVariables[i+1]);
//    }
    
    eps = 10E-07;
}

//GEN Minima::getPARIPolynomial(Polynomial3 p)
//{
//    pari_sp ltop = avma;
//    GEN polynomial = gen_0;//dbltor(p(0,0,0));
//    
//    // GEN powgi(GEN x, GEN y)
//    // mkintn
//    
//    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
//    {
//        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
//        {
//            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
//            {
//                if( !iszeroeps(p(i,j,k), p.eps))
//                    polynomial = gadd(polynomial, gmul( gmael(monomials, 1,i+1),
//                                                   gmul(gmael(monomials, 2,j+1),
//                                                        gmul(gmael(monomials, 3,k+1),
//                                                             dbltor(p(i,j,k)) ))));
//            }
//        }
//    }
//    
//    return gerepileupto(ltop, polynomial);
//}
//
//void Minima::getPolynomial3fromPARI(Polynomial3& result, GEN polynomial)
//{
//    int v[3] = {0,0,0};
//    getPolynomial3fromPARI(result, polynomial, v);
//}
//
//void Minima::getPolynomial3fromPARI(Polynomial3& result, GEN polynomial, int v[3])
//{
//    if( is_scalar_t(typ(polynomial)) )
//    {
//        result(0,0,0) = gtodouble(polynomial);
//        return;
//    }
//    
//    assert( typ(polynomial) == t_POL );
//    
//    
//    long variables[3];
//    
//    variables[0] = varn(polynomial);
//    
//    int w[] = {v[0], v[1], v[2] };
//    for(int i=2; i < lg(polynomial); ++i)
//    {
//        if( t_POL == typ(polynomial[i]) )
//        {
//            
//            getPolynomial3fromPARI(result, gel(polynomial,i), w);
//        }
//        else if( is_scalar_t(typ(gel(polynomial,i))) )
//        {
//            result(w[0],w[1],w[2]) = gtodouble( gel(polynomial,i) );
//        }
//        w[ variables[0] ]++;
//    }
//}

// sort Polynomials in polynomials[3] with respect to the order in degrees (e.g. (0,1,2) or (2,1,0) etc.) and writes to ordered[]
void sortPolynomials(const Polynomial3 polynomials[3], Polynomial3 ordered[3], int degrees[3])
{
    // we use the inelegant way of if statements here in order to avoid unnecessary copying of polynomials
    if( degrees[0] <= degrees[1] && degrees[0] <= degrees[2])
    {
        ordered[0] = polynomials[0];
        if( degrees[1] <= degrees[2])
        {
            ordered[1] = polynomials[1];    // [0] <= [1] <= [2]
            ordered[2] = polynomials[2];
        }
        else
        {
            ordered[1] = polynomials[2];    // [0] <= [2] <= [1]
            ordered[2] = polynomials[1];
            swapints(&degrees[1], &degrees[2]);
        }
    }
    else if( degrees[1] <= degrees[0] && degrees[1] <= degrees[2])
    {
        ordered[0] = polynomials[1];
        if( degrees[0] <= degrees[2])
        {
            ordered[1] = polynomials[0];    // [1] <= [0] <= [2]
            ordered[2] = polynomials[2];
            swapints(&degrees[0], &degrees[1]);
        }
        else
        {
            ordered[1] = polynomials[2];    // [1] <= [2] <= [0]
            ordered[2] = polynomials[0];
            swapints(&degrees[0],&degrees[1]);
            swapints(&degrees[1],&degrees[2]);
        }
    }
    else if( degrees[2] <= degrees[0] && degrees[2] <= degrees[1])
    {
        ordered[0] = polynomials[2];
        if( degrees[0] <= degrees[1])
        {
            ordered[1] = polynomials[0];    // [2] <= [0] <= [1]
            ordered[2] = polynomials[1];
            swapints(&degrees[0],&degrees[2]);
            swapints(&degrees[2],&degrees[1]);
            
        }
        else
        {
            ordered[1] = polynomials[1];    // [2] <= [1] <= [0]
            ordered[2] = polynomials[0];
            swapints(&degrees[0],&degrees[2]);
        }
    }
}

void Minima::MinimizeDeterminant(Polynomial3& polynomial)
{
    // todo: necessary to have 3 here?
    Polynomial3 derivatives[3];
    int numPolynomials = 3;
    
    for(int i=0; i<3; ++i)
        polynomial.derive(derivatives[i],POL3_VARS[i]);
    
    FindAffineVarieties(derivatives, 3);
}

void Minima::FindAffineVarieties(Polynomial3* polynomials, int numPolynomials)
{
    Polynomial3 p[3];
    int totaldegrees[3];
    int degrees_x[3];
    
    for(int i=0; i<numPolynomials; ++i)
    {
        totaldegrees[i] = polynomials[i].totalDegree();
        degrees_x[i] = polynomials[i].degree(VAR_X);
        if( totaldegrees[i] > 2)
        {
            std::cout << "error in Minima::Solve: total degree exceeds 2 and is" << totaldegrees[i] << std::endl;
            return;
        }
    }

    // delete zero polynomials
    for(int i=0; i<numPolynomials; ++i)
    {
        if( -1 > totaldegrees[i])
        {
            for(int j=i+1; j<numPolynomials; ++j)
            {
                polynomials[i] = polynomials[i+1];
            }
            numPolynomials--;
        }
    }
    
    // trivial case I: all polynomials vanish:
    if( 0 == numPolynomials )
    {
        // todo all solutions
        return;
    }
    
    // trivial case II: there's a constant, non-zero polynomial
    for(int i=0; i<numPolynomials; ++i)
        if( totaldegrees[i] == 0 )
        {
            // todo fill out: no solution
            return;
        }
    
    // sort polynomials with respect to degree of x (ascending order)
    sortPolynomials(polynomials, p, degrees_x);
    
    // if we have a ocurrence of x^2 in the polynomials, we use a gauss-like elimination
    // such that only the polynomial polynomials[numPolynomials-1] contains a 'x^2 - monomial'.
    // see paper...
    if( 2 == degrees_x[numPolynomials-1] )
    {
        for(int i=0; i<numPolynomials-1; ++i)
        {
            polynomials[i] -= ( p[i](2,0,0)/p[numPolynomials-1](2,0,0) )*p[numPolynomials-1];
        }
    }
    
    
}