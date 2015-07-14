//
//  UnivariateRoots.cpp
//  densestlatticepackings
//
//  Created by Sören Lennart Berg on 10/03/14.
//  Copyright (c) 2014 Sören Lennart Berg. All rights reserved.
//

#include "UnivariateRoots.h"

UnivariateRoots::UnivariateRoots()
{
    // intialize PARI
    char ch_x[] = "x", ch_y[] = "y", ch_z[] = "z";
	userVariables = cgetg(4, t_VECSMALL);
	userVariables[1] = fetch_user_var(ch_x);
	userVariables[2] = fetch_user_var(ch_y);
	userVariables[3] = fetch_user_var(ch_z);
    
    monomials = cgetg(4, t_VECSMALL);
    for(int i=0; i<3; ++i)
    {
        gel(monomials,i+1) = pol_x_powers(DEFAULT_POLYNOMIAL_LENGTH, userVariables[i+1]);
    }
}

void UnivariateRoots::Roots(const Polynomial3& f, Monomial3 var, LinearSolution3& solutions)
{
    pari_sp av = avma;
    GEN polynomial = getLinearPARIPolynomial(f, var);
    pari_printf("pol be: %Ps \n",polynomial);
    GEN pari_roots = roots(polynomial, DEFAULT_PARI_PRECISION);
    pari_printf("roots be: %Ps \n",pari_roots);
    
    for(int i=1; i < lg(pari_roots); ++i)
    {
        // if root is non-real -> continue
        if( !iszeroeps(PARIscalarToDouble(gmael(pari_roots,i,2)), f.eps))
            continue;
        
        solutions.addSolution(var, PARIscalarToDouble(gmael(pari_roots,i,1)));
    }
    
    avma = av; // clean up PARI stack
}


GEN UnivariateRoots::getPARIPolynomial(const Polynomial3& p)
{
    pari_sp ltop = avma;
    GEN polynomial = mkpoln(1,gen_0); // zero polynomial
    
    for(int i=0; i<DEFAULT_POLYNOMIAL_LENGTH; ++i)
    {
        for(int j=0; j<DEFAULT_POLYNOMIAL_LENGTH; ++j)
        {
            for(int k=0; k<DEFAULT_POLYNOMIAL_LENGTH; ++k)
            {
                if( !iszeroeps(p(i,j,k), p.eps))
                polynomial = gadd(polynomial, gmul( gmael(monomials, 1,i+1),
                                                   gmul(gmael(monomials, 2,j+1),
                                                        gmul(gmael(monomials, 3,k+1),
                                                             dbltor(p(i,j,k)) ))));
            }
        }
    }
    
    return gerepileupto(ltop, polynomial);
}

GEN UnivariateRoots::getLinearPARIPolynomial(const Polynomial3& p, Monomial3 var)
{
    assert( var != VAR_CONST);
    pari_sp ltop = avma;
    GEN polynomial = mkpoln(1,gen_0); // zero polynomial
    
    setvarn(polynomial, userVariables[var+1]);
    
    for(int i=0; i < DEFAULT_POLYNOMIAL_LENGTH; ++i)
    {
        if( !iszeroeps(p.var_permute(i,var,0,0), p.eps) )
            polynomial = gadd(polynomial,
                          gmul( gmael(monomials,var+1,i+1),
                            dbltor( p.var_permute(i,var,0,0) )) );
    }    
    
    return gerepileupto(ltop, polynomial);
}

void UnivariateRoots::getPolynomial3fromPARI(Polynomial3& result, GEN polynomial)
{
    int v[3] = {0,0,0};
    getPolynomial3fromPARI(result, polynomial, v);
}

void UnivariateRoots::getPolynomial3fromPARI(Polynomial3& result, GEN polynomial, int v[3])
{
    if( is_scalar_t(typ(polynomial)) )
    {
        result(0,0,0) = gtodouble(polynomial);
        return;
    }
    
    assert( typ(polynomial) == t_POL );
    
    
    long variables[3];
    
    variables[0] = varn(polynomial);
    
    int w[] = {v[0], v[1], v[2] };
    for(int i=2; i < lg(polynomial); ++i)
    {
        if( t_POL == typ( gel(polynomial, i) ) ) // typ(polynomial[i]) )
        {
            
            getPolynomial3fromPARI(result, gel(polynomial,i), w);
        }
        else if( is_scalar_t(typ(gel(polynomial,i))) )
        {
            result(w[0],w[1],w[2]) = gtodouble( gel(polynomial,i) );
        }
        w[ variables[0] ]++;
    }
}

double UnivariateRoots::PARIscalarToDouble(GEN g)
{
    assert(is_scalar_t(typ(g)));
    long type = typ(g);
    
    switch(type)
    {
        case t_INT:
            return (double)itos(g);
            break;
        case t_FRAC:
            return ((double)itos(gel(g,1)))/((double)itos(gel(g,2)));
            break;
        case t_REAL:
            return rtodbl(g);
            break;
    }
    
    // TODO error handling?
    return 0.0;
}