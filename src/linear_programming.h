//
//  linear_programming.h
//  Linear Programming
//
//  Created by Sören Lennart Berg on 8/17/12.
//  Copyright (c) 2012 Sören Lennart Berg. All rights reserved.
//

#ifndef _DENSEST_LATTICE_PACKINGS_LINEAR_PROGRAMMING_HEADER_
#define _DENSEST_LATTICE_PACKINGS_LINEAR_PROGRAMMING_HEADER_

//#include "soplex.h"
#include "helper.h"
//#include <spxdefines.h>

void print_vector(double const* vec, const int n);
void print_basis(int const* vec, const int n);
void print_matrix(double ** mat, const int rows, const int cols);

void solveLinearProgramNonnegative(double** A, const int m, const int n, double const* b, double const* c, double** x, const int flag_nonopt, const double eps_);
void feasiblePointNonnegative(double** A, const int m, const int n, double const* b, double** x, const double eps_);
void solve_linearsystem(double** mat, const int rows, const int cols, const double eps_, double** sol, double*** ker, int* rnk);

bool feasiblePointNonnegativeSP(double** A, const int m, const int n, double const* b, double** x, const double eps_); // TODO tempori

/*!
 This function writes an inequality to the system of inequalities represented by the matrix matLHS and the right hand
 side vector vecRHS, i.e. \f$ matLHS \cdot x \leq vecRHS \f$
 @param matLHS left hand side matrix of the system
 @param vecRHS right hand side vector of the system
 @param index Index of inequality / row number of matLHS and vecRHS to write to
 @param numVariables Number of variables / number of columns of matLHS
 @param lhs vector describing the lhs of the inequality
 @param rhs right hand side of the inequality
 */
inline void writeInequality(double** matLHS, double* vecRHS, const int index, const int numVariables, double const* lhs, const double rhs)
{
    int j;
    
    for(j=0; j<numVariables; ++j)
        matLHS[index][j] = lhs[j];
    vecRHS[index] = rhs;
}

inline void writeEquality(double** matLHS, double* vecRHS, const int indexPlus, const int indexMinus, const int numVariables, double const* lhs, const double rhs)
{
    int j;
    
    for(j=0; j<numVariables; ++j)
    {
        matLHS[indexPlus][j]  = lhs[j];
        matLHS[indexMinus][j] = lhs[j] * (-1.0);
    }
    vecRHS[indexPlus ] = rhs;
    vecRHS[indexMinus] = rhs * (-1.0);
}
/*!
 This function writes an inequality to the system of inequalities represented by the matrix matLHS and the right hand
 side vector vecRHS, i.e. \f$ matLHS \cdot x \leq vecRHS \f$. This function will create a positive / negative ..... TODO REFACTOR!!
 @param matLHS left hand side matrix of the system
 @param vecRHS right hand side vector of the system
 @param index Index of inequality / row number of matLHS and vecRHS to write to
 @param numVariables Number of variables / number of columns of matLHS
 @param lhs vector describing the lhs of the inequality
 @param rhs right hand side of the inequality
 */
inline void writeInequalityPosNegVariables(double** matLHS, double* vecRHS, const int index, const int numVariables, double const* lhs, const double rhs)
{
    int j;
    
    for(j=0; j<numVariables; ++j)
    {
        matLHS[index][2 * j    ] = lhs[j];
        matLHS[index][2 * j + 1] = lhs[j] * (-1.0);
    }
    vecRHS[index] = rhs;
}

inline void writeEqualityPosNegVariables(double** matLHS, double* vecRHS, const int indexPlus, const int indexMinus, const int numVariables, double const* lhs, const double rhs)
{
    int j;
    
    for(j=0; j<numVariables; ++j)
    {
        matLHS[indexPlus][2 * j    ] = lhs[j];
        matLHS[indexPlus][2 * j + 1] = lhs[j] * (-1.0);
        
        matLHS[indexMinus][2 * j    ] = lhs[j] * (-1.0);
        matLHS[indexMinus][2 * j + 1] = lhs[j];
    }
    vecRHS[indexPlus ] = rhs;
    vecRHS[indexMinus] = rhs * (-1.0);
}

#endif // _DENSEST_LATTICE_PACKINGS_LINEAR_PROGRAMMING_HEADER_

