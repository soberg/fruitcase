//
//  linear_programming.cc
//  Linear Programming
//
//  Created by Sören Lennart Berg on 8/17/12.
//  Copyright (c) 2012 Sören Lennart Berg. All rights reserved.
//

/*! \file linear_programming.cpp
    \brief Simplex Method

    This file provides an implementation for a basic two-phase simplex method to solve linear programs.
    It also includes a simple implementation of the gaussion elimination method using pivots for solving
    system of linear equations.
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "linear_programming.h"

/*!
    Print a vector, i.e. a one dimensional double array, of length n
    @param vec The vector to be printed
    @param n   Length of vec
*/
void print_vector(double const* vec, const int n)
{
    assert(n>0);
    printf("\n");
    for(int i=0; i<n; ++i)
        printf(" %lf",vec[i]);
    printf("\n");
}

/*!
    Print a basis, i.e. a one dimensional integer array, for the simplex method.
    @param vec The vector/basis to be printed
    @param n   Length of vec/basis
*/
void print_basis(int const* vec, const int n)
{
    assert(n>0);
    printf("\n");
    for(int i=0; i<n; ++i)
        printf(" %i",vec[i]);
    printf("\n");
}

/*!
    Output a matrix, i.e. a two dimensional array
    @param Two dimensional array to be printed
    @param rows Number of rows
    @param cols Number of columns
*/
void print_matrix(double** mat, const int rows, const int cols)
{
    int i,j;

    assert( rows>0 && cols>0);

    printf("\n");
    for(i=0;i<rows;++i)
    {
        for(j=0;j<cols;++j)
            printf(" %lf", mat[i][j]);
        printf("\n");
    }
}

// returns index >=k of greatest absolute value in column j of mat
/*!
    Finds the index of the element with greates absolute value of a specific column of matrix
    @param mat Matrix
    @param rows Number of rows of mat
    @param cols Number of columns of mat
    @param j Column index of mat
    @param k Starting index for column j
    @return The index of column j having the greatest absolute value among all elements in column j with index at least k
*/
int argmax_abs(double * const * const mat, const int rows, const int cols, const int j, const int k)
{
    int i, ind_max;
    double max;

    assert( rows>0 && cols>0 && j<cols && k<rows );

    ind_max = k;
    max = fabs(mat[k][j]);

    for(i=k+1; i<rows; ++i)
        if(fabs(mat[i][j]) > max)
        {
            max = fabs(mat[i][j]);
            ind_max = i;
        }
    return ind_max;
}

// swaps i-th and k-th row
/*!
    Swap two rows of a matrix
    @param mat Matrix. Will be modified by this function.
    @param rows Number of rows of mat
    @param cols Number of columns of mat
    @param i Row index
    @param k Row index
*/
void swap_rows(double** mat, const int rows, const int cols, const int i, const int k)
{
    double* temp;

    assert( rows>0 && cols>0 && i<rows && k<rows );

    temp = mat[i];
    mat[i] = mat[k];
    mat[k] = temp;
}

/*!
    Checks if a double should be treated as zero
    @param eps_ Greatest non-negative value, which will be treated as zero
    @return 0 if d is not zero, 1 otherwise
*/
int dzero_eps(const double d, const double eps_)
{
    return fabs(d) <= eps_;
}

// solve system of linear equations
// the matrix mat will be changed. a inhomogenous solution will be written to sol and
// a basis of kernel(mat) will be written to ker, rnk will contain rank(mat) after computations
// allocates memory for sol {cols-1 doubles} and ker {(cols-1-(*rnk))*(cols-1) doubles)}
void solve_linearsystem(double** mat, const int rows, const int cols, const double eps_, double** sol, double*** ker, int* rnk)
{
    int i,j,k,l;
    int piv;
    int lde;
    int var;
    int ker_var;
    double frac;

    (*sol) = (double*)malloc((cols-1)*sizeof(double));
    for(j=0;j<cols-1;++j)
        (*sol)[j] = 0.0;

    assert(mat && cols>1 && rows>0);

    (*rnk) =0;

    // forward elimination
    for(j=0; j<cols-1 && (*rnk)<rows; ++j)
    {
        //print_matrix(mat,3,3);

        // determine pivot element
        piv = argmax_abs(mat, rows, cols, j, (*rnk));
        //printf("j=%i piv=%i rnk=%i pivval=%lf\n",j,piv,*rnk,mat[piv][j]);
        if(dzero_eps(mat[piv][j],eps_)) continue;
        swap_rows(mat, rows, cols, piv, (*rnk));

        for(k=(*rnk)+1; k<rows; ++k)
        {
            frac = mat[k][j] / mat[(*rnk)][j];

            mat[k][j] = 0.0;

            for(l=j+1;l<cols;++l)
                mat[k][l] -= frac*mat[(*rnk)][l];
        }
        (*rnk)++;
        //print_matrix(mat,rows,cols);

    }


    // initialize kernel vectors
    (*ker) = (double**)malloc( (cols-1-(*rnk))*sizeof(double*));
    for(k=0; k<(cols-1-(*rnk)); ++k)
        (*ker)[k] = (double*)malloc( (cols-1)*sizeof(double));
    ker_var = cols-1-(*rnk)-1;
    for(k=0; k < cols-1-(*rnk); ++k)
    {
        for(j=0;j<cols-1;++j)
            (*ker)[k][j] = 0.0;
    }


    // backward elimination
    var = cols-1;
    for(i=rows-1; i>=0; --i)
    {
        //print_matrix(mat,rows,cols);
        piv = 0;    // determine first non zero entry
        while( dzero_eps(mat[i][piv],eps_) && piv<cols)
            piv++;

        if( piv == cols) continue;  // zero row
        if( piv == cols-1) return;  // unsolvable todo here...

        frac = 1.0 / mat[i][piv];
        for(k=i-1; k>=0; --k)
        {
            for(j=piv+1; j<cols;++j)
                mat[k][j] -= mat[k][piv]*frac*mat[i][j];
            mat[k][piv] = 0.0;
        }

        (*sol)[piv] = mat[i][cols-1]*frac;

        //print_matrix(mat,rows,cols);

        var = piv;
    }

    //print_matrix(mat,rows,cols);

    // determine kernel
    var = cols-1;
    for(i=rows-1; i>=0; --i)
    {
        piv = 0;    // determine first non zero entry
        while( dzero_eps(mat[i][piv],eps_) && piv<cols)
            piv++;

        if( piv == cols) continue;  // zero row

        if(var - piv > 1)       // are there free variables?
        {
            for(k=0; k<var-piv-1; ++k)
               (*ker)[ker_var-k][var-k-1] = 1.0;

            for(l=i; l>=0; --l)
            {
                // determine first non zero entry in row l
                lde = 0;
                while( dzero_eps(mat[l][lde],eps_) && lde<cols)
                    lde++;

                frac = (-1.0)/mat[l][lde];

                for(k=0; k<var-piv-1; ++k)
                    (*ker)[ker_var-k][lde] = mat[l][var-k-1]*frac;

            }
            ker_var -= var-piv-1;
        }
        var = piv;
        //ker_var -= var-piv-1;
    }

    //print_matrix(mat,rows,cols);

    //print_matrix(ker, cols-1-rnk, cols-1);

    //print_vector(sol, cols-1);
}

// find a pivot element (j,t) according to Bland's rule
// if j==-1 the tableau is optimal and unbounded for t==-1
void simplexPivotBlandsRule(double** Tab, const int rows, const int cols, const double eps_, int* j, int *t)
{
    int i;
    double min;
    double frac;
    
    // choose smallest (column) index j, with negative costs
    *j=0;
    //while( Tab[0][*j]+eps_ >= 0.0 && *j<cols)
    while( Tab[0][*j] >= 0.0 && *j<cols)
        (*j)++;
    
    if( *j >= cols-1)
    {   *j= -1; // optimal solution found!
        return;
    }
    
    // choose smallest (row) index t, satisfying feasibility
    *t=-1;
    min = -1.0;
    
    for(i=1; i<rows; ++i)
        if( Tab[i][*j]-eps_ > 0.0)
        {
            frac = Tab[i][cols-1] / Tab[i][*j];
            if( frac+eps_ < min || min < 0.0)
            {
                *t = i;
                min = frac;
            }
        }
    
    // if we haven't found any t => unbounded!
    if( -1 == *t)
        return;
    
    return;
}

// find a pivot element (j,t) according to ...
// if j==-1 the tableau is optimal and unbounded for t==-1
void simplexPivotDantzigsRule(double** Tab, const int rows, const int cols, const double eps_, int* j, int *t)
{
    int i;
    double min;
    double frac;
    
    // choose smallest (column) index j, with negative costs
    *j=0;
    //while( Tab[0][*j]+eps_ >= 0.0 && *j<cols)
    while( Tab[0][*j] >= 0.0 && *j<cols)
        (*j)++;
    
    if( *j >= cols-1)
    {   *j= -1; // optimal solution found!
        return;
    }
    
    // choose smallest (row) index t, satisfying feasibility
    *t=-1;
    min = -1.0;
    
    for(i=1; i<rows; ++i)
        //if( Tab[i][*j]-eps_ > 0.0)
        if( Tab[i][*j] > eps_)//0.0)
        {
            frac = Tab[i][cols-1] / Tab[i][*j];
            //if( frac+eps_ < min || min < 0.0)
            if( frac < min || min < 0.0)
            {
                *t = i;
                min = frac;
            }
        }
    
    // if we haven't found any t => unbounded!
    if( -1 == *t)
        return;
    
    //int j;
    double minCost = 0.0;
    bool multChoice = false;
    
    *j = -1;
    for( i = 0; i < cols; i++)
    {
        if( Tab[0][i] >= eps_)//0.0)
            continue;
        
        if( Tab[0][i] <= minCost + eps_ && minCost <= Tab[0][i] + eps_)
            multChoice = true;
        
        if( Tab[0][i] < minCost)
        {
            minCost = Tab[0][i];
            *j = i;
            multChoice = false;
        }
    }
    
    if(multChoice) // use blands rule as fallback to avoid cycles
    {
        simplexPivotBlandsRule(Tab, rows, cols, eps_, j, t);
        return;
    }
    
    if( -1 == *j)   // optimal solution found;
        return;
    
    // otherwise find fitting *t
        
    // choose smallest (row) index t, satisfying feasibility
    *t=-1;
    min = -1.0;
    
    for(i=1; i<rows; ++i)
        //if( Tab[i][*j]-eps_ > 0.0)
        if( Tab[i][*j] > eps_)//0.0)
        {
            frac = Tab[i][cols-1] / Tab[i][*j];
            //if( frac+eps_ < min || min < 0.0)
            if( frac < min || min < 0.0)
            {
                *t = i;
                min = frac;
            }
        }
    
    // if we haven't found any t => unbounded!
    if( -1 == *t)
        return;
    
    return;
}

// return 1 = ok, 0 = unbounded
int tableau_pivoting(double** Tab, const int rows, const int cols, int* bas, const double eps_)
{
    double frac;
    int i,j,t,k;

    //print_matrix(Tab, rows, cols);
    //printf("in\n");
    while(1)
    {
        // find pivot (t,j)
        simplexPivotBlandsRule(Tab, rows, cols, eps_, &j, &t);
        //simplexPivotDantzigsRule(Tab, rows, cols, eps_, &j, &t);

        if( -1 == j)
        {
			//printf("boya\n");
			//break; // optimal
            return 1;
		}

        if( -1 == t)
            return 0;
            //break; // unbounded todo do something

        // update basic
        bas[t-1] = j;

        assert(Tab[t][j] != 0.0);
        frac = (1.0)/Tab[t][j];   
        for(k=0;k<cols;++k)
        {   Tab[t][k] *= frac;
            if( Tab[t][k] > -eps_ && Tab[t][k] < eps_)
                Tab[t][k] = 0.0;
//            if( fabs(Tab[t][k]) < eps_) // numerical stability
//                Tab[t][k] = 0.0;
        }

        //print_matrix(Tab, rows, cols);

        Tab[t][j] = 1.0; // numerical stability

        for(i=0; i<rows; ++i)
        {
            if( i==t)
                continue;

            frac = Tab[i][j] / Tab[t][j];

            for(k=0;k<cols;++k)
            {   Tab[i][k] -= Tab[t][k]*frac;
//                if( fabs(Tab[i][k]) < eps_) // numerical stability
//                    Tab[i][k] = 0.0;
                if( Tab[i][k] > -eps_ && Tab[i][k] < eps_)
                    Tab[i][k] = 0.0;

            }

            Tab[i][j] = 0.0; // numerical stability
        }

        //print_matrix(Tab, rows, cols);
    }
    //printf("out\n");
    return 1;
}

// solve LP: max c*x st Ax<=b, x>=0, where b>=0
// allocates memory for x{n doubles}
void simplex_nonneg_rhs(double** A, const int m, const int n, double const* b, double const* c, double** x, const double eps_)
{
    int i,j;
    double** Tab;
    int* bas;

    assert(m > 0 && n >0 && A && b && c);

    // STEP 1:  build tableau
    Tab = (double**)malloc( (m+1)*sizeof(double*));

    // fill first row of tableau
    Tab[0] = (double*)malloc( (m+n+1)*sizeof(double));
    for(i=0;i<n;++i)
        Tab[0][i] = (-1.0)*c[i];
    for(i=n;i<m+n+1;++i)
        Tab[0][i] = 0.0;

    for(i=1; i<m+1; ++i)
    {
        Tab[i] = (double*)malloc( (m+n+1)*sizeof(double));
        for(j=0;j<n;++j)
            Tab[i][j] = A[i-1][j];
        for(j=n;j<m+n;++j)
            Tab[i][j] = 0.0;
        Tab[i][n+i-1] = 1.0;
        Tab[i][n+m] = b[i-1];
    }
    bas = (int*)malloc( m*sizeof(int));     // initialize basis
    for(i=0; i<m; ++i)
        bas[i] = n + i;


    //print_matrix(Tab, m+1, n+m+1);
    //print_basis(bas, m);

    // STEP 2: pivoting steps
    tableau_pivoting(Tab, m+1, m+n+1, bas, eps_);

    //print_matrix(Tab, m+1, n+m+1);
    //print_basis(bas, m);

    // extract solution
    *x = (double*)malloc(n*sizeof(double));
    for(j=0; j<m; ++j) (*x)[j] = 0.0;
    for(i=0;i<m;i++)
        if( bas[i]<n)
            (*x)[bas[i]] = Tab[i+1][m+n];

    // free memory
    for(i=0;i<m+1;++i)
		free(Tab[i]);
	free(Tab);
	free(bas);
}

// solve LP: max c*x st Ax<=b, x>=0, b \in R^n
// if flag_nonopt is set to nonzero, just a feasible point is returned, i.e. no phase II is performed
// allocates memory for x if there's an solution{n doubles}
void solveLinearProgramNonnegative(double** A, const int m, const int n, double const* b, double const* c, double** x, const int flag_nonopt, const double eps_)
{
    double** Tab, **Tab2;   // tableaus for phase I/II
    int i,j,k,l;
    double d;
    int num_neg;            // # negative rhs entries
    int* bas;               // current basis
    double sum_b2;
    //double* x;              // optimal basic solution

    // phase I: try to find feasible vertex

    //printf("solve lp flag: %i\n",flag_nonopt);
   // print_matrix(A,m,n);

    // build tableau
    num_neg = 0;        // determine number of neg. RHS entries
    for(i=0; i<m; ++i)
        if( b[i] + eps_ < 0.0)
            num_neg++;

    Tab = (double**)malloc( (m+1) * sizeof(double*));

    Tab[0] = (double*)malloc( (m+n+num_neg+1)*sizeof(double));
    for(j=0; j<m+n+num_neg+1;++j)
        Tab[0][j] = 0.0;

    k=1;
    l=m-num_neg+1;
    sum_b2=0.0;
    //printf("solve lp1\n");
    //print_matrix(Tab,m+1,m+n+num_neg+1);
    for(i=0; i<m; ++i)
    {
        //printf("solve lp2\n");
        if( b[i] >= -eps_) //!(b[i] < 0.0) )
        {
            Tab[k] = (double*)malloc( (m+n+num_neg+1)*sizeof(double));
            for(j=0; j<n; ++j)
                Tab[k][j] = A[i][j];
            for(j=n; j<m+n+num_neg; ++j)
                Tab[k][j] = 0.0;
            Tab[k][n+num_neg+k-1] = 1.0;//
            Tab[k][m+n+num_neg] = b[i];
            k++;
        }
        else
        {
            Tab[l] = (double*)malloc( (m+n+num_neg+1)*sizeof(double));
            for(j=0; j<n; ++j)
            {
                Tab[l][j] = (-1.0)*A[i][j];
                Tab[0][j] += A[i][j];
            }
            for(j=n; j<m+n+num_neg; ++j)
                Tab[l][j] = 0.0;
            
            Tab[l][n+l-(m-num_neg+1)] = -1.0;
            Tab[l][n+num_neg+l-1] = 1.0;
           // Tab[0][n+l-(m-num_neg+1)] = 1.0;  ////
            Tab[0][n+l-(m-num_neg+1)] = 1.0;
            Tab[l][m+n+num_neg] = (-1.0)*b[i];
            sum_b2 += (-1.0)*b[i];
            l++;
        }
        // identity matrix for slack vars
        //Tab[i+1][n+num_neg + i] = 1.0;
    }
    
    bas = (int*)malloc( m*sizeof(int));
    for(i=0; i<m; ++i)
        bas[i] = n+num_neg+i;

   //     printf("solve lp333\n");
   // print_matrix(Tab,m+1,m+n+num_neg+1);

    //print_matrix(Tab,m+1,m+n+num_neg+1);
    
    //printf
//    printDoubleMatrix(Tab, m+1, m+n+num_neg+1);

    int isbounded = tableau_pivoting(Tab, m+1, m+n+num_neg+1, bas, eps_);
    ///printf("solve lp4\n");
    //print_matrix(Tab,m+1,m+n+num_neg+1);

    //print_basis(bas, m);

    // test if feasible or infeasible
    
    if( 1==isbounded && Tab[0][m+n+num_neg] > sum_b2 + eps_)
    {    printf(" error solveLinearSystemNonnegative this should never happen \n");
        printf(" error size: %lf \n", Tab[0][m+n+num_neg] - sum_b2);
        print_matrix(A,m,n);
    }
    
    if( !flag_nonopt && (Tab[0][m+n+num_neg] < sum_b2 - eps_ || !isbounded) )  // if we want to optimize and it's unbounded  ... or the problem is infeasible -> done
    {
        // infeasible!
       // printf("infeas!!!\n");
        (*x)=0;
    }
    
    if(flag_nonopt)
    {
       // if(isbounded)
       // printf("isbounded\n");
        if( Tab[0][m+n+num_neg] < sum_b2 - eps_)
        {
            (*x)=0;
          //  printf("ua\n");
        }
        else
        {
        // write a feasible solution to x
        (*x) = (double*)malloc( n*sizeof(double));

        for(j=0;j<n;++j)
			(*x)[j] = 0.0;

        for(i=0; i<m; ++i)
            if(bas[i] < n)
            {    (*x)[ bas[i] ] = Tab[i+1][m+n+num_neg];
				//printf("i: %i\n",i);
			}
        }
    }
    else if(!flag_nonopt)
    {
        (*x) = (double*)malloc( n*sizeof(double));
        for(j=0;j<n;++j)
            (*x)[j] = 0.0;

        // tableau is feasible -> calculate a feasible vertice to start with

        // if 'negative-slack' variables are still in the basis, fix them
        for(i=0; i<m; ++i)
            if( bas[i] >= m+n && bas[i] < m+n+num_neg)
            {   for(j=0; j<m+n+num_neg+1; ++j)
                    Tab[i+1][j] *= (-1.0);
                bas[i] = bas[i] - m;
            }

        // now build the new tableau, for phase II
        Tab2 = (double**)malloc( (m+1)*sizeof(double*));
        for(i=0; i<m+1; ++i)
            Tab2[i] = (double*)malloc( (m+n+1)*sizeof(double));

        for(j=0; j<n; ++j)
            Tab2[0][j] = (-1.0)*c[j];
        for(j=n; j<m+n+1;++j)
            Tab2[0][j] = 0.0;

        for(i=1; i<m+1; ++i)
        {
            for(j=0; j<m+n; ++j)
                Tab2[i][j] = Tab[i][j];
            Tab2[i][m+n] = Tab[i][m+n+num_neg];
        }
        //print_matrix(Tab2, m+1, m+n+1);
        //print_basis(bas, m);
        // fix basic columns
        for(i=0;i<m;++i)
        {
            d = Tab2[0][bas[i]];
            //if(bas[i]>=n)
            //    continue;
            //frac = 1.0/Tab2[0][bas[i]];
            for(j=0;j<m+n+1;++j)
                Tab2[0][j] -= Tab2[i+1][j] * d;
            //Tab2[0][bas[i]] = 0.0; // num. stability

            //print_matrix(Tab2, m+1, m+n+1);
        }

        //print_matrix(Tab2, m+1, m+n+1);

        // phase II
        tableau_pivoting(Tab2, m+1, m+n+1, bas, eps_);

        // write solution to x
        for(i=0; i<m; ++i)
            if(bas[i] < n)
                (*x)[ bas[i] ] = Tab2[i+1][m+n];
    }



    //print_vector(*x, n);

    // free memory
    for(i=0;i<m+1;++i)
        free(Tab[i]);
    if(Tab)
        free(Tab);
    if(bas)
        free(bas);

    if(!flag_nonopt)
    {   for(i=0;i<m+1;++i)
            free(Tab2[i]);
        free(Tab2);
    }
}

// writes a feasible point for A*x<=b, x>=0 to x if there is some, otherwise set x to NULL
/*!
 Finds a feasible point for */


void feasiblePointNonnegative(double** A, const int m, const int n, double const* b, double** x, const double eps_)
{
    solveLinearProgramNonnegative(A, m, n, b, NULL, x, 1, eps_);
    
    // for testing purposes we just soplex
    // using namespace soplex;
    /*using namespace soplex;
    
    SoPlex spl;
    spl.changeSense(SPxLP::MINIMIZE);
    
    int* zero_n = (int*)malloc( n * sizeof(int));
    for(int i=0; i<n; ++i)
        zero_n[i] = i;
    
    for(int i=0; i<m; ++i)
    {
        DSVector vec(n);
        vec.add(n, zero_n, A[i]);
        spl.addRow(LPRow( -infinity, vec, b[i]));
    }
    SPxSolver::Status st = spl.solve();
    
    if(st != SPxSolver::OPTIMAL && st != SPxSolver::UNBOUNDED && st != SPxSolver::INFEASIBLE)
        std::cout << "problem with soplex" << std::endl;*/
}

//bool feasiblePointNonnegativeSP(double** A, const int m, const int n, double const* b, double** x, const double eps_)
//{
//    //solveLinearProgramNonnegative(A, m, n, b, NULL, x, 1, eps_);
//    
//    // for testing purposes we just soplex
//    // using namespace soplex;
//    using namespace soplex;
//    
//    SoPlex spl;
//    spl.changeSense(SPxLP::MINIMIZE);
//   // spl.setDelta(eps_);
//    
//    soplex::Param::setEpsilon(eps_);
//    
//    int* zero_n = (int*)malloc( n * sizeof(int));
//    for(int i=0; i<n; ++i)
//        zero_n[i] = i;
//    
//    for(int i=0; i<m; ++i)
//    {
//        DSVector vec(n);
//        //vec.add(n, zero_n, A[i]);
//        for(int j=0; j<n; ++j) vec.add(j, A[i][j]);
//        spl.addRow(LPRow( -infinity, vec, b[i]));
//                
//    }
//    SPxSolver::Status st = spl.solve();
//    
//    if(st != SPxSolver::OPTIMAL && st != SPxSolver::UNBOUNDED && st != SPxSolver::INFEASIBLE)
//        std::cout << "problem with soplex" << std::endl;
//    
//    if(SPxSolver::INFEASIBLE == st)
//        return false;
//    return true;
//}

/*/ if the i-th position in nonneg is not zero, x_i will be nonnegative, i.e. x_i>=0
void solve_lp(double** A, int m, int n, double* b, int* nonneg, double* c, double** x, int flag_nonopt, double eps_)
{
    int i,k, num_uvar, offset;
    double** A_;

    num_uvar = 0;    // #signed variables

    for(k=0;k<n;++k)
        if(0 == nonneg[k]) ++num_uvar;

    A_ = (double**)malloc( (n+num_uvar) * sizeof(double*));

    offset = 0;
    for(k=0; k<n; ++k)
    {
        if(nonneg[k] != 0)
            continue;

        A_[k+offset] = A[k];
        for(i=0;i<m;++i)
            A_[k+offset+1] = (-1.0)*A[k][i]
    }
}*/


/*/ solve lp in general form Ax<=b x>=0 by two phase Simplex
void solve_lp(double** mat, int m, int n, double* b, double* c, double eps_)
{
    int i,j,k,l;
    //double

    assert(m > 0 && n >0 && mat && b && c);

    // build tableau

}*/



