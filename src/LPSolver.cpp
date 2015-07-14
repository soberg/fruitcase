// http://www.stanford.edu/~liszt90/acm/notebook.html#file16

// Two-phase simplex algorithm for solving linear programs of the form
//
//     maximize     c^T x
//     subject to   Ax <= b
//                  x >= 0
//
// INPUT: A -- an m x n matrix
//        b -- an m-dimensional vector
//        c -- an n-dimensional vector
//        x -- a vector where the optimal solution will be stored
//
// OUTPUT: value of the optimal solution (infinity if unbounded
//         above, nan if infeasible)
//
// To use this code, create an LPSolver object with A, b, and c as
// arguments.  Then, call Solve(x).
#include <cmath>
#include <limits>
#include "LPSolver.h"

LPSolver::LPSolver(const RealMatrix& _A, const RealVector& _b, const RealVector& _c)
{
    N = B = 0;
    
    eps = 10E-07;
    
    m = _b.size();
    n = _c.size();
    D = RealMatrix(m+2, n+2);

    B = new int[m];
    N = new int[n+1];
    
    // copy matrix
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            D(i,j) = _A(i,j);
    for(int i=0; i<m; ++i)
    {
        B[i] = n+i;
        D(i,n) = -1;
        D(i,n+1) = _b(i);
    }
    for(int j=0; j<n; ++j)
    {
        N[j] = j;
        D(m,j) = -_c(j);
    }
    N[n] = -1;
    D(m+1, n) = 1;
}

LPSolver::~LPSolver()
{
    if(N)
    {
        delete[] N;
        N = 0;
    }
    if(B)
    {
        delete[] B;
        B = 0;
    }
}

/*LPSolver::LPSolver(const VVD &A, const VD &b, const VD &c) :
    m(b.size()), n(c.size()), N(n+1), B(m), D(m+2, VD(n+2))
{
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) D[i][j] = A[i][j];
    for (int i = 0; i < m; i++) { B[i] = n+i; D[i][n] = -1; D[i][n+1] = b[i]; }
    for (int j = 0; j < n; j++) { N[j] = j; D[m][j] = -c[j]; }
    N[n] = -1; D[m+1][n] = 1;
}*/

void LPSolver::Pivot(int r, int s)
{
    for(int i=0; i<m+2; ++i)
        if( i != r)
            for(int j=0; j<n+2; ++j)
                if( j != s )
                    D(i,j) -= D(r,j) * D(i,s) / D(r,s);
    for(int j=0; j<n+2; ++j)
        if( j != s )
            D(r,j) /= D(r,s);
    for(int i=0; i<m+2; ++i)
        if( i != r)
            D(i,s) /= -D(r,s);
    D(r,s) = 1.0 / D(r,s);
    
    // swap
    int temp = B[r];
    B[r] = N[s];
    N[s] = temp;
}
    
/*    void LPSolver::Pivot(int r, int s) {
        for (int i = 0; i < m+2; i++) if (i != r)
            for (int j = 0; j < n+2; j++) if (j != s)
                D[i][j] -= D[r][j] * D[i][s] / D[r][s];
        for (int j = 0; j < n+2; j++) if (j != s) D[r][j] /= D[r][s];
        for (int i = 0; i < m+2; i++) if (i != r) D[i][s] /= -D[r][s];
        D[r][s] = 1.0 / D[r][s];
        swap(B[r], N[s]);
    }
 */

bool LPSolver::Simplex(int phase)
{
    int x = phase == 1 ? m+1 : m;
    
    while (true)
    {
        int s = -1;
        for (int j = 0; j <= n; j++)
        {
            if (phase == 2 && N[j] == -1) continue;
            if (s == -1 || D(x,j) < D(x,s) || D(x,j) == D(x,s) && N[j] < N[s]) s = j;
        }
        if (D(x,s) >= -eps) return true; // fertig
        int r = -1;
        for (int i = 0; i < m; i++)
        {
            if (D(i,s) <= 0) continue;
            if (r == -1 || D(i,n+1) / D(i,s) < D(r,n+1) / D(r,s) ||
                D(i,n+1) / D(i,s) == D(r,n+1) / D(r,s) && B[i] < B[r]) r = i;
        }
        if (r == -1) return false; // alle eintraege in der spalte negativ -> unbeschraenkt!
        Pivot(r, s);
    }

}

   /* bool LPSolver::Simplex(int phase) {
        int x = phase == 1 ? m+1 : m;
        while (true) {
            int s = -1;
            for (int j = 0; j <= n; j++) {
                if (phase == 2 && N[j] == -1) continue;
                if (s == -1 || D[x][j] < D[x][s] || D[x][j] == D[x][s] && N[j] < N[s]) s = j;
            }
            if (D[x][s] >= -EPS) return true; // fertig
            int r = -1;
            for (int i = 0; i < m; i++) {
                if (D[i][s] <= 0) continue;
                if (r == -1 || D[i][n+1] / D[i][s] < D[r][n+1] / D[r][s] ||
                    D[i][n+1] / D[i][s] == D[r][n+1] / D[r][s] && B[i] < B[r]) r = i;
            }
            if (r == -1) return false; // alle eintraege in der spalte negativ -> unbeschraenkt!
            Pivot(r, s);
        }
    }*/
double LPSolver::Solve(RealVector &x)
{
    int r = 0;
    for (int i = 1; i < m; i++) if (D(i,n+1) < D(r,n+1)) r = i;
    // D[r][n+1] smalles element
    
    // if < 0
    
    // phase 1
    if (D(r,n+1) <= -eps) {
        Pivot(r, n);
        if (!Simplex(1) || D(m+1,n+1) < -eps) return -std::numeric_limits<double>::infinity();  // infeasible!
        for (int i = 0; i < m; i++)
            if (B[i] == -1)
            {
                int s = -1;
                for (int j = 0; j <= n; j++)
                    if (s == -1 || D(i,j) < D(i,s) || D(i,j) == D(i,s) && N[j] < N[s]) s = j;
                Pivot(i, s);
            }
    }
    // otherwise phase 2...
    if (!Simplex(2)) return std::numeric_limits<double>::infinity();    // unbounded
    x = RealVector(n);
    for (int i = 0; i < m; i++) if (B[i] < n) x(B[i]) = D(i,n+1);
    return D(m,n+1);
}

 /*   DOUBLE LPSolver::Solve(VD &x)
    {
        int r = 0;
        for (int i = 1; i < m; i++) if (D[i][n+1] < D[r][n+1]) r = i;
        // D[r][n+1] smalles element
        
        // if < 0
        
        // phase 1
        if (D[r][n+1] <= -EPS) {
            Pivot(r, n);
            if (!Simplex(1) || D[m+1][n+1] < -EPS) return -numeric_limits<DOUBLE>::infinity();
            for (int i = 0; i < m; i++)
                if (B[i] == -1)
                {
                    int s = -1;
                    for (int j = 0; j <= n; j++)
                        if (s == -1 || D[i][j] < D[i][s] || D[i][j] == D[i][s] && N[j] < N[s]) s = j;
                    Pivot(i, s);
                }
        }
        // otherwise phase 2...
        if (!Simplex(2)) return numeric_limits<DOUBLE>::infinity();
        x = VD(n);
        for (int i = 0; i < m; i++) if (B[i] < n) x[B[i]] = D[i][n+1];
        return D[m][n+1];
    }*/

/*/
 //  main.cpp
 //  testProject
 //
 //  Created by Sören Lennart Berg on 10/1/13.
 //  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
 //
 
 #include <iostream>
 #include <stdio.h>
 #include <getopt.h>
 #include <string.h>
 #include <stdlib.h>
 
 int main(int argc, const char * argv[])
 {
 for(int i=1; i<=72; ++i)
 std::cout << "polytopes/JohnsonPolytopes/J" << i << ".txt ";
 std::cout << std::endl;
 
 double f;
 char aBuffer[256];
 strcpy(aBuffer, "1.5345D0");
 sscanf(aBuffer, "%lf", &f);
 std::cout << f << std::endl;
 std::cout << aBuffer << std::endl;
 aBuffer[ strlen(aBuffer) - 2] = '\0';
 std::cout << aBuffer << std::endl;
 
 return 0;
 
 return 0;
 }*/

