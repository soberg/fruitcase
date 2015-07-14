/*! \file helper.cpp
 \brief Various helper functions
 */

#include "helper.h"

/*!
 For time measuring
 */
void timerGetTime(timeval* timeInterval)
{
    gettimeofday(timeInterval, NULL);
}

/*!
 For time measuring
 */
double getTimeMilliseconds(timeval timeStart, timeval timeEnd)
{
    double elapsedTime = 0.0;
    elapsedTime = (timeEnd.tv_sec - timeStart.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (timeEnd.tv_usec - timeStart.tv_usec) / 1000.0;   // us to ms
    return elapsedTime;
}

void freeDoubleMatrix(double** mat, int m)
{
    int i;
    if(!mat)
        return;
    for(i=0;i<m;++i)
        if(mat[i])
        {
            free(mat[i]);
            mat[i] = 0;
        }

    free(mat);
}

double** allocateDoubleMatrix(int n, int m)
{
    int i;
    double** mat = (double**)malloc(n * sizeof(double*));
    for(i=0; i<n; i++)
        mat[i] = (double*)malloc(m * sizeof(double));
    return mat;
}

void printDoubleMatrix(double const* const* mat, int n, int m)
{
    int i,j;
    printf("\n");
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
            printf("%lf ", mat[i][j]);
        printf("\n");
    }
}

void printDoubleVector(double const* vec, int n)
{
    int i,j;
    printf("\n");
    for(i=0; i<n; ++i)
        printf("%lf ", vec[i]);
    printf("\n");
}

void printSystemInequalities(double const* const* lhs, double const* rhs, int numConst, int numVars)
{
    printf("\n");
    for(int i=0; i<numConst; ++i)
    {
        for(int j=0; j<numVars; ++j)
            printf(" %lf ", lhs[i][j]);
        printf(" <= %lf \n",rhs[i]);
    }
}

void setDoubleMatrixToZero(double** mat, int n, int m)
{
    for(int i=0; i<n; ++i)
        for(int j=0; j<m; ++j)
            mat[i][n] = 0.0;
}

double fastAbs(double f)
{
    if( f < 0.0)
        return -f;
    return f;
}

// true if abs(b) < eps
int absValLT(double v, double eps_)
{
    return ( v < eps_ && -v < eps_);
}

double determinant3x3(double mat[3][3])
{
    double det = 0.0;
    det += mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]);
    det -= mat[1][0] * (mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2]);
    det += mat[2][0] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
    return det;
}

double determinant3vec3(double const* v0, double const* v1, double const* v2)
{
    double det = 0.0;
    det += v0[0] * (v1[1] * v2[2] - v1[2] * v2[1]);
    det -= v0[1] * (v1[0] * v2[2] - v1[2] * v2[0]);
    det += v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);
    return det;
}

double int_pow(double d, int e)
{
    double p = 1.0;
    for(int k=0; k<e; ++k)
        p *= d;
    return p;
}

int int_pow(int i, int e)
{
    int p = 1.0;
    for(int k=0; k<e; ++k)
        p *= i;
    return p;
}

bool iszeroeps(double d, double eps)
{
    if( fastAbs(d) > eps)
        return false;
    return true;
}

double numericzero(double d, double eps)
{
    if( fastAbs(d) > eps)
        return d;
    return 0.0;
}

void swapints(int* a, int* b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

