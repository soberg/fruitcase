#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>                // for gettimeofday()



void timerGetTime(timeval* timeInterval);
double getTimeMilliseconds(timeval timeStart, timeval timeEnd);
void freeDoubleMatrix(double** mat, int m);
double** allocateDoubleMatrix(int n, int m);
void printDoubleMatrix(double const* const* mat, int n, int m);
void printDoubleVector(double const* vec, int n);
void printSystemInequalities(double const* const* lhs, double const* rhs, int numConst, int numVars);
double fastAbs(double f);
int absValLT(double v, double eps_);
void setDoubleMatrixToZero(double** mat, int n, int m);
double determinant3x3(double mat[3][3]);
double determinant3vec3(double const* v0, double const* v1, double const* v2);
double numericzero(double d, double eps);
double int_pow(double d, int e);
int int_pow(int i, int e);
bool iszeroeps(double d, double eps);
void swapints(int* a, int* b);
