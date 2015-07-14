//
//  Matrix.h
//  testProject
//
//  Created by Sören Lennart Berg on 10/2/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#ifndef __testProject__RealMatrix__
#define __testProject__RealMatrix__

#include <iostream>
#include <assert.h>

double** allocateMatrix(int rows, int columns);
void freeMatrix(double** A, int rows);

class RealMatrix
{
public:
    RealMatrix();
    RealMatrix(const RealMatrix& A);
    RealMatrix(int _m);
    RealMatrix(int _m, int _n);
    ~RealMatrix();
    
    RealMatrix& operator= (const RealMatrix& other);
    
    double& operator()(int i, int j) { return entries[ I(i,j) ]; }
    double  operator()(int i, int j) const { return entries[ I(i,j) ]; }
    
    void resize(int _m, int _n);
    
    void zero(void);
    void one(void);
    void set_all(double d);
    void identity(void);
    void identity(int _m);
    
    int rows(void) const { return m; }
    int columns(void) const { return n; }
    
    
    /*/ member access
    operator const double* () const;
    operator double* ();
    const double* operator[] (int row) const;
    double* operator[] (int row);
    double operator() (int iRow, int iCol) const;
    double& operator() (int iRow, int iCol);
    void SetRow (int iRow, const Vector4<Real>& rkV);
    Vector GetRow (int iRow) const;
    void SetColumn (int iCol, const Vector4<Real>& rkV);
    Vector4<Real> GetColumn (int iCol) const;
    void GetColumnMajor (Real* afCMajor) const;
    
    // assignment
    RealMatrix4& operator= (const Matrix4& rkM);
    
    // comparison
    bool operator== (const Matrix4& rkM) const;
    bool operator!= (const Matrix4& rkM) const;
    bool operator<  (const Matrix4& rkM) const;
    bool operator<= (const Matrix4& rkM) const;
    bool operator>  (const Matrix4& rkM) const;
    bool operator>= (const Matrix4& rkM) const;
    
    // arithmetic operations
    Matrix4 operator+ (const Matrix4& rkM) const;
    Matrix4 operator- (const Matrix4& rkM) const;
    Matrix4 operator* (const Matrix4& rkM) const;
    Matrix4 operator* (Real fScalar) const;
    Matrix4 operator/ (Real fScalar) const;
    Matrix4 operator- () const;
    
    // arithmetic updates
    Matrix4& operator+= (const Matrix4& rkM);
    Matrix4& operator-= (const Matrix4& rkM);
    Matrix4& operator*= (Real fScalar);
    Matrix4& operator/= (Real fScalar);
    
    // matrix times vector
    Vector4<Real> operator* (const Vector4<Real>& rkV) const;  // M * v
    
    // other operations
    Matrix4 Transpose () const;  // M^T
    Matrix4 TransposeTimes (const Matrix4& rkM) const;  // this^T * M
    Matrix4 TimesTranspose (const Matrix4& rkM) const;  // this * M^T
    Matrix4 Inverse () const;
    Matrix4 Adjoint () const;
    Real Determinant () const;
    Real QForm (const Vector4<Real>& rkU,
                const Vector4<Real>& rkV) const;  // u^T*M*v */
    
private:
    int I (int _m, int _n) const;
    int I(int _m, int _n, int a, int b) const;
    
    int m, n;
    double* entries;
};


#endif /* defined(__testProject__Matrix__) */
