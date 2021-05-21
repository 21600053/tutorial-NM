/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"



extern	Matrix	addMat(Matrix _A, Matrix _B);

//extern void LUdecomp(Matrix _A, Matrix* _L, Matrix* _U, Matrix* _P);
//
//extern Matrix solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b);

extern Matrix fwdsub(Matrix _U, Matrix _d);

extern Matrix getCols(Matrix _A, int i);

extern void inv(Matrix _A);

extern double L1_Norm(Matrix _A);

extern double L2_Norm(Matrix _A);

extern double inf_Norm(Matrix _A);

extern double vec_Norm(Matrix _v);

extern Matrix multiVec(Matrix _A, Matrix _B);

extern void QR_householder(Matrix _A);

extern void trans(Matrix _A, Matrix _B);

extern Matrix eig(Matrix A);

extern double cond(Matrix _A);

extern Matrix	linearFit(Matrix _x, Matrix _y);

extern Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);

extern Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq);

extern Matrix	gradient(Matrix _x, Matrix _y);

extern Matrix	gradientFunc(double func(const double x), Matrix xin);

extern void gradient1D(double x[], double y[], double dydx[], int m);

extern double NewtonRaphson(double myFunc(const double x), double mydFunc(const double x), double _x0, double _tol);

extern double trapz(double _x[], double _y[], int _m);

extern double integral(double func(const double x), double a, double b, int n);
#endif