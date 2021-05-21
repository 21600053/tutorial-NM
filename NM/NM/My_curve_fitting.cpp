///*-------------------------------------------------------------------------------\
//@ Numerical Methods by Young-Keun Kim - Handong Global University
//
//Author           : [YOUR NAME]
//Created          : 03-05-2021
//Modified         : 03-05-2021
//Language/ver     : C++ in MSVS2019
//
//Description      : [Tutorial]curve_fitting.cpp
//-------------------------------------------------------------------------------*/
//
//#include "myNM.h"
//
//// Returns the parameters of the linear least square function.
//Matrix	linearFit(Matrix _x, Matrix _y);
//
//// Create a matrix from 1D-array
//Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);
//
//Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq);
//
//int main(int argc, char* argv[])
//{
//	int M = 20;
//	double T_array[] = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
//	double P_array[] = { 0.94, 0.96, 1.0, 1.05, 1.07, 1.09, 1.14, 1.17, 1.21, 1.24, 1.28 };
//	double xq_array[] = { 0, 5, 10, 15,  20, 25, 30, 35, 40, 45,  50, 55, 60, 65, 70, 75,  80, 85,  90, 95, 100 };
//	//double xq_array[] = { 5, 15, 25, 35, 45, 55, 65, 75, 85, 95};
//	
//	Matrix T = arr2Mat(T_array, M, 1);
//	Matrix P = arr2Mat(P_array, M, 1);
//	Matrix xq = arr2Mat(xq_array, M, 1);
//	//Matrix z = linearFit(T, P);
//	Matrix z = linearInterp(T, P, xq);
//
//
//	printMat(T, "T");
//	printMat(P, "P");
//	printMat(z, "z");
//	
//	double ans_100 = 0;
//	double ans_150 = 0;
//	ans_100 = z.at[1][0] + z.at[0][0] * 100;
//	ans_150 = z.at[1][0] + z.at[0][0] * 150;
//	printf("------------------------------------------------------------------------------------\n");
//	printf("Predict the pressure at (T = 100C)\n = %f\n", ans_100);
//	printf("Predict the pressure at (T = 150C)\n = %f\n", ans_150);
//	printf("------------------------------------------------------------------------------------\n");
//
//
//
//
//	system("pause");
//	return 0;
//}
//
//// Returns the parameters of the linear least square function.
//Matrix	linearFit(Matrix _x, Matrix _y) {
//	int mx = _x.rows;
//	int my = _y.rows;
//
//	double a1 = 0;
//	double a0 = 0;
//
//	double Sx = 0;
//	double Sxx = 0;
//	double Sxy = 0;
//	double Sy = 0;
//	
//	if (mx != my || mx == 1)
//	{
//		printf("Error: The number of elements in x must be equal to y and more than 1\n");
//	}
//	else
//	{
//		int m = mx;
//
//		for (int k = 0; k < m; k++)
//		{
//			Sx = Sx + _x.at[k][0];
//			Sxx = Sxx + pow(_x.at[k][0], 2);
//			Sxy = Sxy + _x.at[k][0] * _y.at[k][0];
//			Sy = Sy + _y.at[k][0];
//		}// end
//		double den = m * Sxx - pow(Sx, 2);
//		a1 = (m * Sxy - Sx * Sy) / den;
//		a0 = (Sxx * Sy - Sxy * Sx) / den;
//	}
//
//
//	double z_array[] = { a1, a0 };
//	return arr2Mat(z_array, 2, 1);
//}
//
//// Create a matrix from 1D-array
//Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
//{
//	Matrix Output = createMat(_rows, _cols);
//
//	for (int i = 0; i < _rows; i++)
//		for (int j = 0; j < _cols; j++)
//			Output.at[i][j] = _1Darray[i * _cols + j];
//
//	return Output;
//}
//
//
//Matrix	linearInterp(Matrix _x, Matrix _y, Matrix _xq)
//{
//	int temp = 0;
//	int mx = _x.rows;
//	int my = _y.rows;
//	int m = _xq.rows;
//	Matrix yq = createMat(_xq.rows, _xq.cols);
//	initMat(yq, 0);
//	if ((mx != my) || mx == 1 || my == 1)
//	{
//		printf("Error : The number of elements in x must be equal to y and more than 1");
//	}
//	
//	else if (m == 0)
//	{
//		printf("Error : Input Query is not zero");
//	}
//	// lagrange form for 1st order polynomial
//	else {
//		for (int i = 0; i < mx - 1; i++)
//		{
//			for (int j = temp; j < m; j++)
//			{
//				if ((_x.at[i][0] <= _xq.at[j][0]) && (_x.at[i + 1][0] >= _xq.at[j][0]))
//				{
//					yq.at[j][0] = _y.at[i][0] * (_xq.at[j][0] - _x.at[i + 1][0]) / (_x.at[i][0] - _x.at[i + 1][0]) + _y.at[i + 1][0] * (_xq.at[j][0] - _x.at[i][0]) / (_x.at[i + 1][0] - _x.at[i][0]);
//					temp++;
//				}
//			}
//		}
//	}
//	return yq;
//}