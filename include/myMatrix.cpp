/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 03-04-2021
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix.h"
#include "myNM.h"


// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	return Out;
}

// Free a memory allocated matrix
void	freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}

// Print matrix
void	printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.6f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}


// initialization of Matrix elements
void	initMat(Matrix _A, double _val)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = _val;

}

// Create matrix of all zeros
Matrix	zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out, 0);

	return Out;
}

void	gaussElim(Matrix _A, Matrix _b)
{
	int M = _A.rows;
	int N = _A.cols;
	if (M == N) // if matrix is not squar, print error message 
	{
		double m = 0;
		for (int k = 0; k < M - 1; k++) // 행 개수 -1 개의 pivot
		{
			for (int i = k + 1; i < M; i++) // 두번째 행 부터 계산
			{
				m = (_A.at[i][k] / _A.at[k][k]);
				for (int j = 0; j < N; j++) // 첫번째 열부터 열의 갯수만큼 계산
				{
					_A.at[i][j] = _A.at[i][j] - m * _A.at[k][j];

				}
				_b.at[i][0] = _b.at[i][0] - m * _b.at[k][0];

			}
		}
		printf("------------------------------------------------------------------------------------\n");
		printf("				gaussElim_matU result             \n");
		printf("------------------------------------------------------------------------------------\n");
		printMat(_A, "gaussElim_U");
		printf("------------------------------------------------------------------------------------\n");
		printf("				gaussElim_vecd result             \n");
		printf("------------------------------------------------------------------------------------\n");

		printMat(_b, "gaussElim_d");
	}
	else
	{
		printf("------------------------------------------------------------------------------------\n");
		printf("			     Matrix is not square			  \n");
		printf("------------------------------------------------------------------------------------\n");
	}
}

Matrix backsub(Matrix _U, Matrix _d)
{
	int M = _U.rows;
	Matrix _x = zeros(_U.cols, 1); // 결과를 표시할 행렬을 만들어줌.
	if (_U.rows = _U.cols)
	{

		for (int k = M - 1; k >= 0; k--)
		{
			double sum = 0;
			for (int i = k + 1; i < M; i++)
			{
				sum = sum + _U.at[k][i] * _x.at[i][0];

			}
			_x.at[k][0] = (_d.at[k][0] - sum) / _U.at[k][k]; // 수식 표현
		}

	}
	else
	{
		printf("------------------------------------------------------------------------------------\n");
		printf("			     Matrix is not square			  \n");
		printf("------------------------------------------------------------------------------------\n");
	}
	return _x;
}


//Matrix multiMat(Matrix _A, Matrix _B)
//{
//	Matrix Out = zeros(_A.rows, _B.cols);
//	for (int i = 0; i < _A.rows; i++)
//	{
//		for (int j = 0; j < _A.rows; j++)
//		{
//			Out.at[i][j] = 0;
//			for (int k = 0; k < _A.rows; k++)
//			{
//				Out.at[i][j] += _A.at[i][k] * _B.at[k][j];
//			}
//		}
//	}
//	
//	return Out;
//}
void multiMat(Matrix _A, Matrix _b, Matrix _Out)
{
	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _b.cols; j++)
		{
			_Out.at[i][j] = 0;
			for (int k = 0; k < _A.cols; k++)
			{
				_Out.at[i][j] += _A.at[i][k] * _b.at[k][j];
			}
		}
	}

}


Matrix   copyMat(Matrix _A)
{
	Matrix Out = createMat(_A.rows, _A.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] = _A.at[i][j];

	return Out;
}

Matrix diagonal(int _rows, int _cols)
{
	Matrix Out = zeros(_rows, _cols);
	for (int i = 0; i < _rows; i++)
	{
		for (int j = 0; j < _cols; j++)
		{
			if (i == j)
			{
				Out.at[i][j] = +1;
			}
		}
	}

	return Out;
}

//double deflection(double _x)
//{
//    double L = 4;
//    double E = 70 * pow(10, 9);
//    double I = 52.9 * pow(10, -6);
//    double w_0 = 20 * pow(10, 3);
//    double F = (w_0 * _x * (7 * pow(L, 4) - 10 * pow(L, 2) * pow(_x, 2) + 3 * pow(_x, 4))) / (360 * E * I * L); // 주어진 함수 y
//    return F;
//}
//
//double func(double _x)
//{
//    double L = 4;
//    double E = 70 * pow(10, 9);
//    double I = 52.9 * pow(10, -6);
//    double w_0 = 20 * pow(10, 3);
//    double F = (w_0 * (7 * pow(L, 4) - 10 * pow(L, 2) * pow(_x, 2) + 3 * pow(_x, 4))) / (360 * E * I * L) - (w_0 * _x * (20 * pow(L, 2) * _x - 12 * pow(_x, 3))) / (360 * E * I * L); // 주어진 함수 y를 미분한 f(x)
//    return F;
//}
//
//
//double dfunc(double _x)
//{
//    double L = 4;
//    double E = 70 * pow(10, 9);
//    double I = 52.9 * pow(10, -6);
//    double w_0 = 20 * pow(10, 3);
//    double F = -(w_0 * (20 * pow(L, 2) * _x - 12 * pow(_x, 3))) / (180 * E * I * L) - (w_0 * _x * (20 * pow(L, 2) - 36 * pow(_x, 2))) / (360 * E * I * L); //f(x)를 미분
//    return F;
//
//}

//double NewtonRaphson(double _x0, double _tol)
//{
//    int k = 0;
//    int Nmax = 1000;
//    double _x = _x0;  //초기값
//    double xn = _x0;
//    double ep = 0;
//    double hk = 0;
//    printf("------------------------------------------------------------------------------------\n");
//    printf("                 Newton Method Results             \n");
//    printf("------------------------------------------------------------------------------------\n");
//    do {
//        printf("Iteration:%d \t", k);
//        printf("X(n): %f \t", xn);
//        printf("Tolerance: %.10f\n", ep);
//        hk = -(func(_x) / dfunc(_x));
//        xn = _x + hk;
//
//        if (func(xn) == 0) // 기울기가 0이 될 경우 발산하기 때문에 멈춰준다
//        {
//            break;
//        }
//        ep = fabs(hk);
//
//        _x = xn;
//        k++;
//
//    } while (k<Nmax && fabs(ep)>_tol);
//    printf("Solution: %f \n", xn);
//    printf("Max deflection: %f \n", deflection(xn));
//    return xn;
//}
//
//double bisectionNL(double _a0, double _b0, double _tol)
//{
//    int k = 0;
//    int Nmax = 1000;
//    double a = _a0;
//    double b = _b0;
//    double xn = 0;
//    double ep = 1000;
//    printf("------------------------------------------------------------------------------------\n");
//    printf("               Bisection Method Results             \n");
//    printf("------------------------------------------------------------------------------------\n");
//
//    do {
//        xn = (a + b) / 2;
//        ep = fabs(func(xn));
//        printf("Iteration:%d \t", k);
//        printf("X(n): %f \t", xn);
//        printf("Tolerance: %.10f\n", ep);
//
//        if (func(a) * func(xn) < 0)
//            b = xn;
//        else
//            a = xn;
//        k++;
//    } while (k<Nmax && ep>_tol);
//    printf("Solution: %f \n", xn);
//    printf("Max deflection: %f \n", deflection(xn));
//
//    return xn;
//}
//
//double B_NewtonRaphson(double _x0, double _tol)
//{
//    int k = 0;
//    int Nmax = 1000;
//    double _x = _x0;  //초기값
//    double xn = _x0;
//    double ep = 0;
//    double hk = 0;
//    printf("------------------------------------------------------------------------------------\n");
//    printf("                 Newton Method Results             \n");
//    printf("------------------------------------------------------------------------------------\n");
//    do {
//        printf("Iteration:%d \t", k);
//        printf("X(n): %f \t", xn);
//        printf("Tolerance: %.10f\n", ep);
//        hk = -(B_func(_x) / B_dfunc(_x));
//        xn = _x + hk;
//
//        if (B_func(xn) == 0)
//        {
//            break;
//        }
//        ep = fabs(hk);
//
//        _x = xn;
//        k++;
//
//    } while (k<Nmax && fabs(ep)>_tol);
//    printf("Solution: %f \n", xn);
//    printf("Max deflection: %f \n", deflection(xn));
//    return xn;
//}
//
//double B_func(double _x)
//{
//    double F = (1 / _x) - 2;
//    return F;
//}
//
//double B_dfunc(double _x)
//{
//    double F = -(1 / pow(_x, 2));
//    return F;
//}
//
//double hybrid(double _a0, double _b0, double _x0, double _tol)
//{
//    int k = 0;
//    int Nmax = 1000;
//    double a = _a0;
//    double b = _b0;
//    double ep = 0;
//    double _x = _x0;  //init
//    double xn = 0;
//    double hk = 0;
//
//    do {
//
//        if (xn < a || xn > b) // 범위를 벗어날경우 bisection
//        {
//            xn = (a + b) / 2;
//            ep = fabs(B_func(xn));
//            printf("Iteration:%d \t", k);
//            printf("X(n): %f \t", xn);
//            printf("Tolerance: %.10f\n", ep);
//
//            if (B_func(a) * B_func(xn) < 0)
//                b = xn;
//            else
//                a = xn;
//            printf("bisection\n\n");
//        }
//        else if (fabs(B_dfunc(xn)) < 0.01)  //발산할 경우 bisection
//        {
//            xn = (a + b) / 2;
//            ep = fabs(B_func(xn));
//            printf("Iteration:%d \t", k);
//            printf("X(n): %f \t", xn);
//            printf("Tolerance: %.10f\n", ep);
//
//            if (B_func(a) * B_func(xn) < 0)
//                b = xn;
//            else
//                a = xn;
//            printf("bisection\n\n");
//        }
//        else                      // otherwise use Newton
//        {
//            printf("Iteration:%d \t", k);
//            printf("X(n): %f \t", xn);
//            printf("Tolerance: %.10f\n", ep);
//            hk = -(B_func(_x) / B_dfunc(_x));
//            xn = _x + hk;
//
//            if (B_func(xn) == 0)
//            {
//                break;
//            }
//            ep = fabs(hk);
//
//            _x = xn;
//            printf("newton\n\n");
//        }
//        k++;
//
//    } while (k<Nmax && ep>_tol);
//    printf("Solution: %f \n", xn);
//
//    return xn;
//}