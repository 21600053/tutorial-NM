/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 22-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"
#include "myMatrix.h"

// Matrix addition
Matrix   addMat(Matrix _A, Matrix _B)
{
    if (_A.rows != _B.rows || _A.cols != _B.cols) {
        printf("\n*************************************************");
        printf("\n  ERROR!!: dimension error at 'addMat' function");
        printf("\n*************************************************\n");
        return createMat(0, 0);
    }

    Matrix Out = createMat(_A.rows, _B.cols);
    for (int i = 0; i < _A.rows; i++)
        for (int j = 0; j < _B.cols; j++)
            Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

    return Out;
}

//void LUdecomp(Matrix _A, Matrix* L, Matrix* U, Matrix* P)
//{
//   int M = _A.rows;
//   int N = _A.cols;
//   Matrix _P = diagonal(M, N); // 모든 pivoting을 저장할 P
//   Matrix _D = diagonal(M, N);
//   Matrix _L = zeros(M, N);
//   
//   if (M == N) // if matrix is not square, print error message 
//   {
//      double m = 0;
//      for (int k = 0; k < M - 1; k++)
//      {
//         printf("k = %d\n\n", k);
//         int m_row = 0; // maximum row num
//         double MP = 0; // maximum pivot
//         Matrix Pmat = diagonal(M, N); //A에 곱해질 p
//
//         for (int x = k; x < M; x++)      // 행에서 가장 큰 값 찾기
//         {
//            double MV = 0; // max value
//            int m_col = 0; // maximum col num
//
//            for (int y = 0; y < M; y++)   //열에서 가장 큰 값
//            {
//               if (abs(_A.at[x][y]) > abs(MV))
//               {
//                  MV = _A.at[x][y]; // 절대값이 가장 큰 값을 MV에 저장
//                  m_col = y; // 그 때의 col num을 저장
//               }
//               else
//                  MV = MV;
//            }
//            //열을 최대값으로 나누고 그때의 최대 열 값을 찾아옴
//            if (abs(_A.at[x][k] / _A.at[x][m_col]) > abs(MP))
//            {
//               MP = abs(_A.at[x][k] / _A.at[x][m_col]); //나눈 값의 절대값이 가장 큰 값을 MP에 저장
//               m_row = x; // 그 때의 row num을 저장
//            }
//            else
//               MP = MP;
//         }
//
//         for (int t = 0; t < M; t++)
//         {
//            Pmat.at[k][t] = _D.at[m_row][t]; // 각 단계 K에서 pivoting 시킬 P를 구현
//            Pmat.at[m_row][t] = _D.at[k][t];
//         }
//         _A = multiMat(Pmat, _A);
//         _P = multiMat(Pmat, _P);
//         _L = multiMat(Pmat, _L);
//         printMat(Pmat, "Pmat");
//         printMat(_A, "Pivot A");
//         printMat(_P, "P");
//
//         for (int i = k + 1; i < M; i++) // row reduction
//         {
//            m = (_A.at[i][k] / _A.at[k][k]);
//            _L.at[i][k] = m;
//
//            for (int j = 0; j < N; j++)
//            {
//               _A.at[i][j] = _A.at[i][j] - m * _A.at[k][j];
//            }
//         }
//         printMat(_A, "PA");
//         printMat(_L, "L");
//         printf("------------------------------------------------------------------------------------\n");
//      }
//      *L = addMat(_L, _D);
//      *U = _A;
//      *P = _P;
//
//   }
//   else
//   {
//      printf("------------------------------------------------------------------------------------\n");
//      printf("              Matrix is not square           \n");
//      printf("------------------------------------------------------------------------------------\n");
//   }
//}

Matrix fwdsub(Matrix _U, Matrix _d)
{
    int M = _U.rows;
    Matrix _x = zeros(_d.rows, _d.cols); // 결과를 표시할 행렬을 만들어줌.
    if (_U.rows = _U.cols)
    {
        for (int k = 0; k < M; k++)
        {
            double sum = 0;
            for (int i = 0; i < k; i++)
            {
                sum = sum + _U.at[k][i] * _x.at[i][0];

            }
            _x.at[k][0] = (_d.at[k][0] - sum) / _U.at[k][k]; // 수식 표현
        }
    }
    else
    {
        printf("------------------------------------------------------------------------------------\n");
        printf("              Matrix is not square           \n");
        printf("------------------------------------------------------------------------------------\n");
    }
    return _x;
}

//Matrix solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b)
//{
//   _b = multiMat(_P, _b);  // LUx = Pb에서 Pb를 만들어줌
//   Matrix y = fwdsub(_L, _b); // Ly = Pb 형태로 만듬
//   Matrix x = backsub(_U, y); // y = Ux
//   printMat(x, "solution x");   
//   return x;
//}


Matrix getCols(Matrix _A, int i)
{
    Matrix out = zeros(_A.cols, 1);
    for (int j = 0; j < _A.cols; j++)
    {
        out.at[j][0] = _A.at[j][i];
    }

    return out;
}

//void inv(Matrix _A)
//{
//   Matrix _P = diagonal(_A.rows, _A.cols);
//   Matrix _L = zeros(_A.rows, _A.cols);
//   Matrix _U = copyMat(_A);
//   Matrix invA = zeros(_A.rows, _A.cols);
//   Matrix _D = diagonal(_A.rows, _A.cols);
//   //Matrix _b = zeros(_A.rows, 1);
//
//   int M = _A.cols; 
//   LUdecomp(_A, &_L, &_U, &_P); // LU decomposition을 통해 L, U, P의 값을 얻어낸다
//   
//   for (int j = 0; j < M; j++)
//   {
//      Matrix _b = getCols(_D, j); // 해당 열을 불러온다
//
//      Matrix x = solveLU(_L, _U, _P, _b); // 각 단계별로 solveLU를 통해 inverse matrix의 열 값을 찾아낸다
//
//      //printMat(x, "x");
//
//      for (int k = 0; k < M; k++) // 각 열에 inverse matrix 값을 채워준다.
//      {
//         invA.at[k][j] = x.at[k][0];
//      }
//      
//   }
//   printMat(invA, "invA");
//}

double vec_Norm(Matrix _v)      //vector L2_Norm
{
    int n = _v.rows;
    double sum = 0;
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        sum += pow(_v.at[i][0], 2);
    }
    norm = sqrt(sum);
    return norm;
}


double L1_Norm(Matrix _A)   //L1 norm - 열의 합 중 가장 큰 값
{
    int m = _A.rows;
    int n = _A.cols;
    double norm = 0;
    for (int i = 0; i < m; i++)
    {
        double sum = 0;

        for (int j = 0; j < n; j++)
        {
            sum += fabs(_A.at[j][i]);
        }

        if (norm < sum)
        {
            norm = sum;
        }
    }
    return norm;
}

double L2_Norm(Matrix _A)   // Euclidean Norm
{
    int m = _A.rows;
    int n = _A.cols;
    double sum = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sum += pow(_A.at[i][j], 2);
        }
    }
    double norm = sqrt(sum);
    return norm;
}

double inf_Norm(Matrix _A)   //infinite norm 행의 합중 가장 큰값
{
    int m = _A.rows;
    int n = _A.cols;
    double norm = 0;
    for (int i = 0; i < m; i++)
    {
        double sum = 0;

        for (int j = 0; j < n; j++)
        {
            sum += fabs(_A.at[i][j]);
        }

        if (norm < sum)
        {
            norm = sum;
        }
    }
    return norm;
}

Matrix multiVec(Matrix _A, Matrix _B)
{
    Matrix Out = zeros(_A.rows, _B.cols);
    for (int i = 0; i < _A.rows; i++)
    {
        for (int j = 0; j < _A.rows; j++)
        {
            Out.at[i][j] = _A.at[i][0] * _B.at[0][j];

        }
    }

    return Out;
}

void trans(Matrix _A, Matrix _B)
{
    for (int i = 0; i < _A.cols; i++)
    {
        for (int j = 0; j < _A.rows; j++)
        {
            _B.at[i][j] = _A.at[j][i];
        }
    }
}


void QR_householder(Matrix _A)
{
    int n = _A.rows;      //A: nxn matrix
    Matrix R = copyMat(_A);
    Matrix Q = diagonal(n, n);

    Matrix v_trans = zeros(1, n);
    Matrix H = zeros(n, n);
    Matrix q = diagonal(n, n);
    Matrix r = copyMat(_A);

    for (int k = 0; k < n; k++)
    {
        Matrix I = diagonal(n, n);
        Matrix c = zeros(n, 1);
        Matrix e = zeros(n, 1);
        Matrix v = zeros(n, 1);
        for (int i = k; i < n; i++)
        {
            c.at[i][0] = R.at[i][k];
        }

        if (c.at[k][0] > 0)
            e.at[k][0] = 1;
        else
            e.at[k][0] = -1;

        for (int t = 0; t < n; t++)
        {
            v.at[t][0] = c.at[t][0] + vec_Norm(c) * e.at[t][0];
        }

        trans(v, v_trans);

        Matrix vtv = zeros(1, 1);
        Matrix vvt = zeros(n, n);

        multiMat(v_trans, v, vtv);
        multiMat(v, v_trans, vvt);
        if (vtv.at[0][0] == 0)
        {
            vtv.at[0][0] = 0.0001;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    H.at[i][j] = I.at[i][j] - ((2 / vtv.at[0][0]) * vvt.at[i][j]);

                }
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    H.at[i][j] = I.at[i][j] - ((2 / vtv.at[0][0]) * vvt.at[i][j]);

                }
            }
        }

        multiMat(Q, H, q);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Q.at[i][j] = q.at[i][j];
            }
        }

        multiMat(H, R, r);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                R.at[i][j] = r.at[i][j];
            }
        }

    }
    multiMat(r, q, _A);

    printMat(q, "Q");
    printMat(r, "R");
}

Matrix eig(Matrix _A)
{
    int n = _A.rows;

    Matrix lamda = zeros(n, 1);

    for (int i = 0; i < 50; i++)
    {
        QR_householder(_A);

    }

    for (int k = 0; k < n; k++)
    {
        for (int j = 0; j < n; j++)
        {
            if (k == j)
            {
                lamda.at[k][0] = _A.at[k][j];
            }
        }
    }
    //printMat(lamda, "lamda");
    return lamda;
}

double cond(Matrix _A)
{
    int m = _A.rows;
    int n = _A.cols;
    Matrix U = zeros(n, n);
    Matrix trans_A = zeros(n, m);

    trans(_A, trans_A);
    multiMat(trans_A, _A, U);
    printMat(U, "U");
    Matrix lamda = zeros(U.rows, 1);
    lamda = eig(U);
    printMat(lamda, "lamda");

    double max_temp = lamda.at[0][0];
    double min_temp = lamda.at[0][0];
    double max_sig = 0;
    double min_sig = 0;
    double sigma = 0;

    for (int i = 0; i < n; i++)
    {

        if (max_temp < lamda.at[i][0])
        {
            if (lamda.at[i][0] != 0)
            {
                max_temp = lamda.at[i][0];
            }
        }
        else
        {
            max_temp = max_temp;
        }

        if (min_temp > lamda.at[i][0])
        {
            if (lamda.at[i][0] != 0)
            {
                min_temp = lamda.at[i][0];
            }
        }
        else
        {
            min_temp = min_temp;
        }
    }

    max_sig = sqrt(max_temp);
    min_sig = sqrt(min_temp);
    sigma = max_sig / min_sig;

    printf("\n\nsigma = %f\n", sigma);
    return sigma;
}

Matrix   linearFit(Matrix _x, Matrix _y) {
    int mx = _x.rows;
    int my = _y.rows;

    double a1 = 0;
    double a0 = 0;

    double Sx = 0;
    double Sxx = 0;
    double Sxy = 0;
    double Sy = 0;

    if (mx != my || mx == 1)
    {
        printf("Error: The number of elements in x must be equal to y and more than 1\n");
    }
    else
    {
        int m = mx;

        for (int k = 0; k < m; k++)
        {
            Sx = Sx + _x.at[k][0];
            Sxx = Sxx + pow(_x.at[k][0], 2);
            Sxy = Sxy + _x.at[k][0] * _y.at[k][0];
            Sy = Sy + _y.at[k][0];
        }// end
        double den = m * Sxx - pow(Sx, 2);
        a1 = (m * Sxy - Sx * Sy) / den;
        a0 = (Sxx * Sy - Sxy * Sx) / den;
    }


    double z_array[] = { a1, a0 };
    return arr2Mat(z_array, 2, 1);
}

Matrix   arr2Mat(double* _1Darray, int _rows, int _cols)
{
    Matrix Output = createMat(_rows, _cols);

    for (int i = 0; i < _rows; i++)
        for (int j = 0; j < _cols; j++)
            Output.at[i][j] = _1Darray[i * _cols + j];

    return Output;
}
Matrix   linearInterp(Matrix _x, Matrix _y, Matrix _xq)
{
    int temp = 0;
    int mx = _x.rows;
    int my = _y.rows;
    int m = _xq.rows;
    Matrix yq = createMat(_xq.rows, _xq.cols);
    initMat(yq, 0);
    if ((mx != my) || mx == 1 || my == 1)
    {
        printf("Error : The number of elements in x must be equal to y and more than 1");
    }

    else if (m == 0)
    {
        printf("Error : Input Query is not zero");
    }
    // lagrange form for 1st order polynomial
    else {
        for (int i = 0; i < mx - 1; i++)
        {
            for (int j = temp; j < m; j++)
            {
                if ((_x.at[i][0] <= _xq.at[j][0]) && (_x.at[i + 1][0] >= _xq.at[j][0]))
                {
                    yq.at[j][0] = _y.at[i][0] * (_xq.at[j][0] - _x.at[i + 1][0]) / (_x.at[i][0] - _x.at[i + 1][0]) + _y.at[i + 1][0] * (_xq.at[j][0] - _x.at[i][0]) / (_x.at[i + 1][0] - _x.at[i][0]);
                    temp++;
                }
            }
        }
    }
    return yq;
}


Matrix   gradient(Matrix _x, Matrix _y)
{
    int n = _x.rows;
    Matrix df = createMat(n, 1);

    double h = _x.at[1][0] - _x.at[0][0];

    //assumption n>2
    if (n > 2)
    {
        // 1. 3-point FW
        df.at[0][0] = ((-3 * _y.at[0][0]) + (4 * _y.at[1][0]) - _y.at[2][0]) / (2 * h);

        //2. 2-point cetnral diff 
        for (int i = 1; i < n - 1; i++)
        {
            df.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);
        }
        //3 3-point BW
        df.at[n - 1][0] = ((_y.at[n - 3][0]) + (-4 * _y.at[n - 2][0]) + (3 * _y.at[n - 1][0])) / (2 * h);
    }
    else
    {
        df.at[0][0] = 0; //2point fwd
        df.at[1][0] = 1; // 2point bwd
    }
    return df;
}

Matrix   gradientFunc(double func(const double x), Matrix xin) {

    int n = xin.rows;
    Matrix y = createMat(n, 1);

    // define y[0] to y[n-1]
    for (int i = 0; i < n; i++)
    {
        y.at[i][0] = func(xin.at[i][0]);
    }
    //numerical differentiation

    return gradient(xin, y);

}

void gradient1D(double x[], double y[], double dydx[], int m)
{

    double h = x[1] - x[0];
    //assumption n>2
    if (m > 2)
    {
        // 1. 3-point FW
        dydx[0] = ((-3 * y[0]) + (4 * y[1]) - y[2]) / (2 * h);

        //2. 2-point cetnral diff
        for (int i = 1; i < m - 1; i++)
        {
            dydx[i] = (y[i + 1] - y[i - 1]) / (2 * h);
        }
        //3 3-point BW
        dydx[m - 1] = ((y[m - 3]) + (-4 * y[m - 2]) + (3 * y[m - 1])) / (2 * h);
    }
    else
    {
        dydx[0] = 0; //2point fwd
        dydx[1] = 1; // 2point bwd
    }

}



double NewtonRaphson(double myFunc(const double x), double mydFunc(const double x), double _x0, double _tol)
{
    int k = 0;
    int Nmax = 1000;
    double _x = _x0;  //초기값
    double xn = _x0;
    double ep = 0;
    double hk = 0;
    printf("------------------------------------------------------------------------------------\n");
    printf("                 Newton Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");
    do {
        printf("Iteration:%d \t", k);
        printf("X(n): %f \t", xn);
        printf("Tolerance: %.10f\n", ep);
        hk = -(myFunc(_x) / mydFunc(_x));
        xn = _x + hk;

        if (myFunc(xn) == 0)
        {
            break;
        }
        ep = fabs(hk);

        _x = xn;
        k++;

    } while (k<Nmax && fabs(ep)>_tol);

    return xn;
}

double trapz(double _x[], double _y[], int _m)
{
    int N = _m - 1;
    double I = 0;
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum += (_y[i] + _y[i + 1]);
        I = 0.5 * sum * (_x[i + 1] - _x[i]);
    }
    return I;
}

double integral(double func(const double x), double a, double b, int n)
{
    double I = 0;
    double h = (b - a) / n;
    double sum1 = 0;
    double sum2 = 0;
    for (int i = 1; i < n; i += 2)
    {
        double xi = a + h * i;
        sum1 += func(xi);
        /*	printf("i = %d\n", i);*/
    }
    for (int k = 2; k < n - 1; k += 2)
    {
        double xk = a + h * k;
        sum2 += func(xk);
        //printf("k = %d\n", k);
    }
    I = (h / 3) * (func(a) + 4 * sum1 + 2 * sum2 + func(b));
    return I;
}