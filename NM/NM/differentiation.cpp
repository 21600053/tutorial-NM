/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park
Created          : 10-05-2021
Modified         : 10-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Differentiation_student.cpp
-------------------------------------------------------------------------------*/

//#include "myNM.h"
#include "../../include/myNM.h"
// Return the dy/dx results for the input data. (truncation error: O(h^2))
Matrix   gradient(Matrix _x, Matrix _y);

// Define a function that defines the target equation.
double myFunc(const double x);
double mydFunc(const double x);
// Return the dy/dx results for the target equation. (truncation error: O(h^2))
Matrix   gradientFunc(double func(const double x), Matrix xin);

void gradient1D(double x[], double y[], double dydx[], int m);

double NewtonRaphson(double myFunc(const double x), double mydFunc(const double x), double _x0, double _tol);

int main(int argc, char* argv[])
{

    // PART 1
    printf("\n**************************************************");
    printf("\n|                     PART 1.                    |");
    printf("\n**************************************************\n");

    Matrix t = txt2Mat("", "Q1_vect");
    Matrix pos = txt2Mat("", "Q1_vecx");  //position x

    Matrix vel = gradient(t, pos); //velocity
    Matrix acc = gradient(t, vel); // acceleration

    /*printMat(t, "t");

    printMat(pos, "pos");*/
    printMat(vel, "vel");
    printMat(acc, "acc");

    double x[] = { 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4 };
    double y[] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2, 1.1, 0.23, -0.59 };
    double _dydx[] = { 0 };
    double _dydx2[] = { 0 };
    int len = sizeof(x) / sizeof(double);

    gradient1D(x, y, _dydx, len);
    
    printf("\n**************************************************\n");
    printf("velocity : 1D gradient\n");
    for (int i = 0; i < len; i++)
    {
       printf("%f\n", _dydx[i]);
    }  
    printf("\n**************************************************\n");
    gradient1D(x, _dydx, _dydx2, len);
    printf("acceleration : 1D gradient\n");
    for (int i = 0; i < len; i++)
    {
        printf("%f\n", _dydx2[i]);
    }
    printf("\n**************************************************\n");

    // PART 2
    printf("\n**************************************************");
    printf("\n|                     PART 2.                    |");
    printf("\n**************************************************\n");

    Matrix xin = txt2Mat("", "Q2_vecxin");
    Matrix dydx = gradientFunc(myFunc, xin);

    printMat(xin, "xin");
    printMat(dydx, "dydx");

    double result;
    double x0 = 2;
    double tol = 0.00001;
    result = NewtonRaphson(myFunc, mydFunc, x0, tol);


    system("pause");
    return 0;


}

//Gradient

// Return the dy/dx results for the input data. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
//Matrix   gradient(Matrix _x, Matrix _y)
//{
//   int m = _x.rows;
//   Matrix df = createMat(m, 1);
//
//   double h = _x.at[1][0] - _x.at[0][0];
//
//   //assumption n>2
//   if (m > 2)
//   {
//      // 1. 3-point FW
//      df.at[0][0] = ((-3 * _y.at[0][0]) + (4 * _y.at[1][0]) - _y.at[2][0]) / (2 * h);
//
//      //2. 2-point cetnral diff
//      for (int i = 1; i < m - 1; i++)
//      {
//         df.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);
//      }
//      //3 3-point BW
//      df.at[m - 1][0] = ((_y.at[m - 3][0]) + (-4 * _y.at[m - 2][0]) + (3 * _y.at[m - 1][0])) / (2 * h);
//   }
//   else
//   {
//      df.at[0][0] = 0; //2point fwd
//      df.at[1][0] = 1; // 2point bwd
//   }
//   return df;
//}

// Define a function that defines the target equation.
double myFunc(const double x) {
    return  x * x * x;
}

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp

/*
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
*/

/*
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
*/

double mydFunc(const double x)
{
    return 3 * x * x;
}



//double NewtonRaphson(double myFunc(const double x), double mydFunc(const double x), double _x0, double _tol)
//{
//    int k = 0;
//    int Nmax = 1000;
//    double _x = _x0;  //�ʱⰪ
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
//        hk = -(myFunc(_x) / mydFunc(_x));
//        xn = _x + hk;
//
//        if (myFunc(xn) == 0)
//        {
//            break;
//        }
//        ep = fabs(hk);
//
//        _x = xn;
//        k++;
//
//    } while (k<Nmax && fabs(ep)>_tol);
//
//    return xn;
//}