///*-------------------------------------------------------------------------------\
//@ Numerical Methods by Young-Keun Kim - Handong Global University
//
//Author           : [YOUR NAME]
//Created          : 26-03-2018
//Modified         : 03-04-2021
//Language/ver     : C++ in MSVS2019
//
//Description      : NM_MainTemplate.cpp
//-------------------------------------------------------------------------------*/
//
//#define Assignment	4		// enter your assignment number
//#define eval		0		// set 0
//
//#include "myNM.h"
//double L1_Norm(Matrix _A);
//double L2_Norm(Matrix _A);
//double inf_Norm(Matrix _A);
//double vec_Norm(Matrix _v);
//void QR_householder(Matrix _A);
//Matrix eig(Matrix _A);
//double cond(Matrix _A);
//
//int main(int argc, char* argv[])
//{
//	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
//	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";
//
//#if eval
//	path += "eval/";
//#endif
//
//	/*==========================================================================*/
//	/*					Variables declaration & initialization					*/
//	/*--------------------------------------------------------------------------*/
//	/*   - You can change the variable names									*/
//	/*   - However, you must use the specified file name						*/
//	/*	   : For each assignment, the file name will be notified on HISNET		*/
//	/*==========================================================================*/
//	Matrix matA = txt2Mat(path, "prob1_matA");
//	Matrix matC = txt2Mat(path, "prob1_matC");
//	//Matrix vecb = txt2Mat(path, "prob1_vecb");
//	//Matrix matd = txt2Mat(path, "prob1_matA");
//	//Matrix vecd = txt2Mat(path, "prob1_vecd");
//	//Matrix vecx_true = txt2Mat(path, "prob1_vecx_true");
//	
//	
//	/*==========================================================================*/
//	/*							  Print your results							*/
//	/*==========================================================================*/
//	
//	printf("------------------------------------------------------------------------------------\n");
//	//printMat(matA, "A");
//	printMat(matC, "C");
//
//	/*QR_householder(matA);
//	eig(matA);*/
//	//cond(matA);
//	cond(matC);
//	printf("------------------------------------------------------------------------------------\n");
//	
//
//	/*==========================================================================*/
//	/*							  Deallocate memory 							*/
//	/*==========================================================================*/
//	freeMat(matA);		/*freeMat(vecb);*/
//
//	//freeMat(matU);		
//	freeMat(matC);		//freeMat(vecx_true);
//	//freeMat(matP);
//	//freeMat(matL);
//	system("pause");
//	return 0;
//}