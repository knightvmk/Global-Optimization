#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <omp.h>
#include "TGlobal.h"

int main(int argv, char* argc[])
{
	//system("color 0B"); // черный-голубой
	//system("color 0A"); // черный-зеленый
	//system("color 02"); // черно-зеленый
	system("color 1B"); // сине-голубой
	//system("color F1"); // бело-черный
	system("title Optimizer - AGS");
	double *A = new double[2];
	double *B = new double[2];
	A[0] = -4; A[1] = -5;
	B[0] = 2; B[1] = 0;
	double eps = 0.00001;
	int n_max = 1000;
	//Parallel_2D GrishaginTest0(A, B, 2);
	//int fails0 = GrishaginTest0.RunTest(omp_get_max_threads());
	//printf("\n Failed tests = [ %i ] \n", fails0);
	////system("pause");
	//Parallel_2D GrishaginTest1(A, B, 2);
	//int fails1 = GrishaginTest1.RunTestGPU();
	//printf("\n Failed tests = [ %i ] \n", fails1);
	////system("pause");
	//Sequental_2D GrishaginTest2(A, B, 2);
	//int fails2 = GrishaginTest2.RunTest();
	//printf("\n Failed tests = [ %i ] \n", fails2);
	//system("pause");

	std::string exp2D = "2 * cos(Y) + 3 * sin(X)";

	Parallel_2D test03(A, B, 2);
	test03.SetExpression(exp2D);
	test03.OptimizeGPU(eps, n_max);
	printf("\n Parallel 2-Dimension optimization by evolving with GPGPU\n");
	printf("\n Min X = [ %.15lf ]", test03.GetSolutionX());
	printf("\n Min Y = [ %.15lf ]", test03.GetSolutionY());
	printf("\n Min Z = [ %.15lf ]     [ N_max = %i , Eps = %.8lf ]\n Step count = [ %i ]\n", test03.GetSolutionZ(), n_max, eps, test03.GetSteps());
	printf("\n Time = [ %lf ] sec\n", test03.Time());
	printf("\n * * * * * * * * * * * * * * * * * * \n");

	Parallel_2D test01(A, B, 2);
	test01.SetExpression(exp2D);
	test01.Optimize(eps, n_max, omp_get_max_threads());
	printf("\n Parallel 2-Dimension optimization by evolving\n");
	printf("\n Threads = [ %i ]", omp_get_max_threads());
	printf("\n Min X = [ %.15lf ]", test01.GetSolutionX());
	printf("\n Min Y = [ %.15lf ]", test01.GetSolutionY());
	printf("\n Min Z = [ %.15lf ]     [ N_max = %i , Eps = %.8lf ]\n Step count = [ %i ]\n", test01.GetSolutionZ(), n_max, eps, test01.GetSteps());
	printf("\n Time = [ %lf ] sec\n", test01.Time());
	printf("\n * * * * * * * * * * * * * * * * * * \n");

	Sequental_2D test0(A, B, 2);
	test0.SetExpression(exp2D);
	//eps = 0.00001;
	//n_max = 1200;
	test0.Optimize(A, B, n_max, eps);
	printf("\n Sequental 2-Dimension optimization by evolving\n");
	printf("\n Min X = [ %.15lf ]", test0.GetSolutionX());
	printf("\n Min Y = [ %.15lf ]", test0.GetSolutionY());
	printf("\n Min Z = [ %.15lf ]     [ N_max = %i , Eps = %.8lf ]\n Step count = [ %i ]\n", test0.GetSolutionZ(), n_max, eps, test0.GetSteps());
	printf("\n Time = [ %lf ] sec\n", test0.Time());
	printf("\n * * * * * * * * * * * * * * * * * * \n");
	//system("pause");

	std::string exp1D = "2*sin(3*X)+3*cos(5*X)";

	Parallel_1D test12(0, 8, 2);
	test12.SetExpression(exp1D);
	eps = 0.000001;
	n_max = 1000;
	test12.OptimizeGPU(4, 8, n_max, eps);
	printf("\n Parallel optimization on Accelerators with GPGPU\n");
	printf("\n Min X = [ %.15lf ]", test12.GetSolutionX());
	printf("\n Min Z = [ %.15lf ]     [ N_max = %i , Eps = %.8lf ]\n Step count = [ %i ]\n", test12.GetSolutionZ(), n_max, eps, test12.GetSteps());
	printf("\n Time = [ %lf ] sec\n", test12.Time());
	printf("\n * * * * * * * * * * * * * * * * * * \n");

	Parallel_1D test1(0, 8, 2);
	test1.SetExpression(exp1D);
	eps = 0.000001;
	n_max = 1000;
	test1.Optimize(eps, n_max, omp_get_max_threads());
	printf("\n Parallel optimization by R_max intervals\n");
	printf("\n Min X = [ %.15lf ]", test1.GetSolutionX());
	printf("\n Min Z = [ %.15lf ]     [ N_max = %i , Eps = %.8lf ]\n Step count = [ %i ]\n", test1.GetSolutionZ(), n_max, eps, test1.GetSteps());
	printf("\n Threads = [ %i ]", omp_get_max_threads());
	printf("\n Time = [ %lf ] sec\n", test1.Time());
	printf("\n * * * * * * * * * * * * * * * * * * \n");

	Parallel_1D test2(0, 8, 2);
	test2.SetExpression(exp1D);
	eps = 0.000001;
	n_max = 1000;
	test2.Optimize(0, 8, n_max, eps, omp_get_max_threads());
	printf("\n Parallel optimization by division to intervals\n");
	printf("\n Min X = [ %.15lf ]", test2.GetSolutionX());
	printf("\n Min Z = [ %.15lf ]     [ N_max = %i , Eps = %.8lf ]\n Step count = [ %i ]\n", test2.GetSolutionZ(), n_max, eps, test2.GetSteps());
	printf("\n Threads = [ %i ]", omp_get_max_threads());
	printf("\n Time = [ %lf ] sec\n", test2.Time());
	printf("\n * * * * * * * * * * * * * * * * * * \n");

	Sequental_1D test4(0, 8, 2);
	test4.SetExpression(exp1D);
	eps = 0.000001;
	n_max = 1000;
	test4.Optimize(0, 8, n_max, eps);
	printf("\n Sequental optimization\n");
	printf("\n Min X = [ %.15lf ]", test4.GetSolutionX());
	printf("\n Min Z = [ %.15lf ]     [ N_max = %i , Eps = %.8lf ]\n Step count = [ %i ]\n", test4.GetSolutionZ(), n_max, eps, test4.GetSteps());
	printf("\n Time = [ %lf ] sec\n", test4.Time());
	printf("\n * * * * * * * * * * * * * * * * * * \n");

	system("pause");
	return 0;
}