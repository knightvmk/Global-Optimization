#ifndef __TGLOBAL_H__
#define __TGLOBAL_H__

#include "TGlobal_defines.h"

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <list>

#include <omp.h>
#include <amp.h>
#include <amp_math.h>
//#include <mpi.h>
#include "evolvent.h"
#include <iostream>
//#include <CL\opencl.h>
#include "grishagin\include\grishagin_function.hpp"
//#include "math_parcer\exprtk.hpp"
#include <memory>
#include "math_parcer\MParcer.h"

struct block_1D
{
	double R; // оценка константы Липшица
	double x_left; // координата слева
	double x_right; // координата справа
	double z_left; // значение слева
	double z_right; // значение справа
	double M; // параметр оценки для интервала
};

struct block_2D
{
	double R; // оценка константы Липшица
	double x_img_left;
	double x_img_right;
	double x_left; // координата слева
	double x_right; // координата справа
	double y_left; // координата слева
	double y_right; // координата справа
	double z_left; // значение слева
	double z_right; // значение справа
	double M; // параметр оценки для интервала
};

//struct math_expression
//{
//	exprtk::symbol_table<double> symbol_table;
//	exprtk::expression<double> expression;
//	exprtk::parser<double> parser;
//};

void Quick(concurrency::array_view<block_1D> &s_arr, int first, int last) restrict(amp, cpu);
void Quick(concurrency::array_view<block_2D> &s_arr, int first, int last) restrict(amp, cpu);

void Bubble(concurrency::array_view<block_1D> &arr, int size) restrict(amp, cpu);
void Bubble(concurrency::array_view<block_2D> &arr, int size) restrict(amp, cpu);

bool operator<(const block_1D &i1, const block_1D &i2) restrict(amp, cpu);
bool operator>(const block_1D &i1, const block_1D &i2) restrict(amp, cpu);

bool operator<(const block_2D &i1, const block_2D &i2) restrict(amp, cpu);
bool operator>(const block_2D &i1, const block_2D &i2) restrict(amp, cpu);

class TGlobal
{
protected:
#ifndef NOT_USE_CL
	std::vector<cl::Platform> avaible_platforms;
	std::vector<std::vector<cl::Device>> avaible_devices;
#endif
	std::vector<std::string> avaible_platform_names;
	std::vector<std::string> avaible_devices_names;
	std::vector<concurrency::accelerator> accelerators;
	bool can_use_cl = true;
	void InitCL();
	//math_expression *matexp; // структура для параллельной версии
	TCALC *matexp;
	//std::string expression_string; // мат.выражение
	char *math_exp_char = new char;;
public:
	TGlobal() 
	{
#ifdef NOT_USE_CL
		can_use_cl = false;
#else
		can_use_cl = true;
#endif
		accelerators = concurrency::accelerator::get_all();
		//delete matexp;
		matexp = new TCALC[omp_get_max_threads()]();
		//InitCL();
	}
	void SetExpression(std::string str_exp) 
	{ 
		//expression_string = str_exp; 
		delete math_exp_char;
		math_exp_char = new char[str_exp.length() + 1];
		strcpy(math_exp_char, str_exp.c_str());
	};
	~TGlobal() 
	{
#ifndef NOT_USE_CL
		avaible_platforms.clear(); avaible_platforms.~vector();
		avaible_devices.clear(); avaible_devices.~vector();
#endif
		avaible_platform_names.clear(); avaible_platform_names.~vector();
		avaible_devices_names.clear(); avaible_devices_names.~vector();
		delete[]matexp;
		delete[]math_exp_char;
		//symbol_table.clear(); symbol_table.~symbol_table();
		//expression.release(); expression.~expression();
		//parser.~parser();
		//expression_string.clear(); expression_string.~basic_string();
	}
};

class OneDimension : public TGlobal
{

protected:
	double left; // левая граница
	double right; // правая граница
	double r; // параметр
	double solutionX = 0;
	double solutionZ = 0;
	int procs = 1; // процессоров
	long double time; // время работы
	int steps; // кол-во шагов

public:
	OneDimension(const double _left, const double _right, const double _r) : left(_left), right(_right), r(_r) { TGlobal(); }
	virtual ~OneDimension() {};
	virtual long double Time() { return time; };
	virtual double GetSolutionX() { return solutionX; };
	virtual double GetSolutionZ() { return solutionZ; };
	virtual int GetSteps() { return steps; };
	virtual void Set(const double _left, const double _right, const double _r) { left = _right; right = _right; r = _r; };

protected:
	inline double R(const double &_m_small, const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev);
	inline double M(const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev);
	inline double Func(const double &_x, int &thread);
	//inline double Calc(double _x) 
	//{
	//	int thread = omp_get_thread_num();
	//	double x_copy = _x;
	//	matexp[thread].symbol_table.add_variable("x", x_copy);
	//	matexp[thread].symbol_table.add_constants();
	//	matexp[thread].expression.register_symbol_table(matexp[thread].symbol_table);
	//	matexp[thread].parser.compile(expression_string, matexp[thread].expression);
	//	return matexp[thread].expression.value();
	//};
};

class Sequental_1D : public OneDimension
{
private:

public:
	std::list<double> Values_Z;
	std::list<double> Arguments_X;

	Sequental_1D(const double _left, const double _right, const double _r) : OneDimension(_left, _right, _r) {}
	~Sequental_1D() { Values_Z.clear(); Arguments_X.clear(); }

	double Optimize(const double _left, const double _right, const int _N_max, const double _Eps);
};

class Parallel_1D : public OneDimension
{
private:
	int parallel_steps = 0;
	std::vector<std::list<double>> parallel_z;
	std::vector<std::list<double>> parallel_x;
	double Optimize(const double _left, const double _right, const int _N_max, const double _Eps);
public:
	std::list<double> Values_Z;
	std::list<double> Arguments_X;
	Parallel_1D(const double _left, const double _right, const double _r) : OneDimension(_left, _right, _r) {}
	~Parallel_1D() { parallel_z.clear(); parallel_x.clear(); Values_Z.clear(); Arguments_X.clear(); }

	double Optimize(const double _left, const double _right, const int _N_max, const double _Eps, const int _threads); //оптимизация делений по отрезку
	double Optimize(const double _Epsilon, const int _Steps, const int _thread_count); // оптимизация по характеристикам R отрезков с R_max приоритетом
	double OptimizeGPU(const double _left, const double _right, const int _N_max, const double _Eps); // оптимизация на с привлечением GPU на C++ AMP
};

class TwoDimension : public TGlobal
{
protected:
	double* Left; // ограничение области слева
	double* Right; // ограничение области справа
	double r; // параметр
	double solutionX = 0;
	double solutionY = 0;
	double solutionZ = 0;
	int procs = 1;
	long double time;
	int steps;
public:
	virtual long double Time() { return time; };
	virtual double GetSolutionX() { return solutionX; };
	virtual double GetSolutionY() { return solutionY; };
	virtual double GetSolutionZ() { return solutionZ; };
	virtual int GetSteps() { return steps; };
	virtual void Set(const double *_Left, const double *_Right, const double _r) 
	{ 
		delete[] Left;
		delete[] Right;
		Left = new double[2];
		Right = new double[2];
		for (int i = 0; i < 2; i++)
		{
			Left[i] = _Left[i];
			Right[i] = _Right[i];
		}
	};
protected:
	inline double R(const double &_m_small, const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev);
	inline double M(const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev);
	inline double Func(const double *_y, int &thread);
	//inline double Calc(const double *_y)
	//{
	//	double y_arr[2] = { _y[0], _y[1] };
	//	matexp[thread].symbol_table.add_variable("x", y_arr[0]);
	//	matexp[thread].symbol_table.add_variable("y", y_arr[1]);
	//	matexp[thread].symbol_table.add_constants();
	//	matexp[thread].expression.register_symbol_table(matexp[thread].symbol_table);
	//	matexp[thread].parser.compile(expression_string, matexp[thread].expression);
	//	return matexp[thread].expression.value();
	//};
	inline double Test(const double *_y, vagrish::GrishaginFunction *func, int num);
	
private:

public:
	TwoDimension(const double *_Left, const double *_Right, const double _r) : r(_r)
	{
		TGlobal();
		Left = new double[2];
		Right = new double[2];
		for (int i = 0; i < 2; i++)
		{
			Left[i] = _Left[i];
			Right[i] = _Right[i];
		}

	}
	~TwoDimension() { delete[]Left; delete[]Right; }
};

class Sequental_2D : public TwoDimension
{
public:
	std::list<double> Values_Z;
	std::list<double> Arguments_Y;
	std::list<double> Arguments_X;
	Sequental_2D(const double *_Left, const double *_Right, const double _r) : TwoDimension(_Left, _Right, _r) {};
	~Sequental_2D() { Values_Z.clear(); Arguments_X.clear(); Arguments_Y.clear(); }

	double Optimize(const double *_Left, const double *_Right, const int _N_max, const double _Eps);
	int RunTest();
private:
	double Optimize(const double *_Left, const double *_Right, const int _N_max, const double _Eps, int num, vagrish::GrishaginFunction *func);
};

class Parallel_2D : public TwoDimension
{
private:
	int parallel_steps = 0;
public:
	std::list<double> Values_Z;
	std::list<double> Arguments_Y;
	std::list<double> Arguments_X;
	Parallel_2D(const double *_Left, const double *_Right, const double _r) : TwoDimension(_Left, _Right, _r) {};
	~Parallel_2D() { Values_Z.clear(); Arguments_X.clear(); Arguments_Y.clear(); }

	double Optimize(const double _Epsilon, const int _Steps, const int _thread_count); // оптимизация по характеристикам R отрезков с R_max приоритетом
	double OptimizeGPU(const double _Epsilon, const int _Steps); // оптимизация на с привлечением GPU на C++ AMP
	int RunTest(int _threads);
	int RunTestGPU();
private:
	double Optimize(const double _Epsilon, const int _Steps, const int _thread_count, int num, vagrish::GrishaginFunction *func); // оптимизация по характеристикам R отрезков с R_max приоритетом
	double OptimizeGPU(const double _Epsilon, const int _Steps, int num, vagrish::GrishaginFunction *func); // оптимизация на с привлечением GPU на C++ AMP
};

// no realization

class NDimension : public TGlobal
{
protected:
	//
private:
	//
public:
	// 
	NDimension() { throw("NO CODE"); }
};

class Sequental_ND : public NDimension
{
	//
	Sequental_ND() { throw("NO CODE"); }
};

class Parallel_ND : public NDimension
{
	//
	Parallel_ND() { throw("NO CODE"); }
};

#endif