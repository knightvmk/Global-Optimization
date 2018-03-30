#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <stdio.h>
#include <vector>
#include <queue>
#include <list>

#include <omp.h>

#include <math.h>

struct block
{
	double left;
	double right;
};


class global
{
private:
	double left_x; // левая граница
	double right_x; // правая граница
	double r; // параметр метода

	std::list<double> Values_Z;
	std::list<double> Arguments_X;

	std::vector<std::list<double>> parallel_z;
	std::vector<std::list<double>> parallel_x;

	double solutionX = 0;
	double solutionZ = 0;
	int procs = 1;
	long double time = 0;
	__int64 steps = 0;
	__int64 parallel_steps = 0;

public:
	global(const double _left, const double _right, const double _r) { left_x = _left, right_x = _right, r = _r; }; // использовать данный конструктор
	global() {left_x = 0, right_x = 0, r = 0,time = 0.L;}; // по-умолч
	double Optimize(const double _left, const double _right, const int _N_max, const double _Eps); // последовательный АГП с возвратом x_min
	double Optimize(const double _left, const double _right, const int _N_max, const double _Eps, const int _threads); // параллельный по отрезку АГП с возвратом x_min
	inline double R(const double &_m_small, const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev);
	inline double M(const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev);
	inline double func(const double &_x);
	long double Time();
	double GetSolutionX() { return solutionX; };
	double GetSolutionZ() { return solutionZ; };
	int GetSteps() { return steps; };
	void Set(const double _left, const double _right, const double _r) { left_x = _right; right_x = _right; r = _r; };
};

#endif
