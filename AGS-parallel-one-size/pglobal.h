
#ifndef __PGLOBAL_H__
#define __PGLOBAL_H__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <queue>
#include <list>

#include <omp.h>

struct block
{
	double x_left; // координата слева
	double x_right; // координата справа
	double z_left; // значение слева
	double z_right; // значение справа
	double R; // оценка константы Липшица
};

bool operator<(const block &i1, const block &i2)
{
	return (i1.R < i2.R) ? true : false;
}


class parallel_global
{
private:
	double left; // левая граница
	double right; // правая граница
	double r; // параметр
	std::priority_queue < block, std::vector<block>, decltype(&operator<)> arr;
	double solutionX = 0;
	double solutionZ = 0;
	int procs = 1;
	long double time;
	__int64 steps;


public:
	std::list<block> list_arr_x;
	std::list<block> list_arr_z;

	parallel_global();
	parallel_global(const double _left, const double _right, const double _r) { left = _left; right = _right; r = _r; }
	parallel_global(const parallel_global &_A) { left = _A.left; right = _A.right; arr = _A.arr; }
	~parallel_global() { while (!arr.empty()) arr.top(); }

	void Set(const double _left, const double _right, const double _r) { left = _left; right = _right; r = _r; }
	int Optimize(const double _Epsilon, const int _Steps);
	inline double R(const double _m_small, const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev);
	inline double M(const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev);
	inline double Func(const double _x);
	long double Time() { return time; }
	double GetSolutionX() { return solutionX; };
	double GetSolutionZ() { return solutionZ; };
	int GetSteps() { return steps; };


};


#endif

