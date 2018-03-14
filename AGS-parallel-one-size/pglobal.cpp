#include "pglobal.h"

int Parallel_global::Optimize(const double _Epsilon, const int _Steps)
{
	time = 0.L;
	int procs = omp_get_max_threads();
	int k = 0; // количество шагов
	double eps = _Epsilon;
	double curr_left, curr_right; // текущее значение слева и справа для оценки разности

	curr_left = left; // инициализация границ
	curr_right = right; // инициализация границ
	block *work = new block();
	work->x_left = curr_left;
	work->x_right = curr_right;
	work->z_left = Func(curr_left);
	work->z_right = Func(curr_right);
	arr.push(*work);
	block *parallel_block = new block[procs];
	while (k < _Steps && (fabs(curr_left - curr_right) > eps))
	{

		for (int i = 0; i < procs; i++)
		{

		}
	}


	return 0;
}

inline double Parallel_global::R(const double _m_small, const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev)
{
	return _m_small*(_x_curr - _x_prev) + pow(_z_curr - _z_prev, 2) / (_m_small*(_x_curr - _x_prev)) - 2 * (_z_curr + _z_prev);
}

inline double Parallel_global::M(const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev)
{
	return fabs((_z_curr - _z_prev) / (_x_curr - _x_prev));
}

inline double Parallel_global::Func(const double _x)
{
	return sin(_x) + sin((10. * _x) / 3.);
	//return (2 * powl((_x - 3), 2) + exp((powl(_x, 2) / 2)));
	//return ((3 * _x - 1.4)*sin(18 * _x));
	//return (sin(_x) + sin((10. * _x) / 3.) + log(_x) - 0.84*_x + 3);
	//return ((pow(_x, 2) - 5 * _x + 6) / (pow(_x, 2) + 1));
}