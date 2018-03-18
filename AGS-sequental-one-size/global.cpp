#include "global.h"

global::global()
{
	left_x = 0;
	right_x = 0;
	r = 0;
	time = 0.L;
}

double global::solve(const double _left, const double _right, const double _r, const int _N_max, const double _Eps)
{
	time = 0.L;

	procs = omp_get_num_procs();
	left_x = _left;
	double curr_left = _left;
	right_x = _right;
	double curr_right = _right;
	r = _r;
	// очистка списка
	Values_Z.clear();
	Arguments_X.clear();
	
	double start, end;

	/////
	double z_begin = func(left_x); // левая графица
	double z_end = func(right_x); // правая граница
	//// прочие переменные
	double temp = 0.; // переменная для подсчёта 
	int k = 1; // количество шагов в Настоящем алгоритме
	int new_size; // тоже шаги
	double M_big, m_small = 1.; // Перменные для хранения итоговой константы "М" и константы "m"
	int R_max_index = 0; // Переменная индекса следующей точки испытания

	double left, right, new_point; // переменные для рабочего отрезка и точки нового испытания
	////
	//// контейнеры для хранения значений и аргументов функции
	std::list<double> z; //значения Z
	std::list<double> x; //значения рабочих X
	std::list<double> M_vector; // вектор неких коэфф M для выбора m = max(M[i])
	std::list<double> R_vector; // вектор значений R для принятия решения точки следующего испытания
	//std::queue<double> M_queue; 
	//std::queue<double> R_queue;


	//// блок итераторов для работы со списками
	std::list<double>::iterator place_x, place_z, place_M, place_R; // итераторы-указатели КУДА положить новые значения
	std::list<double>::iterator iter_x, iter_z, iter_R, iter_M;
	////
	left = curr_left;
	right = curr_right;
	new_point = temp = (curr_left + curr_right) / 2;

	x.push_back(left);
	x.push_back(new_point);
	x.push_back(right);

	iter_x = x.begin();
	for (int i = 0; i < 3; i++)
	{
		z.push_back(func(*iter_x++));
	}

	iter_z = z.begin();
	iter_x = x.begin();

	//// начало алгоритма // первый шаг, новая точка - средняя //

	block val_z; // структура для хранения левой и правой границы отрезка с Z
	block arg_x; // структура для хранения левой и правой границы отрезка с Х

	//Начало времени отсчёта работы алгоритма
	start = omp_get_wtime();

#pragma omp single
	{
		for (int i = 1; i < 3; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;

			M_vector.push_back(M(val_z.right, val_z.left, arg_x.right, arg_x.left));
		}
		double max = *M_vector.begin();
		if (max < *M_vector.begin() + 1) max = *M_vector.begin() + 1;
		M_big = max;
		if (M_big > 0.) m_small = r*M_big;
		else if (M_big == 0.) m_small = 1.;

		R_max_index = 0;
		max = 0.;

		iter_z = z.begin();
		iter_x = x.begin();

		for (int i = 1; i < 3; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;

			temp = R(m_small, arg_x.right, arg_x.left, val_z.right, val_z.left);
			R_vector.push_back(temp);

			if (max < temp)
			{
				max = temp;
				R_max_index = i;
			}
		}

		iter_z = z.begin();
		iter_x = x.begin();

		for (int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}

		solutionZ = *iter_z;
		solutionX = *iter_x;
		curr_left = *--iter_x;
		curr_right = *++iter_x;

		iter_z = z.begin();
		iter_x = x.begin();
	}

	//// начало работы АГП с испытанием в новой точке по правилу //
	while (k < _N_max && (fabs(curr_left - curr_right) > _Eps))
	{
		new_size = 2 + k;

		place_z = z.begin();
		place_x = x.begin();
		for (int i = 0; i < R_max_index; i++)
		{
			place_z++;
			place_x++;
		}
		val_z.right = *place_z--;
		val_z.left = *place_z++;

		left = curr_left;
		right = curr_right;
		new_point = 0.5*(right + left) - (val_z.right - val_z.left) / (2 * m_small);

		z.insert(place_z, func(new_point));
		x.insert(place_x, new_point);

		iter_z = z.begin();
		iter_x = x.begin();
		place_M = M_vector.begin();
		place_R = R_vector.begin();
		for (int i = 0; i < R_max_index - 1; i++)
		{
			iter_z++;
			iter_x++;
			if (i == R_max_index - 1) break;
			place_M++;
			place_R++;
		}

		for (int i = 0; i < 2; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;
			if (i == 0)
			{
				M_vector.insert(place_M, M(val_z.right, val_z.left, arg_x.right, arg_x.left));
			}
			else
			{
				*place_M = M(val_z.right, val_z.left, arg_x.right, arg_x.left);
			}
		}

		iter_M = M_vector.begin();
		double max = *iter_M++;
		temp = *iter_M;
		for (int i = 0; i < new_size - 1; i++)
		{
			if (max < temp) max = temp;
		}
		M_big = max;
		if (M_big > 0.) m_small = r*M_big;
		else if (M_big == 0.) m_small = 1.;
		R_max_index = 1;
		max = 0.;

		--iter_z; --iter_z;
		--iter_x; --iter_x;

		for (int i = 0; i < 2; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;
			if (i == 0)	R_vector.insert(place_R, R(m_small, arg_x.right, arg_x.left, val_z.right, val_z.left));
			else
				*place_R = R(m_small, arg_x.right, arg_x.left, val_z.right, val_z.left);
		}

		iter_R = R_vector.begin();
		max = *iter_R;
		for (int i = 0; i < new_size - 1; i++)
		{
			temp = *iter_R++;
			if (max < temp)
			{
				max = temp;
				R_max_index = i + 1;
			}
		}

		iter_x = x.begin();
		iter_z = z.begin();
		for (int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}
		solutionZ = *iter_z;
		solutionX = *iter_x;
		curr_left = *--iter_x;
		curr_right = *++iter_x;
		++k;
	}

	end = omp_get_wtime();
	time = end - start;

	--new_size;
	steps = new_size;
	iter_x = x.begin();
	iter_z = z.begin();
	for (register int i = 0; i < new_size; i++)
	{
		Arguments_X.push_back(*iter_x++);
		Values_Z.push_back(*iter_z++);
	}
	return solutionZ;
}

inline double global::R(const double _m_small, const double _x_curr, const double _x_prev, const double _z_curr, const double _z_prev)
{
	return _m_small*(_x_curr - _x_prev) + pow(_z_curr - _z_prev, 2) / (_m_small*(_x_curr - _x_prev)) - 2 * (_z_curr + _z_prev);
}

inline double global::M(const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev)
{
	return fabs((_z_curr - _z_prev) / (_x_curr - _x_prev));
}

inline double global::func(const double _x)
{
	return sin(_x) + sin((10. * _x) / 3.);
	//return (2 * powl((_x - 3), 2) + exp((powl(_x, 2) / 2)));
	//return ((3 * _x - 1.4)*sin(18 * _x));
	//return (sin(_x) + sin((10. * _x) / 3.) + log(_x) - 0.84*_x + 3);
	//return ((pow(_x, 2) - 5 * _x + 6) / (pow(_x, 2) + 1));
}

long double global::Time()
{
	return time;
}
