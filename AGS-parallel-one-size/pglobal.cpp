#include "pglobal.h"

bool operator<(const block& i1, const block& i2)
{
	return (i1.R < i2.R) ? true : false;
}

bool operator>(const block& i1, const block& i2)
{
	return (i1.R > i2.R) ? true : false;
}


double parallel_global::Optimize(const double _Epsilon, const int _Steps)
{
	time = 0.L;
	int procs = omp_get_max_threads();
	int k = 0; // количество шагов
	double eps = _Epsilon;
	double curr_left, curr_right; // текущее значение слева и справа для оценки разности
	double start, end; // время работы
	double m_small = 1., M_big;
	double *x_new = new double[procs];

	curr_left = left; // инициализация границ
	curr_right = right; // инициализация границ

	std::vector<block> InitVector(procs);
	std::vector<block> WorkVector(procs * 2);
	std::priority_queue<block> Queue;
	std::list<block> Intervals;

	int size_bank_intervals; // размер list<block> Intervals

	std::list<block>::iterator Pointer;


	// начало первичной инициализации параллельного алгоритма
	// деление отрезка на звенья по количеству ЦП
	// вычисление R характеристик для звеньев - интервалов
	double len = right - left; // расстояние поиска
	double li = len / procs; // длина i-го участка
	double tmp = 0.;

	omp_set_dynamic(0);

	start = omp_get_wtime();
	
	omp_set_num_threads(procs);
	{
#pragma omp parallel for firstprivate(tmp)
		for (int i = 0; i < procs; i++)
		{
			tmp = (right - len) + (i * li); // _b-(b-a)_ + _расстояние_участка_ * _номер_участка_;
			InitVector[i].x_left = tmp;
			InitVector[i].x_right = tmp + li;
			InitVector[i].z_left = Func(tmp);
			InitVector[i].z_right = Func(tmp + li);
		}

		// считаем величины М каждого интервала
#pragma omp parallel for
		for (int i = 0; i < procs; i++)
		{
			InitVector[i].M = M(InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_right, InitVector[i].x_left);
		}

		M_big = InitVector[0].M;

		for (register int i = 1; i < procs; i++) // это надо б как-нить распараллелить (сомнительная выгода) // нахождение макса в массиве длины CPU_COUNT
		{
			if (M_big < InitVector[i].M) M_big = InitVector[i].M;
		}
		m_small = (M_big == 0) ? 1 : r * M_big;

#pragma omp parallel for 
		for (int i = 0; i < procs; i++)
		{
			InitVector[i].R = R(m_small, InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_right, InitVector[i].x_left);
		}
		for (register int i = 0; i < procs; i++)
		{
			Intervals.push_back(InitVector[i]);
		}
		Intervals.sort();
		
		curr_left = Intervals.back().x_left;
		curr_right = Intervals.back().x_right;

		++k;
		// начало параллельного строгого алгоритма
		///////////////////////////////////////////
		//////////// Параллельный блок ////////////
		///////////////////////////////////////////
		while (k < _Steps && (abs(curr_left - curr_right) > eps))
		{
			Pointer = --Intervals.end();
			for (register int i = 0; i < procs * 2-2; i += 2) // забрали на работу CPU_COUNT участков с макс. хар-кой R в алгоритм из очереди
			{
				WorkVector[i] = *Pointer;
				Intervals.erase(Pointer--);
			}
			WorkVector[procs * 2 - 2] = *Pointer; Intervals.erase(Pointer);
			size_bank_intervals = Intervals.size();
			// ищем на каждом участке точку нового испытания x_new
#pragma omp parallel for shared(x_new)
			for (int i = 0; i < procs; i++)
			{
				x_new[i] = 0.5*(WorkVector[i * 2].x_right + WorkVector[i * 2].x_left)
					- (WorkVector[i * 2].z_right - WorkVector[i * 2].z_left) / (2 * m_small);
				WorkVector[i * 2 + 1].z_left = Func(x_new[i]);
			}
			for (register int i = 0; i < procs * 2; i += 2)
			{
				WorkVector[i + 1].x_right = WorkVector[i].x_right;
				WorkVector[i + 1].z_right = WorkVector[i].z_right;
				WorkVector[i + 1].x_left = WorkVector[i].x_right = x_new[i / 2];
				WorkVector[i].z_right = WorkVector[i + 1].z_left;
			}

			// всё скопировали, ничего не потеряли, пора считать характеристики
#pragma omp parallel for
			for (int i = 0; i < procs * 2; i++)
			{
				WorkVector[i].M  = M(WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_right, WorkVector[i].x_left);
			}

			double M_max_array = WorkVector[0].M;
#pragma omp parallel sections num_threads(2)
			{

#pragma omp section
				{
					for (register int i = 1; i < procs * 2; i++) // это надо б как-нить распараллелить (сомнительная выгода) // нахождение макса в массиве длины CPU_COUNT * 2
					{
						if (M_max_array < WorkVector[i].M) M_max_array = WorkVector[i].M;
					}
				}
#pragma omp section
				{
					// теперь ищем M_max среди отрезоков, лежащих в приоритетной очереди
					if (size_bank_intervals != 0)
					{
						Pointer = Intervals.begin();
						double tmp_M; //= (*Pointer).M;
						double M_max_list = (*Pointer++).M;

						for (; Pointer != Intervals.end(); Pointer++)
						{
							tmp_M = (*Pointer).M;
							if (M_max_list < tmp_M) M_max_list = tmp_M;
						}

						// сравниваем M_max_array и M_max_list, выбираем большее
						M_big = (M_max_array < M_max_list) ? M_max_list : M_max_array;
					}
					else M_big = M_max_array;
				}
			}
			m_small = (M_big == 0) ? 1 : r * M_big;

#pragma omp parallel for 
			for (int i = 0; i < procs * 2; i++)
			{
				WorkVector[i].R = R(m_small, WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_right, WorkVector[i].x_left);
			}
			// рассчитали R характеристики для рабочих отрезков
			// необходимо пересчитать R для отрезков в очереди
			if (size_bank_intervals != 0)
			{
				std::vector<block> R_vec(Intervals.begin(), Intervals.end());
#pragma omp parallel for
				for (register int i = 0; i < size_bank_intervals; i++)
				{
					R_vec[i].R = R(m_small, R_vec[i].z_right, R_vec[i].z_left, R_vec[i].x_right, R_vec[i].x_left);
				}
				register int j = 0;
				for (auto i = Intervals.begin(); i != Intervals.end(); i++)
				{
					i->R = R_vec[j++].R;
				}
				R_vec.clear();
			}
			for (register int i = 0; i < procs * 2; i++)
			{
			 // возвращаем CPU_COUNT отрезков с подсчитанными характеристиками в приоритетную очередь для строгого алгоритма
				Intervals.push_back(WorkVector[i]);
			}
			Intervals.sort();
			++k;
			curr_left = Intervals.back().x_left;
			curr_right = Intervals.back().x_right;
		}
	}

	end = omp_get_wtime();
	time = end - start;

	omp_set_dynamic(1);

	solutionX = curr_right;
	solutionZ = Intervals.back().z_right;

	steps = k;
	// очистка
	InitVector.clear();
	WorkVector.clear();
	delete[] x_new;

	return solutionZ;
}

double parallel_global::Optimize(const double _Epsilon, const int _Steps, const int _thread_count)
{
	time = 0.L;
	int procs = omp_get_max_threads();
	int k = 0; // количество шагов
	double eps = _Epsilon;
	double curr_left, curr_right; // текущее значение слева и справа для оценки разности
	double start, end; // время работы
	double m_small = 1., M_big;
	double *x_new = new double[procs];


	curr_left = left; // инициализация границ
	curr_right = right; // инициализация границ

	std::vector<block> InitVector(procs);
	std::vector<block> WorkVector(procs * 2);
	std::priority_queue<block> Queue;
	std::list<block> Intervals;

	int size_bank_intervals; // размер list<block> Intervals

	std::list<block>::iterator Pointer;


	// начало первичной инициализации параллельного алгоритма
	// деление отрезка на звенья по количеству ЦП
	// вычисление R характеристик для звеньев - интервалов
	double len = right - left; // расстояние поиска
	double li = len / procs; // длина i-го участка
	double tmp = 0.;

	omp_set_dynamic(0);

	start = omp_get_wtime();

	omp_set_num_threads(procs);
	{
#pragma omp parallel for firstprivate(tmp) num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			tmp = (right - len) + (i * li); // _b-(b-a)_ + _расстояние_участка_ * _номер_участка_;
			InitVector[i].x_left = tmp;
			InitVector[i].x_right = tmp + li;
			InitVector[i].z_left = Func(tmp);
			InitVector[i].z_right = Func(tmp + li);
		}

		// считаем величины М каждого интервала
#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			InitVector[i].M = M(InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_right, InitVector[i].x_left);
		}

		M_big = InitVector[0].M;

		for (register int i = 1; i < procs; i++) // это надо б как-нить распараллелить (сомнительная выгода) // нахождение макса в массиве длины CPU_COUNT
		{
			if (M_big < InitVector[i].M) M_big = InitVector[i].M;
		}
		m_small = (M_big == 0) ? 1 : r * M_big;

#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			InitVector[i].R = R(m_small, InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_right, InitVector[i].x_left);
		}
		for (register int i = 0; i < procs; i++)
		{
			Intervals.push_back(InitVector[i]);
		}
		Intervals.sort();

		curr_left = Intervals.back().x_left;
		curr_right = Intervals.back().x_right;

		++k;
		// начало параллельного строгого алгоритма
		///////////////////////////////////////////
		//////////// Параллельный блок ////////////
		///////////////////////////////////////////
		while (k < _Steps && (abs(curr_left - curr_right) > eps))
		{
			Pointer = --Intervals.end();
			for (register int i = 0; i < procs * 2 - 2; i += 2) // забрали на работу CPU_COUNT участков с макс. хар-кой R в алгоритм из очереди
			{
				WorkVector[i] = *Pointer;
				Intervals.erase(Pointer--);
			}
			WorkVector[procs * 2 - 2] = *Pointer; Intervals.erase(Pointer);
			size_bank_intervals = Intervals.size();
			// ищем на каждом участке точку нового испытания x_new
#pragma omp parallel for shared(x_new) num_threads(procs)
			for (int i = 0; i < procs; i++)
			{
				x_new[i] = 0.5*(WorkVector[i * 2].x_right + WorkVector[i * 2].x_left)
					- (WorkVector[i * 2].z_right - WorkVector[i * 2].z_left) / (2 * m_small);
				WorkVector[i * 2 + 1].z_left = Func(x_new[i]);
			}

			for (register int i = 0; i < procs * 2; i += 2)
			{
				WorkVector[i + 1].x_right = WorkVector[i].x_right;
				WorkVector[i + 1].z_right = WorkVector[i].z_right;
				WorkVector[i + 1].x_left = WorkVector[i].x_right = x_new[i / 2];
				WorkVector[i].z_right = WorkVector[i + 1].z_left;
			}

			// всё скопировали, ничего не потеряли, пора считать характеристики
#pragma omp parallel for num_threads(procs)
			for (int i = 0; i < procs * 2; i++)
			{
				WorkVector[i].M = M(WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_right, WorkVector[i].x_left);
			}
			double M_max_array = WorkVector[0].M;
#pragma omp parallel sections num_threads(2)
			{

#pragma omp section
				{
					for (register int i = 1; i < procs * 2; i++) // это надо б как-нить распараллелить (сомнительная выгода) // нахождение макса в массиве длины CPU_COUNT * 2
					{
						if (M_max_array < WorkVector[i].M) M_max_array = WorkVector[i].M;
					}
				}
#pragma omp section
				{
					// теперь ищем M_max среди отрезоков, лежащих в приоритетной очереди
					if (size_bank_intervals != 0)
					{
						Pointer = Intervals.begin();
						double tmp_M; //= (*Pointer).M;
						double M_max_list = (*Pointer++).M;

						for (; Pointer != Intervals.end(); Pointer++)
						{
							tmp_M = (*Pointer).M;
							if (M_max_list < tmp_M) M_max_list = tmp_M;
						}

						// сравниваем M_max_array и M_max_list, выбираем большее
						M_big = (M_max_array < M_max_list) ? M_max_list : M_max_array;
					}
					else M_big = M_max_array;
				}
			}
			m_small = (M_big == 0) ? 1 : r * M_big;

#pragma omp parallel for num_threads(procs)
			for (int i = 0; i < procs * 2; i++)
			{
				WorkVector[i].R = R(m_small, WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_right, WorkVector[i].x_left);
			}
			// рассчитали R характеристики для рабочих отрезков
			// необходимо пересчитать R для отрезков в очереди
			if (size_bank_intervals != 0)
			{
				std::vector<block> R_vec(Intervals.begin(), Intervals.end());
#pragma omp parallel for num_threads(procs)
				for (register int i = 0; i < size_bank_intervals; i++)
				{
					R_vec[i].R = R(m_small, R_vec[i].z_right, R_vec[i].z_left, R_vec[i].x_right, R_vec[i].x_left);
				}
				register int j = 0;
				for (auto i = Intervals.begin(); i != Intervals.end(); i++)
				{
					i->R = R_vec[j++].R;
				}
				R_vec.clear();
			}
			for (register int i = 0; i < procs * 2; i++)
			{
				// возвращаем CPU_COUNT отрезков с подсчитанными характеристиками в приоритетную очередь для строгого алгоритма
				Intervals.push_back(WorkVector[i]);
			}
			Intervals.sort();
			++k;
			curr_left = Intervals.back().x_left;
			curr_right = Intervals.back().x_right;
		}
	}

	end = omp_get_wtime();
	time = end - start;

	omp_set_dynamic(1);

	solutionX = curr_right;
	solutionZ = Intervals.back().z_right;

	steps = k;
	// очистка
	InitVector.clear();
	WorkVector.clear();
	delete[] x_new;

	return solutionZ;
}


inline double parallel_global::R(const double &_m_small, const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev)
{
	return _m_small*(_x_curr - _x_prev) + pow(_z_curr - _z_prev, 2) / (_m_small*(_x_curr - _x_prev)) - 2 * (_z_curr + _z_prev);
}

inline double parallel_global::M(const double &_z_curr, const double &_z_prev, const double &_x_curr, const double &_x_prev)
{
	return abs((_z_curr - _z_prev) / (_x_curr - _x_prev));
}

inline double parallel_global::Func(const double &_x)
{
	//return sin(_x) + sin((10. * _x) / 3.); // (2.7 , 7.5 , r = 4.29 )
	//return (2 * pow((_x - 3), 2) + exp((pow(_x, 2) / 2))); // (-3, 3, r = 85)
	//return ((3 * _x - 1.4)*sin(18 * _x)); // //(0, 1.2 , r = 36)
	//return (sin(_x) + sin((10. * _x) / 3.) + log(_x) - 0.84*_x + 3); //(2.7 , 7.5 , r = 6)
	//return ((pow(_x, 2) - 5 * _x + 6) / (pow(_x, 2) + 1)); //(-5, 5 , r = 6.5)
	return 2 * sin(3 * _x) + 3 * cos(5 * _x); // (0,8 r = 2);
}