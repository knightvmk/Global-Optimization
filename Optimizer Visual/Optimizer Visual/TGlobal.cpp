#include "TGlobal.h"



//inline void Bubble(concurrency::array_view<block_1D> &arr, int _size) restrict(amp, cpu)
//{
//	block_1D tmp;
//
//	bool flag = true;
//	for (register int i = 0; i < _size; i++) 
//	{
//		flag = true;
//		for (register int j = 0; j < _size - i; j++) 
//		{
//			if (arr[j] > arr[j + 1]) 
//			{
//				flag = false;
//				tmp = arr[j]; 
//				arr[j] = arr[j + 1]; 
//				arr[j + 1] = tmp; 
//			}
//		}
//		if (flag) 	break;
//	}
//}
//
//inline void Bubble(concurrency::array_view<block_2D> &arr, int _size) restrict(amp, cpu)
//{
//	block_2D tmp;
//
//	bool flag = true;
//	for (register int i = 0; i < _size; i++)
//	{
//		flag = true;
//		for (register int j = 0; j < _size - i; j++)
//		{
//			if (arr[j] > arr[j + 1])
//			{
//				flag = false;
//				tmp = arr[j];
//				arr[j] = arr[j + 1];
//				arr[j + 1] = tmp;
//			}
//		}
//		if (flag) 	break;
//	}
//}
//
//inline void Quick(concurrency::array_view<block_1D> &s_arr, int first, int last) restrict(amp, cpu)
//{
//	int i = first, j = last;
//	block_1D x = s_arr[(first + last) / 2];
//
//	while (i <= j)
//	{
//		while (s_arr[i] < x) ++i;
//		while (s_arr[j] > x) --j;
//
//		if (i <= j)
//		{
//			if (s_arr[i] > s_arr[j])
//			{
//				block_1D tmp = s_arr[i];
//				s_arr[i] = s_arr[j];
//				s_arr[j] = tmp;
//			}
//			++i;
//			--j;
//		}
//	}
//
//	if (i < last)
//		Quick(s_arr, i, last);
//	if (first < j)
//		Quick(s_arr, first, j);
//}
//
//inline void Quick(concurrency::array_view<block_2D> &s_arr, int first, int last) restrict(amp, cpu)
//{
//	int i = first, j = last;
//	block_2D x = s_arr[(first + last) / 2];
//
//	while (i <= j)
//	{
//		while (s_arr[i] < x) ++i;
//		while (s_arr[j] > x) --j;
//
//		if (i <= j) 
//		{
//			if (s_arr[i] > s_arr[j])
//			{
//				block_2D tmp = s_arr[i];
//				s_arr[i] = s_arr[j];
//				s_arr[j] = tmp;
//			}
//			++i;
//			--j;
//		}
//	}
//
//	if (i < last)
//		Quick(s_arr, i, last);
//	if (first < j)
//		Quick(s_arr, first, j);
//}


bool operator<(const block_1D& i1, const block_1D& i2) 
{
	return (i1.R < i2.R) ? true : false;
}

bool operator>(const block_1D& i1, const block_1D& i2) 
{
	return (i1.R > i2.R) ? true : false;
}

bool operator<(const block_2D& i1, const block_2D& i2) 
{
	return (i1.R < i2.R) ? true : false;
}

bool operator>(const block_2D& i1, const block_2D& i2) 
{
	return (i1.R > i2.R) ? true : false;
}

inline double OneDimension::R(const double & _m_small, const double & _z_curr, const double & _z_prev, const double & _x_curr, const double & _x_prev)
{
	register double x_tmp = _x_curr - _x_prev;
	return _m_small*(x_tmp) + (pow(_z_curr - _z_prev, 2) / (_m_small*x_tmp)) - (2 * (_z_curr + _z_prev));
}

inline double OneDimension::M(const double & _z_curr, const double & _z_prev, const double & _x_curr, const double & _x_prev)
{
	return abs((_z_curr - _z_prev) / (_x_curr - _x_prev));
}

inline double OneDimension::Func(const double & _x, int &thread)
{
	//return sin(_x) + sin((10. * _x) / 3.); // (2.7 , 7.5 , r = 4.29 )
	//return (2 * pow((_x - 3), 2) + exp((pow(_x, 2) / 2))); // (-3, 3, r = 85)
	//return ((3 * _x - 1.4)*sin(18 * _x)); // //(0, 1.2 , r = 36)
	//return (sin(_x) + sin((10. * _x) / 3.) + log(_x) - 0.84*_x + 3); //(2.7 , 7.5 , r = 6)
	//return ((pow(_x, 2) - 5 * _x + 6) / (pow(_x, 2) + 1)); //(-5, 5 , r = 6.5)
	//return 2 * sin(3 * _x) + 3 * cos(5 * _x); // (0,8 r = 2);
	matexp[thread].Compile(math_exp_char, _x, 1);
	matexp[thread].Evaluate();
	return matexp[thread].GetResult();
}

inline double TwoDimension::R(const double & _m_small, const double & _z_curr, const double & _z_prev, const double & _x_curr, const double & _x_prev)
{
	register double x_tmp = _x_curr - _x_prev;
	return _m_small*(x_tmp)+(pow(_z_curr - _z_prev, 2) / (_m_small*x_tmp)) - (2 * (_z_curr + _z_prev));
}

inline double TwoDimension::M(const double & _z_curr, const double & _z_prev, const double & _x_curr, const double & _x_prev)
{
	return abs((_z_curr - _z_prev) / (_x_curr - _x_prev));
}

inline double TwoDimension::Func(const double *_y, int &thread)
{
	//return 2 * cos(_y[1]) + 3 * sin(_y[0]); //y[-5,0] x[-4,2] examin ~= (x_1.6 , y_-3.1)
	matexp[thread].Compile(math_exp_char, _y[0], _y[1]);
	matexp[thread].Evaluate();
	return matexp[thread].GetResult();
}

void TGlobal::InitCL()
{
#ifndef NOT_USE_CL
	if (can_use_cl)
		try // блок кода
	{
#ifdef USE_CL_200
		cl::string *str = new cl::string();
#elif defined(USE_CL_100)
		cl::STRING_CLASS *str = new cl::STRING_CLASS();
#endif
		// Получение доступных OpenCL платформ
		std::vector <cl::Platform > platforms;
		cl::Platform::get(&platforms);
		if (!platforms.empty())
		{
			for (auto iter = platforms.begin(); iter != platforms.end(); iter++)
			{
				iter->getInfo(CL_PLATFORM_NAME, str);
				avaible_platform_names.push_back(*str);
			}
#ifndef NOT_USE_CL
			avaible_platforms = platforms;
#endif

		}
		else return;
		//
		// Выбор устройства для запуска
		std::vector <cl::Device> device;
		for (int i = 0; i < platforms.size(); i++)
		{
			platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &device);
			for (auto iter = device.begin(); iter != device.end(); iter++)
			{
				iter->getInfo(CL_DEVICE_NAME, str);
				avaible_devices_names.push_back(*str);
			}
#ifndef NOT_USE_CL
			avaible_devices.push_back(device);
			device.clear();
#endif
		}
	}
	catch (cl::Error err) // блок ошибок
	{
		std::cerr << " ERROR : " << err.what() << "(" << err.err() << ")" << std::endl;
	}
	else return;
#endif
}

double Sequental_1D::Optimize(const double _left, const double _right, const int _N_max, const double _Eps)
{
	time = 0.L;

	double curr_left = _left;

	double curr_right = _right;

	// очистка списка
	Values_Z.clear();
	Arguments_X.clear();

	double start, end;

	int thread_num = omp_get_thread_num();

	/////
	double z_begin = Func(curr_left, thread_num); // левая граница
	double z_end = Func(curr_right, thread_num); // правая граница
									 //// прочие переменные
	double temp = 0.; // переменная для подсчёта 
	int k = 1; // количество шагов в Настоящем алгоритме
	int new_size; // тоже шаги
	double M_big, m_small = 1.; // Перменные для хранения итоговой константы "М" и константы "m"
	int R_max_index = 0; // Переменная индекса следующей точки испытания

	double new_point; // переменные для рабочего отрезка и точки нового испытания
					  ////
					  //// контейнеры для хранения значений и аргументов функции
	std::list<double> z; //значения Z
	std::list<double> x; //значения рабочих X
	std::list<double> M_vector; // вектор неких коэфф M для выбора m = max(M[i])
	std::list<double> R_vector; // вектор значений R для принятия решения точки следующего испытания

	block_1D val_z; // структура для хранения левой и правой границы отрезка с Z
	block_1D arg_x; // структура для хранения левой и правой границы отрезка с Х

				 //// блок итераторов для работы со списками
	std::list<double>::iterator place_M, place_R; // итераторы-указатели КУДА положить новые значения
	std::list<double>::iterator iter_x, iter_z, iter_R, iter_M;
	////
	M_big = M(z_end, z_begin, curr_right, curr_left);
	m_small = (M_big == 0) ? 1 : r * M_big;

	new_point = (curr_right + curr_left) / 2 - (z_end - z_begin) / (2 * m_small);

	x.push_back(curr_left);
	x.push_back(new_point);
	x.push_back(curr_right);

	iter_x = x.begin();
	for (int i = 0; i < 3; i++)
	{
		z.push_back(Func(*iter_x++, thread_num));
	}
	
	iter_z = z.begin();
	iter_x = x.begin();

	//// начало алгоритма // первый шаг, новая точка - средняя //

	//Начало времени отсчёта работы алгоритма
	start = omp_get_wtime();

	{
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			M_vector.push_back(M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // вставка новой М_big
		}

		double max = M_vector.front();

		if (max < M_vector.back()) max = M_vector.back();
		M_big = max; // поиск максимальной M_big

		m_small = (M_big == 0) ? 1 : r * M_big;

		R_max_index = 1;

		iter_z = z.begin();
		iter_x = x.begin();
		// считаем и вставляем новую характеристику для рабочего интервала и пересчитываем для получившегося "старого" интервала
		// старый интервал получился делением x(i-1)<x_new<x(i)
		// новый интервал (рабочий) x(i-1)<x_new
		// старый интервал x_new<x(i)
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;

			temp = R(m_small, val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left);
			R_vector.push_back(temp);
			if (i == 1) max = temp;

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

		curr_right = *iter_x;
		curr_left = *--iter_x;
	}

	//// начало работы АГП с испытанием в новой точке по правилу //
	while (k < _N_max && (abs(curr_left - curr_right) > _Eps))
	{
		new_size = 2 + k; // + 2 для быстрого алгоритма, потом будет корректировка

		iter_z = z.begin(); // найдём интервалы с максимальной характеристикой
		iter_x = x.begin();
		// ищем интервалы сдвигая итераторы
		for (register int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}
		// забираем значения для хранения и работы нового цикла
		val_z.z_right = *iter_z--;
		val_z.z_left = *iter_z++;

		//left = curr_left;
		//right = curr_right;
		// вычисление новой точки испытания x_new
		new_point = (curr_right + curr_left) / 2 - (val_z.z_right - val_z.z_left) / (2 * m_small);
		// сохраняем x_new и z_new
		z.insert(iter_z, Func(new_point, thread_num));
		x.insert(iter_x, new_point);

		iter_z = z.begin();
		iter_x = x.begin();
		place_M = M_vector.begin();
		place_R = R_vector.begin();

		// ищем интервал с лучшей характеристикой сдвигом итераторов
		for (int i = 0; i < R_max_index - 1; i++)
		{
			iter_z++;
			iter_x++;
			if (i == R_max_index - 1) break;
			place_M++;
			place_R++;
		}
		// старый интервал получился делением x(i-1)<x_new<x(i)
		// новый интервал (рабочий) x(i-1)<x_new
		// старый интервал x_new<x(i)
		for (int i = 0; i < 2; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			if (i == 0) M_vector.insert(place_M, M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // x(i - 1)<x_new
			else *place_M = M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left); // x_new<x(i)
		}

		iter_M = M_vector.begin();

		double max = *iter_M;
		temp = M_vector.back();

		for (int i = 0; i < new_size - 1; i++)
		{
			if (max < temp) max = temp;
			temp = *++iter_M;
		}
		M_big = max;

		m_small = (M_big == 0) ? 1 : r * M_big;

		// старый интервал получился делением [ x(i-1) < x_new < x(i) ]
		// новый интервал (рабочий) [ x(i-1) < x_new ]
		// старый интервал [ x_new < x(i) ]

		iter_z = z.begin(); iter_x = x.begin();

		R_max_index = 1;
		place_R = R_vector.begin();

		for (int i = 0; i < new_size; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			if (i == 0)	R_vector.insert(place_R, R(m_small, val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // x(i - 1)<x_new
			else *place_R++ = R(m_small, val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left); // x_new<x(i)
		}

		iter_R = R_vector.begin();
		max = *iter_R;
		for (register int i = 1; i < new_size; i++)
		{
			temp = *++iter_R;
			if (max < temp)
			{
				max = temp;
				R_max_index = i + 1; // нашли R(t) // R(i)
			}
		}

		iter_x = x.begin();
		iter_z = z.begin();

		for (register int i = 0; i < R_max_index; i++) // поиск точки глоб.минимума (интервала)
		{
			iter_z++;
			iter_x++;
		}

		// сохранение
		curr_right = *iter_x;
		curr_left = *--iter_x;
		++k;
	}

	end = omp_get_wtime();
	time = end - start;

	--new_size; // корректировка после работы алгоритма
	steps = new_size;

	// сохраняем для последовательной однопроцессорной версии

	solutionZ = *iter_z;
	solutionX = curr_right;
	Arguments_X = x;
	Values_Z = z;

	//return solutionX;
	return curr_right; // для параллельной версии
}

double Parallel_1D::Optimize(const double _left, const double _right, const int _N_max, const double _Eps)
{
	time = 0.L;

	double curr_left = _left;

	double curr_right = _right;

	// очистка списка
	//Values_Z.clear();
	//Arguments_X.clear();

	//double start, end;

	int thread_num = omp_get_thread_num();

	/////
	double z_begin = Func(curr_left, thread_num); // левая граница
	double z_end = Func(curr_right, thread_num); // правая граница
									 //// прочие переменные
	double temp = 0.; // переменная для подсчёта 
	int k = 1; // количество шагов в Настоящем алгоритме
	int new_size; // тоже шаги
	double M_big, m_small = 1.; // Перменные для хранения итоговой константы "М" и константы "m"
	int R_max_index = 0; // Переменная индекса следующей точки испытания

	double new_point; // переменные для рабочего отрезка и точки нового испытания
					  ////
					  //// контейнеры для хранения значений и аргументов функции
	std::list<double> z; //значения Z
	std::list<double> x; //значения рабочих X
	std::list<double> M_vector; // вектор неких коэфф M для выбора m = max(M[i])
	std::list<double> R_vector; // вектор значений R для принятия решения точки следующего испытания

	block_1D val_z; // структура для хранения левой и правой границы отрезка с Z
	block_1D arg_x; // структура для хранения левой и правой границы отрезка с Х


				 //// блок итераторов для работы со списками
	std::list<double>::iterator place_M, place_R; // итераторы-указатели КУДА положить новые значения
	std::list<double>::iterator iter_x, iter_z, iter_R, iter_M;
	////
	M_big = M(z_end, z_begin, curr_right, curr_left);
	m_small = (M_big == 0) ? 1 : r * M_big;

	new_point = (curr_right + curr_left) / 2 - (z_end - z_begin) / (2 * m_small);

	x.push_back(curr_left);
	x.push_back(new_point);
	x.push_back(curr_right);

	iter_x = x.begin();
	for (int i = 0; i < 3; i++)
	{
		z.push_back(Func(*iter_x++, thread_num));
	}

	iter_z = z.begin();
	iter_x = x.begin();

	//// начало алгоритма // первый шаг, новая точка - средняя //

	//Начало времени отсчёта работы алгоритма
	//start = omp_get_wtime();

	{
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			M_vector.push_back(M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // вставка новой М_big
		}

		double max = M_vector.front();

		if (max < M_vector.back()) max = M_vector.back();
		M_big = max; // поиск максимальной M_big

		m_small = (M_big == 0) ? 1 : r * M_big;

		R_max_index = 1;

		iter_z = z.begin();
		iter_x = x.begin();
		// считаем и вставляем новую характеристику для рабочего интервала и пересчитываем для получившегося "старого" интервала
		// старый интервал получился делением x(i-1)<x_new<x(i)
		// новый интервал (рабочий) x(i-1)<x_new
		// старый интервал x_new<x(i)
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;

			temp = R(m_small, val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left);
			R_vector.push_back(temp);
			if (i == 1) max = temp;

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

		curr_right = *iter_x;
		curr_left = *--iter_x;
	}

	//// начало работы АГП с испытанием в новой точке по правилу //
	while (k < _N_max && (abs(curr_left - curr_right) > _Eps))
	{
		new_size = 2 + k; // + 2 для быстрого алгоритма, потом будет корректировка

		iter_z = z.begin(); // найдём интервалы с максимальной характеристикой
		iter_x = x.begin();
		// ищем интервалы сдвигая итераторы
		for (register int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}
		// забираем значения для хранения и работы нового цикла
		val_z.z_right = *iter_z--;
		val_z.z_left = *iter_z++;

		//left = curr_left;
		//right = curr_right;
		// вычисление новой точки испытания x_new
		new_point = (curr_right + curr_left) / 2 - (val_z.z_right - val_z.z_left) / (2 * m_small);
		// сохраняем x_new и z_new
		z.insert(iter_z, Func(new_point, thread_num));
		x.insert(iter_x, new_point);

		iter_z = z.begin();
		iter_x = x.begin();
		place_M = M_vector.begin();
		place_R = R_vector.begin();

		// ищем интервал с лучшей характеристикой сдвигом итераторов
		for (int i = 0; i < R_max_index - 1; i++)
		{
			iter_z++;
			iter_x++;
			if (i == R_max_index - 1) break;
			place_M++;
			place_R++;
		}
		// старый интервал получился делением x(i-1)<x_new<x(i)
		// новый интервал (рабочий) x(i-1)<x_new
		// старый интервал x_new<x(i)
		for (int i = 0; i < 2; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			if (i == 0) M_vector.insert(place_M, M(val_z.z_right, val_z.z_left, arg_x.z_right, arg_x.x_left)); // x(i - 1)<x_new
			else *place_M = M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left); // x_new<x(i)
		}

		iter_M = M_vector.begin();

		double max = *iter_M;
		temp = M_vector.back();

		for (int i = 0; i < new_size - 1; i++)
		{
			if (max < temp) max = temp;
			temp = *++iter_M;
		}
		M_big = max;

		m_small = (M_big == 0) ? 1 : r * M_big;

		// старый интервал получился делением [ x(i-1) < x_new < x(i) ]
		// новый интервал (рабочий) [ x(i-1) < x_new ]
		// старый интервал [ x_new < x(i) ]

		iter_z = z.begin(); iter_x = x.begin();

		R_max_index = 1;
		place_R = R_vector.begin();

		for (int i = 0; i < new_size; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			if (i == 0)	R_vector.insert(place_R, R(m_small, val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // x(i - 1)<x_new
			else *place_R++ = R(m_small, val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left); // x_new<x(i)
		}

		iter_R = R_vector.begin();
		max = *iter_R;
		for (register int i = 1; i < new_size; i++)
		{
			temp = *++iter_R;
			if (max < temp)
			{
				max = temp;
				R_max_index = i + 1; // нашли R(t) // R(i)
			}
		}

		iter_x = x.begin();
		iter_z = z.begin();

		for (register int i = 0; i < R_max_index; i++) // поиск точки глоб.минимума (интервала)
		{
			iter_z++;
			iter_x++;
		}

		// сохранение
		curr_right = *iter_x;
		curr_left = *--iter_x;
		++k;
	}

	//end = omp_get_wtime();
	//time = end - start;

	--new_size; // корректировка после работы алгоритма
	steps = new_size;

#pragma omp critical
		{
			//short thread = omp_get_thread_num();
			parallel_x[thread_num] = x;
			parallel_z[thread_num] = z;
			parallel_steps += new_size;
		}

	//return solutionX;
	return curr_right; // для параллельной версии
}

double Parallel_1D::Optimize(const double _left, const double _right, const int _N_max, const double _Eps, const int _threads)
{
	time = 0.L;

	if (_threads < 0 || _threads > omp_get_max_threads())
	{
		throw("Incorrect number of CPU_THREADS. Must be > 0 and < MAX_THREADS");
		procs = omp_get_max_threads();
	}
	else procs = _threads;
	// выделяем память для хранения точек испытания работы алгоритма
	parallel_x.resize(procs);
	parallel_z.resize(procs);

	std::vector<block_1D> optimized_values(procs); // left = x, right = z;
	omp_set_dynamic(0);

	double start = omp_get_wtime();
	omp_set_num_threads(procs);

#pragma omp parallel for num_threads(procs)
	for (int i = 0; i < procs; i++)
	{
		optimized_values[i].x_right = Optimize(_left + (_right - _left)*i / procs, _left + (_right - _left)*(i + 1) / procs, _N_max / procs, _Eps);
		optimized_values[i].z_right = Func(optimized_values[i].x_right, i);
	}

	// найдём минимум
	double min_x = optimized_values[0].x_right;
	double min_z = optimized_values[0].z_right;

	double tmp;

#pragma omp parallel sections num_threads(2)
	{
#pragma omp section
		{
			for (short i = 1; i < procs; i++)
			{
				// находим минимум по x и по z
				tmp = optimized_values[i].z_right;
				if (min_z > tmp)
				{
					min_x = optimized_values[i].x_right;
					min_z = optimized_values[i].z_right;
				}
				// сохраняем все точки испытаний всех отрезков в один результирующий вектор
			}
		}
#pragma omp section
		{
			// очищаем от мусора предыдущего эксперимента
			if (!Arguments_X.empty())
			{
				Arguments_X.clear();
				Values_Z.clear();
			}
		}
	}


	// сохраняем в чистый список результатов все списки результатов работы потоков
	for (short i = 0; i < procs; i++)
	{
		parallel_x[i].sort();
		Arguments_X.merge(parallel_x[i]);
		Values_Z.splice(Values_Z.end(), parallel_z[i]);
	}

	omp_set_dynamic(1);

	double end = omp_get_wtime();
	time = end - start;

	steps = parallel_steps;

	solutionX = min_x;
	solutionZ = min_z;
	return min_x;
}

double Parallel_1D::Optimize(const double _Epsilon, const int _Steps, const int _thread_count)
{
	time = 0.L;
	if (_thread_count < 0 || _thread_count > omp_get_max_threads())
	{
		throw("Incorrect number of CPU_THREADS. Must be > 0 and < MAX_THREADS");
		procs = omp_get_max_threads();
	}
	else procs = _thread_count;
	//int procs = omp_get_max_threads();
	int k = 0; // количество шагов
	double eps = _Epsilon;
	double curr_left, curr_right; // текущее значение слева и справа для оценки разности
	double start, end; // время работы
	double m_small = 1., M_big;
	double *x_new = new double[procs];

	curr_left = left; // инициализация границ
	curr_right = right; // инициализация границ

	std::vector<block_1D> InitVector(procs);
	std::vector<block_1D> WorkVector(procs * 2);
	std::priority_queue<block_1D> Queue;
	std::list<block_1D> Intervals;

	int size_bank_intervals; // размер list<block_1D> Intervals

	std::list<block_1D>::iterator Pointer;


	// начало первичной инициализации параллельного алгоритма
	// деление отрезка на звенья по количеству ЦП
	// вычисление R характеристик для звеньев - интервалов
	double len = right - left; // расстояние поиска
	double li = len / procs; // длина i-го участка
	register double tmp = 0.;

	omp_set_dynamic(0);

	start = omp_get_wtime();

	omp_set_num_threads(procs);
	{
#pragma omp parallel num_threads(procs)
		{
#pragma omp for firstprivate(tmp)
			for (int i = 0; i < procs; i++)
			{
				tmp = (right - len) + (i * li); // _b-(b-a)_ + _расстояние_участка_ * _номер_участка_;
				InitVector[i].x_left = tmp;
				InitVector[i].x_right = tmp + li;
				InitVector[i].z_left = Func(tmp , i);
				InitVector[i].z_right = Func(tmp + li , i);
			}

			// считаем величины М каждого интервала
#pragma omp for 
			for (int i = 0; i < procs; i++)
			{
				InitVector[i].M = M(InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_right, InitVector[i].x_left);
			}
#pragma omp single
			{
				M_big = InitVector[0].M;

				for (register int i = 1; i < procs; i++) // это надо б как-нить распараллелить (сомнительная выгода) // нахождение макса в массиве длины CPU_COUNT
				{
					M_big = (M_big < InitVector[i].M) ? InitVector[i].M : M_big;
				}
				m_small = (M_big == 0) ? 1 : r * M_big;
			}
#pragma omp for 
			for (int i = 0; i < procs; i++)
			{
				InitVector[i].R = R(m_small, InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_right, InitVector[i].x_left);
			}
#pragma omp single
			for (register int i = 0; i < procs; i++)
			{
				Intervals.push_back(InitVector[i]);
			}
		}
		Intervals.sort();

		curr_left = Intervals.back().x_left;
		curr_right = Intervals.back().x_right;

		++k;
		// начало параллельного строгого алгоритма
		///////////////////////////////////////////
		//////////// Параллельный блок ////////////
		///////////////////////////////////////////
		size_bank_intervals = procs;
		while (k < _Steps/procs && (abs(curr_left - curr_right) > eps))
		{
			//Pointer = --Intervals.end();
			size_bank_intervals -= procs;
			for (register int i = 0; i < procs * 2 - 2; i += 2) // забрали на работу CPU_COUNT участков с макс. хар-кой R в алгоритм из очереди
			{
				WorkVector[i] = Intervals.back();
				Intervals.pop_back();
				//Intervals.erase(Pointer--);
			}
			WorkVector[procs * 2 - 2] = Intervals.back(); Intervals.pop_back();
			//size_bank_intervals = Intervals.size();
			
			// ищем на каждом участке точку нового испытания x_new
#pragma omp parallel num_threads(procs)
			{
#pragma omp for 
				for (int i = 0; i < procs; i++)
				{
					x_new[i] = 0.5*(WorkVector[i * 2].x_right + WorkVector[i * 2].x_left)
						- (WorkVector[i * 2].z_right - WorkVector[i * 2].z_left) / (2 * m_small);
					WorkVector[i * 2 + 1].z_left = Func(x_new[i] , i);
				}
#pragma omp single
				for (register int i = 0; i < procs * 2; i += 2)
				{
					WorkVector[i + 1].x_right = WorkVector[i].x_right;
					WorkVector[i + 1].z_right = WorkVector[i].z_right;
					WorkVector[i + 1].x_left = WorkVector[i].x_right = x_new[i / 2];
					WorkVector[i].z_right = WorkVector[i + 1].z_left;
				}

				// всё скопировали, ничего не потеряли, пора считать характеристики
#pragma omp for 
				for (int i = 0; i < procs * 2; i++)
				{
					WorkVector[i].M = M(WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_right, WorkVector[i].x_left);
				}
			}
			register double M_max_array = WorkVector[0].M;
#pragma omp parallel sections num_threads(2)
			{
#pragma omp section
				{
					for (register int i = 1; i < procs * 2; i++) // это надо б как-нить распараллелить (сомнительная выгода) // нахождение макса в массиве длины CPU_COUNT * 2
					{
						M_max_array = (M_max_array < WorkVector[i].M) ? WorkVector[i].M : M_max_array;
					}
				}
#pragma omp section
				{
					// теперь ищем M_max среди отрезоков, лежащих в приоритетной очереди
					if (size_bank_intervals != 0)
					{
						double tmp_M = -100; 
						for (auto Pointer = Intervals.begin(); Pointer != Intervals.end(); Pointer++)
						{
							tmp_M = (tmp_M < Pointer->M) ? Pointer->M : tmp_M;
						}

						// сравниваем M_max_array и M_max_list, выбираем большее
						M_big = (M_max_array < tmp_M) ? tmp_M : M_max_array;
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
				std::vector<block_1D> R_vec(Intervals.begin(), Intervals.end());
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
			size_bank_intervals += 2 * procs;
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

	Values_Z.clear();
	Arguments_X.clear();

	// сохранение результатов испытаний
	//auto iter_vec = InitVector.begin();
	//Values_Z.push_back(iter_vec->z_left);
	//Arguments_X.push_back(iter_vec->x_left);

	//for (; iter_vec != InitVector.end(); iter_vec++)
	//{
	//	Values_Z.push_back(iter_vec->z_right);
	//	Arguments_X.push_back(iter_vec->x_right);
	//}

	auto iter = Intervals.begin();
	//Values_Z.push_back(iter->z_left);
	//Arguments_X.push_back(iter->x_left); iter++;

	for (; iter != Intervals.end(); iter++)
	{
		Values_Z.push_back(iter->z_left);
		Values_Z.push_back(iter->z_right);
		Arguments_X.push_back(iter->x_left);
		Arguments_X.push_back(iter->x_right);
	}
	Values_Z.unique(); Arguments_X.unique();

	steps = k*procs;
	// очистка
	InitVector.clear();
	WorkVector.clear();
	Intervals.clear();
	delete[] x_new;
	//InitVector.~vector();
	//WorkVector.~vector();
	//Intervals.~list();

	return solutionZ;
}

double Sequental_2D::Optimize(const double *_Left, const double *_Right, const int _N_max, const double _Eps)
{
	time = 0.L;
	procs = 1;

	double curr_left_y = Left[0], curr_left_x = Left[1];
	double curr_right_y = Right[0], curr_right_x = Right[1];

	// новая метрика |X(i) - X(i-1)|^1/N , N = 2

	// очистка списка
	Values_Z.clear();
	Arguments_X.clear();
	Arguments_Y.clear();

	double start, end; // time
	double curr_left, curr_right; // x[0;1]

	double *A = new double[2]; // ограничения слева
	double *B = new double[2]; // ограничения справа

	A[0] = Left[0]; // x
	A[1] = Left[1]; // y
	B[0] = Right[0]; // x
	B[1] = Right[1]; // y


	int k = 0; // количество шагов
	double eps = _Eps;
	double m_small = 1., M_big;
	double *x_img = new double;
	double *y_arr = new double[2];
	block_2D WorkVector[2];
	std::list<block_2D> Intervals;

	std::list<block_2D>::iterator Pointer;

	TEvolvent evolve(2, 10);
	evolve.SetBounds(A, B);

	double tmp[2];

	int thread_num = omp_get_thread_num();

	start = omp_get_wtime();
	{

		tmp[0] = curr_left_x;
		tmp[1] = curr_left_y;

		WorkVector[0].x_left = curr_left_x;
		WorkVector[0].x_right = curr_right_x;

		WorkVector[0].x_img_left = 0;
		WorkVector[0].x_img_right = 1;

		WorkVector[0].y_left = curr_left_y;
		WorkVector[0].y_right = curr_right_y;

		evolve.GetImage(0, tmp);

		WorkVector[0].z_left = Func(tmp, thread_num);

		evolve.GetImage(1, tmp);

		WorkVector[0].z_right = Func(tmp, thread_num);


		// считаем величины М каждого интервала

		WorkVector[0].M = M(WorkVector[0].z_right, WorkVector[0].z_left, WorkVector[0].x_img_right, WorkVector[0].x_img_left);

		M_big = WorkVector[0].M;

		m_small = (M_big == 0) ? 1 : r * M_big;

		WorkVector[0].R = R(m_small, WorkVector[0].z_right, WorkVector[0].z_left, WorkVector[0].x_img_right, WorkVector[0].x_img_left);
		Intervals.push_back(WorkVector[0]);

		curr_left = Intervals.back().x_img_left;
		curr_right = Intervals.back().x_img_right;

		++k;

		while (k < _N_max && (abs(curr_left - curr_right) > eps))
		{
			WorkVector[0] = Intervals.back(); Intervals.pop_back();
			// ищем на каждом участке точку нового испытания x_new

			*x_img = 0.5*(WorkVector[0].x_img_right + WorkVector[0].x_img_left)
				- (WorkVector[0].z_right - WorkVector[0].z_left) / (2 * m_small);
			evolve.GetImage(*x_img, y_arr);
			WorkVector[1].x_left = y_arr[0];
			WorkVector[1].y_left = y_arr[1];
			WorkVector[1].z_left = Func(y_arr, thread_num);

			WorkVector[1].x_right = WorkVector[0].x_right;
			WorkVector[1].y_right = WorkVector[0].y_right;
			WorkVector[1].x_img_right = WorkVector[0].x_img_right;
			WorkVector[1].z_right = WorkVector[0].z_right;
			// перенос совершён
			WorkVector[1].x_img_left = WorkVector[0].x_img_right = *x_img;

			WorkVector[0].x_right = WorkVector[1].x_left;
			WorkVector[0].y_right = WorkVector[1].y_left;
			WorkVector[0].z_right = WorkVector[1].z_left;


			// всё скопировали, ничего не потеряли, пора считать характеристики

			for (int i = 0; i < 2; i++)
			{
				WorkVector[i].M = M(WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_img_right, WorkVector[i].x_img_left);
			}
			M_big = WorkVector[0].M;

			for (register int i = 1; i < 2; i++) // это надо б как-нить распараллелить (сомнительная выгода) // нахождение макса в массиве длины CPU_COUNT * 2
			{
				M_big = (M_big < WorkVector[i].M) ? WorkVector[i].M : M_big;
			}

			// теперь ищем M_max среди отрезков, лежащих в приоритетной очереди
			Pointer = Intervals.begin();
			for (; Pointer != Intervals.end(); Pointer++)
			{
				M_big = (M_big < Pointer->M) ? Pointer->M : M_big;
			}

			m_small = (M_big == 0) ? 1 : r * M_big;

			for (register int i = 0; i < 2; i++)
			{
				Intervals.push_back(WorkVector[i]);
			}
			Pointer = Intervals.begin();
			for (; Pointer != Intervals.end(); Pointer++)
			{
				Pointer->R = R(m_small, Pointer->z_right, Pointer->z_left, Pointer->x_img_right, Pointer->x_img_left);
			}
			Intervals.sort();
			++k;
			curr_left = Intervals.back().x_img_left;
			curr_right = Intervals.back().x_img_right;
		}
	}

	end = omp_get_wtime();
	time = end - start;

	solutionX = Intervals.back().x_right;
	solutionY = Intervals.back().y_right;
	solutionZ = Intervals.back().z_right;

	Values_Z.clear();
	Arguments_X.clear();
	Arguments_Y.clear();

	// сохранение результатов испытаний

	auto iter = Intervals.begin();

	for (; iter != Intervals.end(); iter++)
	{
		Values_Z.push_back(iter->z_left);
		Values_Z.push_back(iter->z_right);
		Arguments_X.push_back(iter->x_left);
		Arguments_X.push_back(iter->x_right);
		Arguments_Y.push_back(iter->y_left);
		Arguments_Y.push_back(iter->y_right);
	}
	Values_Z.unique(); Arguments_X.unique(); Arguments_Y.unique();

	steps = k;
	// очистка


	Intervals.clear();
	delete x_img;
	delete[]A;
	delete[]B;
	delete[]y_arr;

	return solutionZ;
}

double Parallel_2D::Optimize(const double _Epsilon, const int _Steps, const int _thread_count)
{
	time = 0.L;

	if (_thread_count < 0 || _thread_count > omp_get_max_threads())
	{
		throw("Incorrect number of CPU_THREADS. Must be > 0 and < MAX_THREADS");
		procs = omp_get_max_threads();
	}
	else procs = _thread_count;
	//int procs = omp_get_max_threads();

	double curr_left_y = Left[0], curr_left_x = Left[1];
	double curr_right_y = Right[0], curr_right_x = Right[1];

	// новая метрика |X(i) - X(i-1)|^1/N , N = 2

	// очистка списка
	Values_Z.clear();
	Arguments_X.clear();
	Arguments_Y.clear();

	double start, end; // time
	double curr_left, curr_right; // x[0;1]

	double *A = new double[2]; // ограничения слева
	double *B = new double[2]; // ограничения справа

	A[0] = Left[0]; // x
	A[1] = Left[1]; // y
	B[0] = Right[0]; // x
	B[1] = Right[1]; // y


	int k = 0; // количество шагов
	double eps = _Epsilon;
	double m_small = 1., M_big;
	double *x_img = new double[procs];
	double *y_arr = new double[2];

	std::vector<block_2D> InitVector(procs);
	std::vector<block_2D> WorkVector(procs * 2);
	std::list<block_2D> Intervals;

	int size_bank_intervals; // размер list<block_2D> Intervals

	std::list<block_2D>::iterator Pointer;

	double right_x, left_x, right_x_img, left_x_img, right_y, left_y;

	right_x = B[0];
	right_y = B[1];
	left_x = A[0];
	left_y = A[1];

	omp_set_dynamic(0);

	std::vector<TEvolvent*> evolve(procs);
//#pragma omp parallel for num_threads(procs)
	for (int i = 0; i < procs; i++)
	{
		evolve[i] = new TEvolvent(2, 20);
		evolve[i]->SetBounds(A, B);
	}
	left_x_img = x_img[0] = 0;
	right_x_img = x_img[procs - 1] = 1;

	// начало первичной инициализации параллельного алгоритма
	// деление отрезка на звенья по количеству ЦП
	// вычисление R характеристик для звеньев - интервалов
	double len_x = right_x - left_x; // расстояние поиска x
	double li_x = len_x / procs; // длина i-го участка x

	double len_x_img = right_x_img - left_x_img; // расстояние поиска x
	double li_x_img = len_x_img / procs; // длина i-го участка x

	double len_y = right_y - left_y; // расстояние поиска y
	double li_y = len_y / procs; // длина i-го участка y
	double tmp[2];
	double tmp_img = 0;

	start = omp_get_wtime();

	omp_set_num_threads(procs);
	{
#pragma omp parallel for firstprivate(tmp, tmp_img) num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			tmp[0] = (right_x - len_x) + (i * li_x); // _b-(b-a)_ + _расстояние_участка_ * _номер_участка_;
			tmp[1] = (right_y - len_y) + (i * li_y); // _b-(b-a)_ + _расстояние_участка_ * _номер_участка_;
			tmp_img = (right_x_img - len_x_img) + (i * li_x_img); // _b-(b-a)_ + _расстояние_участка_ * _номер_участка_;

			InitVector[i].x_left = tmp[0];
			InitVector[i].x_right = tmp[0] + li_x;

			InitVector[i].x_img_left = tmp_img;
			InitVector[i].x_img_right = tmp_img + li_x_img;

			InitVector[i].y_left = tmp[1];
			InitVector[i].y_right = tmp[1] + li_y;

			evolve[i]->GetImage(tmp_img, tmp);

			InitVector[i].z_left = Func(tmp, i);

			evolve[i]->GetImage(tmp_img + li_x_img, tmp);

			InitVector[i].z_right = Func(tmp, i);
		}
		double tmp1 = procs;
		// считаем величины М каждого интервала
#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			InitVector[i].M = M(InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_img_right, InitVector[i].x_img_left);
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
			InitVector[i].R = R(m_small, InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_img_right, InitVector[i].x_img_left);
		}
		for (register int i = 0; i < procs; i++)
		{
			Intervals.push_back(InitVector[i]);
		}
		Intervals.sort();

		curr_left = Intervals.back().x_img_left;
		curr_right = Intervals.back().x_img_right;

		++k;
		// начало параллельного строгого алгоритма
		///////////////////////////////////////////
		//////////// Параллельный блок ////////////
		///////////////////////////////////////////
		size_bank_intervals = procs;
		while (k < _Steps/procs && (abs(curr_left - curr_right) > eps))
		{
			size_bank_intervals -= procs;
			for (register int i = 0; i < procs * 2 - 2; i += 2) // забрали на работу CPU_COUNT участков с макс. хар-кой R в алгоритм из очереди
			{
				WorkVector[i] = Intervals.back();
				Intervals.pop_back();
			}
			WorkVector[procs * 2 - 2] = Intervals.back(); Intervals.pop_back();

			// ищем на каждом участке точку нового испытания x_new
#pragma omp parallel for shared(x_img)  num_threads(procs)
			for (int i = 0; i < procs; i++)
			{
				register double y_arr[2];
				x_img[i] = 0.5*(WorkVector[i * 2].x_img_right + WorkVector[i * 2].x_img_left)
					- (WorkVector[i * 2].z_right - WorkVector[i * 2].z_left) / (2 * m_small);
				evolve[i]->GetImage(x_img[i], y_arr);
				WorkVector[i * 2 + 1].x_left = y_arr[0];
				WorkVector[i * 2 + 1].y_left = y_arr[1];
				WorkVector[i * 2 + 1].z_left = Func(y_arr, i);
			}

			for (register int i = 0; i < procs * 2; i += 2)
			{
				WorkVector[i + 1].x_right = WorkVector[i].x_right;
				WorkVector[i + 1].y_right = WorkVector[i].y_right;
				WorkVector[i + 1].x_img_right = WorkVector[i].x_img_right;
				WorkVector[i + 1].z_right = WorkVector[i].z_right;
				// перенос совершён
				WorkVector[i + 1].x_img_left = WorkVector[i].x_img_right = x_img[i / 2];

				WorkVector[i].x_right = WorkVector[i + 1].x_left;
				WorkVector[i].y_right = WorkVector[i + 1].y_left;
				WorkVector[i].z_right = WorkVector[i + 1].z_left;
			}

			// всё скопировали, ничего не потеряли, пора считать характеристики
#pragma omp parallel for num_threads(procs)
			for (int i = 0; i < procs * 2; i++)
			{
				WorkVector[i].M = M(WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_img_right, WorkVector[i].x_img_left);
			}
			register double M_max_array = WorkVector[0].M;
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
						register double tmp_M; //= (*Pointer).M;
						register double M_max_list = (*Pointer++).M;

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
				WorkVector[i].R = R(m_small, WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_img_right, WorkVector[i].x_img_left);
			}
			// рассчитали R характеристики для рабочих отрезков
			// необходимо пересчитать R для отрезков в очереди
			if (size_bank_intervals != 0)
			{
				std::vector<block_2D> R_vec(Intervals.begin(), Intervals.end());
#pragma omp parallel for num_threads(procs)
				for (register int i = 0; i < size_bank_intervals; i++)
				{
					R_vec[i].R = R(m_small, R_vec[i].z_right, R_vec[i].z_left, R_vec[i].x_img_right, R_vec[i].x_img_left);
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
			curr_left = Intervals.back().x_img_left;
			curr_right = Intervals.back().x_img_right;
		}
	}

	end = omp_get_wtime();
	time = end - start;

	omp_set_dynamic(1);

	solutionX = Intervals.back().x_right;
	solutionY = Intervals.back().y_right;
	solutionZ = Intervals.back().z_right;

	Values_Z.clear();
	Arguments_X.clear();
	Arguments_Y.clear();

	// сохранение результатов испытаний

	auto iter = Intervals.begin();

	for (; iter != Intervals.end(); iter++)
	{
		Values_Z.push_back(iter->z_left);
		Values_Z.push_back(iter->z_right);
		Arguments_X.push_back(iter->x_left);
		Arguments_X.push_back(iter->x_right);
		Arguments_Y.push_back(iter->y_left);
		Arguments_Y.push_back(iter->y_right);
	}
	Values_Z.unique(); Arguments_X.unique(); Arguments_Y.unique();

	steps = k*procs;
	// очистка

	InitVector.clear();
	WorkVector.clear();
	Intervals.clear();
	evolve.clear();
	delete[]x_img;
	delete[]A;
	delete[]B;
	delete[]y_arr;
	//InitVector.~vector();
	//WorkVector.~vector();
	//Intervals.~list();

	return solutionZ;
}


//////////////////////////// OLD //////////////////////////////////////
