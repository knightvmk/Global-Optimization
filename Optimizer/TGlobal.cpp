#include "TGlobal.h"

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

inline double OneDimension::Func(const double & _x)
{
	//return sin(_x) + sin((10. * _x) / 3.); // (2.7 , 7.5 , r = 4.29 )
	//return (2 * pow((_x - 3), 2) + exp((pow(_x, 2) / 2))); // (-3, 3, r = 85)
	//return ((3 * _x - 1.4)*sin(18 * _x)); // //(0, 1.2 , r = 36)
	//return (sin(_x) + sin((10. * _x) / 3.) + log(_x) - 0.84*_x + 3); //(2.7 , 7.5 , r = 6)
	//return ((pow(_x, 2) - 5 * _x + 6) / (pow(_x, 2) + 1)); //(-5, 5 , r = 6.5)
	return 2 * sin(3 * _x) + 3 * cos(5 * _x); // (0,8 r = 2);
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

inline double TwoDimension::Func(const double *_y)
{
	// _y[1] = y ; _y[0] = x;
	return 2 * cos(_y[1]) + 3 * sin(_y[0]); //y[-5,0] x[-4,2] examin ~= (x_1.6 , y_-3.1)
}

void TGlobal::InitCL()
{
	if (can_use_cl)
		try // ���� ����
	{
#ifdef USE_CL_200
		cl::string *str = new cl::string();
#elif defined(USE_CL_100)
		cl::STRING_CLASS *str = new cl::STRING_CLASS();
#endif
		// ��������� ��������� OpenCL ��������
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
		// ����� ���������� ��� �������
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
	catch (cl::Error err) // ���� ������
	{
		std::cerr << " ERROR : " << err.what() << "(" << err.err() << ")" << std::endl;
	}
	else return;
}

double Sequental_1D::Optimize(const double _left, const double _right, const int _N_max, const double _Eps)
{
	time = 0.L;

	double curr_left = _left;

	double curr_right = _right;

	// ������� ������
	Values_Z.clear();
	Arguments_X.clear();

	double start, end;

	/////
	double z_begin = Func(curr_left); // ����� �������
	double z_end = Func(curr_right); // ������ �������
									 //// ������ ����������
	double temp = 0.; // ���������� ��� �������� 
	int k = 1; // ���������� ����� � ��������� ���������
	int new_size; // ���� ����
	double M_big, m_small = 1.; // ��������� ��� �������� �������� ��������� "�" � ��������� "m"
	int R_max_index = 0; // ���������� ������� ��������� ����� ���������

	double new_point; // ���������� ��� �������� ������� � ����� ������ ���������
					  ////
					  //// ���������� ��� �������� �������� � ���������� �������
	std::list<double> z; //�������� Z
	std::list<double> x; //�������� ������� X
	std::list<double> M_vector; // ������ ����� ����� M ��� ������ m = max(M[i])
	std::list<double> R_vector; // ������ �������� R ��� �������� ������� ����� ���������� ���������

	block_1D val_z; // ��������� ��� �������� ����� � ������ ������� ������� � Z
	block_1D arg_x; // ��������� ��� �������� ����� � ������ ������� ������� � �


				 //// ���� ���������� ��� ������ �� ��������
	std::list<double>::iterator place_M, place_R; // ���������-��������� ���� �������� ����� ��������
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
		z.push_back(Func(*iter_x++));
	}

	iter_z = z.begin();
	iter_x = x.begin();

	//// ������ ��������� // ������ ���, ����� ����� - ������� //

	//������ ������� ������� ������ ���������
	start = omp_get_wtime();

	{
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			M_vector.push_back(M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // ������� ����� �_big
		}

		double max = M_vector.front();

		if (max < M_vector.back()) max = M_vector.back();
		M_big = max; // ����� ������������ M_big

		m_small = (M_big == 0) ? 1 : r * M_big;

		R_max_index = 1;

		iter_z = z.begin();
		iter_x = x.begin();
		// ������� � ��������� ����� �������������� ��� �������� ��������� � ������������� ��� ������������� "�������" ���������
		// ������ �������� ��������� �������� x(i-1)<x_new<x(i)
		// ����� �������� (�������) x(i-1)<x_new
		// ������ �������� x_new<x(i)
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

	//// ������ ������ ��� � ���������� � ����� ����� �� ������� //
	while (k < _N_max && (abs(curr_left - curr_right) > _Eps))
	{
		new_size = 2 + k; // + 2 ��� �������� ���������, ����� ����� �������������

		iter_z = z.begin(); // ����� ��������� � ������������ ���������������
		iter_x = x.begin();
		// ���� ��������� ������� ���������
		for (register int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}
		// �������� �������� ��� �������� � ������ ������ �����
		val_z.z_right = *iter_z--;
		val_z.z_left = *iter_z++;

		//left = curr_left;
		//right = curr_right;
		// ���������� ����� ����� ��������� x_new
		new_point = (curr_right + curr_left) / 2 - (val_z.z_right - val_z.z_left) / (2 * m_small);
		// ��������� x_new � z_new
		z.insert(iter_z, Func(new_point));
		x.insert(iter_x, new_point);

		iter_z = z.begin();
		iter_x = x.begin();
		place_M = M_vector.begin();
		place_R = R_vector.begin();

		// ���� �������� � ������ ��������������� ������� ����������
		for (int i = 0; i < R_max_index - 1; i++)
		{
			iter_z++;
			iter_x++;
			if (i == R_max_index - 1) break;
			place_M++;
			place_R++;
		}
		// ������ �������� ��������� �������� x(i-1)<x_new<x(i)
		// ����� �������� (�������) x(i-1)<x_new
		// ������ �������� x_new<x(i)
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

		// ������ �������� ��������� �������� [ x(i-1) < x_new < x(i) ]
		// ����� �������� (�������) [ x(i-1) < x_new ]
		// ������ �������� [ x_new < x(i) ]

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
				R_max_index = i + 1; // ����� R(t) // R(i)
			}
		}

		iter_x = x.begin();
		iter_z = z.begin();

		for (register int i = 0; i < R_max_index; i++) // ����� ����� ����.�������� (���������)
		{
			iter_z++;
			iter_x++;
		}

		// ����������
		curr_right = *iter_x;
		curr_left = *--iter_x;
		++k;
	}

	end = omp_get_wtime();
	time = end - start;

	--new_size; // ������������� ����� ������ ���������
	steps = new_size;

	// ��������� ��� ���������������� ���������������� ������

	solutionZ = *iter_z;
	solutionX = curr_right;
	Arguments_X = x;
	Values_Z = z;

	//return solutionX;
	return curr_right; // ��� ������������ ������
}

double Parallel_1D::Optimize(const double _left, const double _right, const int _N_max, const double _Eps)
{
	time = 0.L;

	double curr_left = _left;

	double curr_right = _right;

	// ������� ������
	//Values_Z.clear();
	//Arguments_X.clear();

	//double start, end;

	/////
	double z_begin = Func(curr_left); // ����� �������
	double z_end = Func(curr_right); // ������ �������
									 //// ������ ����������
	double temp = 0.; // ���������� ��� �������� 
	int k = 1; // ���������� ����� � ��������� ���������
	int new_size; // ���� ����
	double M_big, m_small = 1.; // ��������� ��� �������� �������� ��������� "�" � ��������� "m"
	int R_max_index = 0; // ���������� ������� ��������� ����� ���������

	double new_point; // ���������� ��� �������� ������� � ����� ������ ���������
					  ////
					  //// ���������� ��� �������� �������� � ���������� �������
	std::list<double> z; //�������� Z
	std::list<double> x; //�������� ������� X
	std::list<double> M_vector; // ������ ����� ����� M ��� ������ m = max(M[i])
	std::list<double> R_vector; // ������ �������� R ��� �������� ������� ����� ���������� ���������

	block_1D val_z; // ��������� ��� �������� ����� � ������ ������� ������� � Z
	block_1D arg_x; // ��������� ��� �������� ����� � ������ ������� ������� � �


				 //// ���� ���������� ��� ������ �� ��������
	std::list<double>::iterator place_M, place_R; // ���������-��������� ���� �������� ����� ��������
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
		z.push_back(Func(*iter_x++));
	}

	iter_z = z.begin();
	iter_x = x.begin();

	//// ������ ��������� // ������ ���, ����� ����� - ������� //

	//������ ������� ������� ������ ���������
	//start = omp_get_wtime();

	{
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x++;
			arg_x.x_right = *iter_x;
			M_vector.push_back(M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // ������� ����� �_big
		}

		double max = M_vector.front();

		if (max < M_vector.back()) max = M_vector.back();
		M_big = max; // ����� ������������ M_big

		m_small = (M_big == 0) ? 1 : r * M_big;

		R_max_index = 1;

		iter_z = z.begin();
		iter_x = x.begin();
		// ������� � ��������� ����� �������������� ��� �������� ��������� � ������������� ��� ������������� "�������" ���������
		// ������ �������� ��������� �������� x(i-1)<x_new<x(i)
		// ����� �������� (�������) x(i-1)<x_new
		// ������ �������� x_new<x(i)
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

	//// ������ ������ ��� � ���������� � ����� ����� �� ������� //
	while (k < _N_max && (abs(curr_left - curr_right) > _Eps))
	{
		new_size = 2 + k; // + 2 ��� �������� ���������, ����� ����� �������������

		iter_z = z.begin(); // ����� ��������� � ������������ ���������������
		iter_x = x.begin();
		// ���� ��������� ������� ���������
		for (register int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}
		// �������� �������� ��� �������� � ������ ������ �����
		val_z.z_right = *iter_z--;
		val_z.z_left = *iter_z++;

		//left = curr_left;
		//right = curr_right;
		// ���������� ����� ����� ��������� x_new
		new_point = (curr_right + curr_left) / 2 - (val_z.z_right - val_z.z_left) / (2 * m_small);
		// ��������� x_new � z_new
		z.insert(iter_z, Func(new_point));
		x.insert(iter_x, new_point);

		iter_z = z.begin();
		iter_x = x.begin();
		place_M = M_vector.begin();
		place_R = R_vector.begin();

		// ���� �������� � ������ ��������������� ������� ����������
		for (int i = 0; i < R_max_index - 1; i++)
		{
			iter_z++;
			iter_x++;
			if (i == R_max_index - 1) break;
			place_M++;
			place_R++;
		}
		// ������ �������� ��������� �������� x(i-1)<x_new<x(i)
		// ����� �������� (�������) x(i-1)<x_new
		// ������ �������� x_new<x(i)
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

		// ������ �������� ��������� �������� [ x(i-1) < x_new < x(i) ]
		// ����� �������� (�������) [ x(i-1) < x_new ]
		// ������ �������� [ x_new < x(i) ]

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
				R_max_index = i + 1; // ����� R(t) // R(i)
			}
		}

		iter_x = x.begin();
		iter_z = z.begin();

		for (register int i = 0; i < R_max_index; i++) // ����� ����� ����.�������� (���������)
		{
			iter_z++;
			iter_x++;
		}

		// ����������
		curr_right = *iter_x;
		curr_left = *--iter_x;
		++k;
	}

	//end = omp_get_wtime();
	//time = end - start;

	--new_size; // ������������� ����� ������ ���������
	steps = new_size;

	// ��������� ��� ���������������� ���������������� ������
	if (procs == 1)
	{
		solutionZ = *iter_z;
		solutionX = curr_right;
		Arguments_X = x;
		Values_Z = z;
	} // ����� - ��� �������, �������� � ����� ������ ���� ������� � �������������
	else
	{
#pragma omp critical
		{
			short thread = omp_get_thread_num();
			parallel_x[thread] = x;
			parallel_z[thread] = z;
			parallel_steps += new_size;
		}
	}

	//return solutionX;
	return curr_right; // ��� ������������ ������
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
	// �������� ������ ��� �������� ����� ��������� ������ ���������
	parallel_x.resize(procs);
	parallel_z.resize(procs);

	std::vector<block_1D> optimized_values(procs); // left = x, right = z;
	omp_set_dynamic(0);

	double start = omp_get_wtime();
	omp_set_num_threads(procs);

#pragma omp parallel for num_threads(procs)
	for (short i = 0; i < procs; i++)
	{
		optimized_values[i].x_right = Optimize(_left + (_right - _left)*i / procs, _left + (_right - _left)*(i + 1) / procs, _N_max, _Eps);
		optimized_values[i].z_right = Func(optimized_values[i].x_right);
	}

	// ����� �������
	double min_x = optimized_values[0].x_right;
	double min_z = optimized_values[0].z_right;

	double tmp;

#pragma omp parallel sections num_threads(2)
	{
#pragma omp section
		{
			for (short i = 1; i < procs; i++)
			{
				// ������� ������� �� x � �� z
				tmp = optimized_values[i].z_right;
				if (min_z > tmp)
				{
					min_x = optimized_values[i].x_right;
					min_z = optimized_values[i].z_right;
				}
				// ��������� ��� ����� ��������� ���� �������� � ���� �������������� ������
			}
		}
#pragma omp section
		{
			// ������� �� ������ ����������� ������������
			if (!Arguments_X.empty())
			{
				Arguments_X.clear();
				Values_Z.clear();
			}
		}
	}


	// ��������� � ������ ������ ����������� ��� ������ ����������� ������ �������
	for (short i = 0; i < procs; i++)
	{
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
	int k = 0; // ���������� �����
	double eps = _Epsilon;
	double curr_left, curr_right; // ������� �������� ����� � ������ ��� ������ ��������
	double start, end; // ����� ������
	double m_small = 1., M_big;
	double *x_new = new double[procs];


	curr_left = left; // ������������� ������
	curr_right = right; // ������������� ������

	std::vector<block_1D> InitVector(procs);
	std::vector<block_1D> WorkVector(procs * 2);
	std::priority_queue<block_1D> Queue;
	std::list<block_1D> Intervals;

	int size_bank_intervals; // ������ list<block_1D> Intervals

	std::list<block_1D>::iterator Pointer;


	// ������ ��������� ������������� ������������� ���������
	// ������� ������� �� ������ �� ���������� ��
	// ���������� R ������������� ��� ������� - ����������
	double len = right - left; // ���������� ������
	double li = len / procs; // ����� i-�� �������
	register double tmp = 0.;

	omp_set_dynamic(0);

	start = omp_get_wtime();

	omp_set_num_threads(procs);
	{
#pragma omp parallel for firstprivate(tmp) num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			tmp = (right - len) + (i * li); // _b-(b-a)_ + _����������_�������_ * _�����_�������_;
			InitVector[i].x_left = tmp;
			InitVector[i].x_right = tmp + li;
			InitVector[i].z_left = Func(tmp);
			InitVector[i].z_right = Func(tmp + li);
		}

		// ������� �������� � ������� ���������
#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			InitVector[i].M = M(InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_right, InitVector[i].x_left);
		}

		M_big = InitVector[0].M;

		for (register int i = 1; i < procs; i++) // ��� ���� � ���-���� �������������� (������������ ������) // ���������� ����� � ������� ����� CPU_COUNT
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
		// ������ ������������� �������� ���������
		///////////////////////////////////////////
		//////////// ������������ ���� ////////////
		///////////////////////////////////////////
		while (k < _Steps && (abs(curr_left - curr_right) > eps))
		{
			Pointer = --Intervals.end();
			for (register int i = 0; i < procs * 2 - 2; i += 2) // ������� �� ������ CPU_COUNT �������� � ����. ���-��� R � �������� �� �������
			{
				WorkVector[i] = *Pointer;
				Intervals.erase(Pointer--);
			}
			WorkVector[procs * 2 - 2] = *Pointer; Intervals.erase(Pointer);
			size_bank_intervals = Intervals.size();
			// ���� �� ������ ������� ����� ������ ��������� x_new
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

			// �� �����������, ������ �� ��������, ���� ������� ��������������
#pragma omp parallel for num_threads(procs)
			for (int i = 0; i < procs * 2; i++)
			{
				WorkVector[i].M = M(WorkVector[i].z_right, WorkVector[i].z_left, WorkVector[i].x_right, WorkVector[i].x_left);
			}
			register double M_max_array = WorkVector[0].M;
#pragma omp parallel sections num_threads(2)
			{

#pragma omp section
				{
					for (register int i = 1; i < procs * 2; i++) // ��� ���� � ���-���� �������������� (������������ ������) // ���������� ����� � ������� ����� CPU_COUNT * 2
					{
						if (M_max_array < WorkVector[i].M) M_max_array = WorkVector[i].M;
					}
				}
#pragma omp section
				{
					// ������ ���� M_max ����� ���������, ������� � ������������ �������
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

						// ���������� M_max_array � M_max_list, �������� �������
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
			// ���������� R �������������� ��� ������� ��������
			// ���������� ����������� R ��� �������� � �������
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
				// ���������� CPU_COUNT �������� � ������������� ���������������� � ������������ ������� ��� �������� ���������
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

	Values_Z.clear();
	Arguments_X.clear();

	// ���������� ����������� ���������
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

	steps = k;
	// �������
	InitVector.clear();
	WorkVector.clear();
	Intervals.clear();
	delete[] x_new;
	//InitVector.~vector();
	//WorkVector.~vector();
	//Intervals.~list();

	return solutionZ;
}

double Parallel_1D::OptimizeGPU(const double _left, const double _right, const int _N_max, const double _Eps)
{
	return 0.0;
}

double Sequental_2D::Optimize(const double *_Left, const double *_Right, const int _N_max, const double _Eps)
{
	time = 0.L;

	double curr_left_y = _Left[0], curr_left_x = _Left[1];
	double curr_right_y = _Right[0], curr_right_x = _Right[1];

	// ����� ������� |X(i) - X(i-1)|^1/N , N = 2

	// ������� ������
	Values_Z.clear();
	Arguments_X.clear();
	Arguments_Y.clear();

	double start, end; // time
	double curr_left, curr_right; // x[0;1]

	double *A = new double[2]; // ����������� �����
	double *B = new double[2]; // ����������� ������

	A[0] = _Left[0]; // x
	A[1] = _Left[1]; // y
	B[0] = _Right[0]; // x
	B[1] = _Right[1]; // y

	//// ������ ����������
	double temp = 0.; // ���������� ��� �������� 
	int k = 1; // ���������� ����� � ��������� ���������
	int new_size; // ���� ����
	double M_big, m_small = 1.; // ��������� ��� �������� �������� ��������� "�" � ��������� "m"
	int R_max_index = 0; // ���������� ������� ��������� ����� ���������

	////
	//// ���������� ��� �������� �������� � ���������� �������
	std::list<double> z; //�������� Z
	std::list<double> x; //�������� ������� X
	std::list<double> x_images; //�������� X[0;1]
	std::list<double> y; //�������� ������� Y
	std::list<double> M_vector; // ������ ����� ����� M ��� ������ m = max(M[i])
	std::list<double> R_vector; // ������ �������� R ��� �������� ������� ����� ���������� ���������

	block_1D val_z; // ��������� ��� �������� ����� � ������ ������� ������� � Z
	block_1D arg_x; // ��������� ��� �������� ����� � ������ ������� ������� � �
	//block_1D arg_y; // ��������� ��� �������� ����� � ������ ������� ������� � Y

	// ���� ���������� ��� ������ �� ��������
	std::list<double>::iterator place_M, place_R; // ���������-��������� ���� �������� ����� ��������
	std::list<double>::iterator iter_x, iter_y, iter_z, iter_x_image, iter_R, iter_M;
	//

	TEvolvent evolve(2, 20); // ������ ���������
	evolve.SetBounds(A, B); // ��������� �����������

	// ������ ��������� ��� x(1)
	// y[0] = x ; y[1] = y � �������� ��������� ������� z = f(x,y)

	// transformation from hypercube P to hyperinterval D � GetImage
	// transformation from hyperinterval D to hypercube P � GetPreimages � GetInverseImage
	// ������ ��������� � ����� x[0;1] = 0.5; // ����� �� ���������

	double *y_arr = new double[2]; 
	double *x_img = new double; // x[0;1] 

	start = omp_get_wtime();
	// ������ ������� ����
	{
		*x_img = 0.5; 
		x_images.push_back(0);
		x_images.push_back(0.5);
		x_images.push_back(1);
		
		evolve.GetImage(*x_img, y_arr); // �������� g = y(x[0;1])
		
		// ����� ������� ��������� ��������� 
		double borders[2] = { A[0], A[1] };
		x.push_back(A[0]); y.push_back(A[1]); z.push_back(Func(borders));
		x.push_back(y_arr[0]); y.push_back(y_arr[1]); z.push_back(Func(y_arr)); borders[0] = B[0]; borders[1] = B[1];
		x.push_back(B[0]); y.push_back(B[1]); z.push_back(Func(borders));

		iter_z = z.begin();
		iter_x_image = x_images.begin();
		//iter_x = x.begin();
		//iter_y = y.begin();

		// ������� ����� �
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x_image++;
			arg_x.x_right = *iter_x_image;
			M_vector.push_back(M(val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left)); // ������� ����� �_big
		}

		double max = M_vector.front();

		if (max < M_vector.back()) max = M_vector.back();
		M_big = max; // ����� ������������ M_big

		m_small = (M_big == 0) ? 1 : r * M_big;

		R_max_index = 1;

		iter_z = z.begin();
		iter_x_image = x_images.begin();
		//iter_y = y.begin();
		//iter_x = x.begin();
		// ������� � ��������� ����� �������������� ��� �������� ��������� � ������������� ��� ������������� "�������" ���������
		// ������ �������� ��������� �������� x(i-1)<x_new<x(i)
		// ����� �������� (�������) x(i-1)<x_new
		// ������ �������� x_new<x(i)
		for (int i = 1; i < 3; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x_image++;
			arg_x.x_right = *iter_x_image;

			temp = R(m_small, val_z.z_right, val_z.z_left, arg_x.x_right, arg_x.x_left);
			R_vector.push_back(temp);
			if (i == 1) max = temp;

			if (max < temp)
			{
				max = temp;
				R_max_index = i;
			}
		}

		//iter_z = z.begin();
		//iter_x = x.begin();
		//iter_y = y.begin();
		iter_x_image = x_images.begin();

		for (int i = 0; i < R_max_index; i++)
		{
			//iter_z++;
			//iter_x++;
			iter_x_image++;
			//iter_y++;
		}

		curr_right = *iter_x_image;
		curr_left = *--iter_x_image;
	}
	//// ������ ������ ��� � ���������� � ����� ����� �� ������� //
	while (k < _N_max && (abs(curr_left - curr_right) > _Eps))
	{
		new_size = 2 + k; // + 2 ��� �������� ���������, ����� ����� �������������

		iter_z = z.begin(); // ����� ��������� � ������������ ���������������
		iter_x = x.begin();
		iter_x_image = x_images.begin();
		iter_y = y.begin();
		// ���� ��������� ������� ���������
		for (register int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
			iter_y++;
			iter_x_image++;
		}
		// �������� �������� ��� �������� � ������ ������ �����
		val_z.z_right = *iter_z--;
		val_z.z_left = *iter_z++;

		//left = curr_left;
		//right = curr_right;
		// ���������� ����� ����� ��������� x_new
		*x_img = (curr_right + curr_left) / 2 - (val_z.z_right - val_z.z_left) / (2 * m_small);

		evolve.GetImage(*x_img, y_arr);

		// ��������� x_new, y_new � z_new
		z.insert(iter_z, Func(y_arr));
		x_images.insert(iter_x_image, *x_img);
		x.insert(iter_x, y_arr[0]);
		y.insert(iter_y, y_arr[1]);

		iter_z = z.begin();
		//iter_x = x.begin();
		iter_x_image = x_images.begin();
		//iter_y = y.begin();
		place_M = M_vector.begin();
		place_R = R_vector.begin();

		// ���� �������� � ������ ��������������� ������� ����������
		for (int i = 0; i < R_max_index - 1; i++)
		{
			iter_z++;
			//iter_x++;
			iter_x_image++;
			//iter_y++;
			if (i == R_max_index - 1) break;
			place_M++;
			place_R++;
		}
		// ������ �������� ��������� �������� x(i-1)<x_new<x(i)
		// ����� �������� (�������) x(i-1)<x_new
		// ������ �������� x_new<x(i)
		for (int i = 0; i < 2; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x_image++;
			arg_x.x_right = *iter_x_image;
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

		// ������ �������� ��������� �������� [ x(i-1) < x_new < x(i) ]
		// ����� �������� (�������) [ x(i-1) < x_new ]
		// ������ �������� [ x_new < x(i) ]

		iter_z = z.begin(); 
		//iter_x = x.begin(); 
		//iter_y = y.begin(); 
		iter_x_image = x_images.begin();

		R_max_index = 1;
		place_R = R_vector.begin();

		for (int i = 0; i < new_size; i++)
		{
			val_z.z_left = *iter_z++;
			val_z.z_right = *iter_z;
			arg_x.x_left = *iter_x_image++;
			arg_x.x_right = *iter_x_image;
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
				R_max_index = i + 1; // ����� R(t) // R(i)
			}
		}

		iter_x = x.begin();
		iter_x_image = x_images.begin();
		iter_y = y.begin();
		iter_z = z.begin();

		for (register int i = 0; i < R_max_index; i++) // ����� ����� ����.�������� (���������)
		{
			iter_z++;
			iter_x++;
			iter_y++;
			iter_x_image++;
		}

		// ����������
		curr_right_x = *iter_x--;
		curr_right_y = *iter_y--;
		curr_left_x = *iter_x;
		curr_left_y = *iter_y;

		curr_right = *iter_x_image;
		curr_left = *--iter_x_image;
		++k;
	}
	end = omp_get_wtime();
	time = end - start;

	--new_size; // ������������� ����� ������ ���������
	steps = new_size;

	// ��������� ���������
	solutionX = curr_right_x;
	solutionY = curr_right_y;
	solutionZ = Func(y_arr);

	// �� ����������� ������!!!
	// ����� ����� � ��������, ���� ������������ zi = f(x_real_i, y_real_i) � x_img_i[0;1]

	Values_Z = z;
	Arguments_X = x;
	Arguments_Y = y;

	// �������
	delete[]A;
	delete[]B;
	delete[]y_arr;
	delete x_img;
	z.clear();
	x.clear();
	x_images.clear();
	y.clear();
	R_vector.clear();
	M_vector.clear();

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

	// ����� ������� |X(i) - X(i-1)|^1/N , N = 2

	// ������� ������
	Values_Z.clear();
	Arguments_X.clear();
	Arguments_Y.clear();

	double start, end; // time
	double curr_left, curr_right; // x[0;1]

	double *A = new double[2]; // ����������� �����
	double *B = new double[2]; // ����������� ������

	A[0] = Left[0]; // x
	A[1] = Left[1]; // y
	B[0] = Right[0]; // x
	B[1] = Right[1]; // y


	int k = 0; // ���������� �����
	double eps = _Epsilon;
	double m_small = 1., M_big;
	double *x_img = new double[procs];
	double *y_arr = new double[2];

	std::vector<block_2D> InitVector(procs);
	std::vector<block_2D> WorkVector(procs * 2);
	std::list<block_2D> Intervals;

	int size_bank_intervals; // ������ list<block_2D> Intervals

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

	// ������ ��������� ������������� ������������� ���������
	// ������� ������� �� ������ �� ���������� ��
	// ���������� R ������������� ��� ������� - ����������
	double len_x = right_x - left_x; // ���������� ������ x
	double li_x = len_x / procs; // ����� i-�� ������� x

	double len_x_img = right_x_img - left_x_img; // ���������� ������ x
	double li_x_img = len_x_img / procs; // ����� i-�� ������� x

	double len_y = right_y - left_y; // ���������� ������ y
	double li_y = len_y / procs; // ����� i-�� ������� y
	double tmp[2];
	double tmp_img = 0;

	start = omp_get_wtime();

	omp_set_num_threads(procs);
	{
#pragma omp parallel for firstprivate(tmp, tmp_img) num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			tmp[0] = (right_x - len_x) + (i * li_x); // _b-(b-a)_ + _����������_�������_ * _�����_�������_;
			tmp[1] = (right_y - len_y) + (i * li_y); // _b-(b-a)_ + _����������_�������_ * _�����_�������_;
			tmp_img = (right_x_img - len_x_img) + (i * li_x_img); // _b-(b-a)_ + _����������_�������_ * _�����_�������_;

			InitVector[i].x_left = tmp[0];
			InitVector[i].x_right = tmp[0] + li_x;

			InitVector[i].x_img_left = tmp_img;
			InitVector[i].x_img_right = tmp_img + li_x_img;

			InitVector[i].y_left = tmp[1];
			InitVector[i].y_right = tmp[1] + li_y;

			evolve[i]->GetImage(tmp_img, tmp);

			InitVector[i].z_left = Func(tmp);

			evolve[i]->GetImage(tmp_img + li_x_img, tmp);

			InitVector[i].z_right = Func(tmp);
		}

		// ������� �������� � ������� ���������
#pragma omp parallel for num_threads(procs)
		for (int i = 0; i < procs; i++)
		{
			InitVector[i].M = M(InitVector[i].z_right, InitVector[i].z_left, InitVector[i].x_img_right, InitVector[i].x_img_left);
		}

		M_big = InitVector[0].M;

		for (register int i = 1; i < procs; i++) // ��� ���� � ���-���� �������������� (������������ ������) // ���������� ����� � ������� ����� CPU_COUNT
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
		// ������ ������������� �������� ���������
		///////////////////////////////////////////
		//////////// ������������ ���� ////////////
		///////////////////////////////////////////
		while (k < _Steps && (abs(curr_left - curr_right) > eps))
		{
			Pointer = --Intervals.end();
			for (register int i = 0; i < procs * 2 - 2; i += 2) // ������� �� ������ CPU_COUNT �������� � ����. ���-��� R � �������� �� �������
			{
				WorkVector[i] = *Pointer;
				Intervals.erase(Pointer--);
			}
			WorkVector[procs * 2 - 2] = *Pointer; Intervals.erase(Pointer);
			size_bank_intervals = Intervals.size();
			// ���� �� ������ ������� ����� ������ ��������� x_new
#pragma omp parallel for shared(x_img)  num_threads(procs)
			for (int i = 0; i < procs; i++)
			{
				register double y_arr[2];
				x_img[i] = 0.5*(WorkVector[i * 2].x_img_right + WorkVector[i * 2].x_img_left)
					- (WorkVector[i * 2].z_right - WorkVector[i * 2].z_left) / (2 * m_small);
				evolve[i]->GetImage(x_img[i], y_arr);
				WorkVector[i * 2 + 1].x_left = y_arr[0];
				WorkVector[i * 2 + 1].y_left = y_arr[1];
				WorkVector[i * 2 + 1].z_left = Func(y_arr);
			}

			for (register int i = 0; i < procs * 2; i += 2)
			{
				WorkVector[i + 1].x_right = WorkVector[i].x_right;
				WorkVector[i + 1].y_right = WorkVector[i].y_right;
				WorkVector[i + 1].x_img_right = WorkVector[i].x_img_right;
				WorkVector[i + 1].z_right = WorkVector[i].z_right;
				// ������� ��������
				WorkVector[i + 1].x_img_left = WorkVector[i].x_img_right = x_img[i / 2];

				WorkVector[i].x_right = WorkVector[i + 1].x_left;
				WorkVector[i].y_right = WorkVector[i + 1].y_left;
				WorkVector[i].z_right = WorkVector[i + 1].z_left;
			}

			// �� �����������, ������ �� ��������, ���� ������� ��������������
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
					for (register int i = 1; i < procs * 2; i++) // ��� ���� � ���-���� �������������� (������������ ������) // ���������� ����� � ������� ����� CPU_COUNT * 2
					{
						if (M_max_array < WorkVector[i].M) M_max_array = WorkVector[i].M;
					}
				}
#pragma omp section
				{
					// ������ ���� M_max ����� ���������, ������� � ������������ �������
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

						// ���������� M_max_array � M_max_list, �������� �������
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
			// ���������� R �������������� ��� ������� ��������
			// ���������� ����������� R ��� �������� � �������
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
				// ���������� CPU_COUNT �������� � ������������� ���������������� � ������������ ������� ��� �������� ���������
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

	// ���������� ����������� ���������
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
		Arguments_Y.push_back(iter->y_left);
		Arguments_Y.push_back(iter->y_right);
	}
	Values_Z.unique(); Arguments_X.unique(); Arguments_Y.unique();

	steps = k;
	// �������

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
