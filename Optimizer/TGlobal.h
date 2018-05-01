#ifndef __TGLOBAL_H__
#define __TGLOBAL_H__

#include "TGlobal_defines.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <list>

#include <omp.h>
//#include <mpi.h>
#include "evolvent.h"
#include <iostream>
//#include <CL\opencl.h>

struct block_1D
{
	double R; // ������ ��������� �������
	double x_left; // ���������� �����
	double x_right; // ���������� ������
	double z_left; // �������� �����
	double z_right; // �������� ������
	double M; // �������� ������ ��� ���������
};

struct block_2D
{
	double R; // ������ ��������� �������
	double x_img_left;
	double x_img_right;
	double x_left; // ���������� �����
	double x_right; // ���������� ������
	double y_left; // ���������� �����
	double y_right; // ���������� ������
	double z_left; // �������� �����
	double z_right; // �������� ������
	double M; // �������� ������ ��� ���������
};

bool operator<(const block_1D &i1, const block_1D &i2);
bool operator>(const block_1D &i1, const block_1D &i2);

bool operator<(const block_2D &i1, const block_2D &i2);
bool operator>(const block_2D &i1, const block_2D &i2);

class TGlobal
{
protected:
#ifndef NOT_USE_CL
	std::vector<cl::Platform> avaible_platforms;
	std::vector<std::vector<cl::Device>> avaible_devices;
#endif
	std::vector<std::string> avaible_platform_names;
	std::vector<std::string> avaible_devices_names;
	bool can_use_cl = true;
	void InitCL();
public:
	TGlobal() 
	{
#ifdef NOT_USE_CL
		can_use_cl = false;
#else
		can_use_cl = true;
#endif
		InitCL();
	}
	~TGlobal() 
	{
#ifndef NOT_USE_CL
		avaible_platforms.clear(); avaible_platforms.~vector();
		avaible_devices.clear(); avaible_devices.~vector();
#endif
		avaible_platform_names.clear(); avaible_platform_names.~vector();
		avaible_devices_names.clear(); avaible_devices_names.~vector();
	}
};

class OneDimension : public TGlobal
{

protected:
	double left; // ����� �������
	double right; // ������ �������
	double r; // ��������
	double solutionX = 0;
	double solutionZ = 0;
	int procs = 1;
	long double time;
	__int64 steps;

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
	inline double Func(const double &_x);
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

	double Optimize(const double _left, const double _right, const int _N_max, const double _Eps, const int _threads); //����������� ������� �� �������
	double Optimize(const double _Epsilon, const int _Steps, const int _thread_count); // ����������� �� ��������������� R �������� � R_max �����������
	double OptimizeGPU(const double _left, const double _right, const int _N_max, const double _Eps); // ����������� �� OpenCL
};

class TwoDimension : public TGlobal
{
protected:
	double* Left; // ����������� ������� �����
	double* Right; // ����������� ������� ������
	double r; // ��������
	double solutionX = 0;
	double solutionY = 0;
	double solutionZ = 0;
	int procs = 1;
	long double time;
	__int64 steps;
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
	inline double Func(const double *_y);
	
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
};

class Parallel_2D : public TwoDimension
{
private:
	int parallel_steps = 0;
	std::vector<std::list<double>> parallel_z;
	std::vector<std::list<double>> parallel_x;
	std::vector<std::list<double>> parallel_y;
public:
	std::list<double> Values_Z;
	std::list<double> Arguments_Y;
	std::list<double> Arguments_X;
	Parallel_2D(const double *_Left, const double *_Right, const double _r) : TwoDimension(_Left, _Right, _r) {};
	~Parallel_2D() { parallel_z.clear(); parallel_x.clear(); parallel_y.clear(); Values_Z.clear(); Arguments_X.clear(); Arguments_Y.clear(); }

	double Optimize(const double _Epsilon, const int _Steps, const int _thread_count); // ����������� �� ��������������� R �������� � R_max �����������
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