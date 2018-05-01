#ifndef __TGLOBAL_DEFINES__
#define __TGLOBAL_DEFINES__
#pragma warning (disable: 6255)
#pragma warning (disable: 6386)
#pragma warning (disable: 6385)
//#define USE_CL_200 // ���������� ������������� OpenCL 2.0
#define USE_CL_100 // ���������� ������������� OpenCL 1.2

#define  _CRT_SECURE_NO_WARNINGS
#pragma warning (disable: 4996) // ������ ������ ������������� ���������� �������
								// ��������� ��� OpenCL 1.2

#if defined(USE_CL_200)  // ���� ������������ CL 2.0

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL\cl2.hpp>

#elif defined(USE_CL_100)  // ���� ������������ CL 1.2

#define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL\cl.hpp>
#else  // ����� ��������� ������������� OpenCL

#define NOT_USE_CL

#endif 
#pragma warning (default: 6385)
#pragma warning (default: 6386)
#pragma warning (default: 6255)

#define NO_MPREAL // ������ �� ������������� � Extended � Evolvent ���������� MPFREAL

#define scanf scanf_s



// *********************** ������� // *********************** //
/*1) ����������� �� �������� ��� NO_MPREAL - ��������� ���� ������ ��������� TEvolvent ��� ���������� ������� */
/*2) ���������� �������� ��� Extended � ���������� � ���� ���������� ���������� -> 
  -> ����������� ������ ������������, ��� ���������� ��������� ��� N > 2 ������������ - ���������� ���
  -> ����� �� ����, ��� ���� ������ �� ����� "���" � ���������*/
/*3) --*/
 

// ********************************************************** //


#endif
