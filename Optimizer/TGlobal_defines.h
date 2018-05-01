#ifndef __TGLOBAL_DEFINES__
#define __TGLOBAL_DEFINES__
#pragma warning (disable: 6255)
#pragma warning (disable: 6386)
#pragma warning (disable: 6385)
//#define USE_CL_200 // определить использование OpenCL 2.0
#define USE_CL_100 // определить использование OpenCL 1.2

#define  _CRT_SECURE_NO_WARNINGS
#pragma warning (disable: 4996) // Обойти ошибку использования устаревших функций
								// Конкретно для OpenCL 1.2

#if defined(USE_CL_200)  // Если используется CL 2.0

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_ENABLE_EXCEPTIONS
#include <CL\cl2.hpp>

#elif defined(USE_CL_100)  // Если используется CL 1.2

#define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL\cl.hpp>
#else  // Иначе исключить использование OpenCL

#define NOT_USE_CL

#endif 
#pragma warning (default: 6385)
#pragma warning (default: 6386)
#pragma warning (default: 6255)

#define NO_MPREAL // запрет на использование в Extended и Evolvent библиотеки MPFREAL

#define scanf scanf_s



// *********************** ЗАМЕТКИ // *********************** //
/*1) Глобалайзер не работает при NO_MPREAL - ломатется стек вокруг развертки TEvolvent при разрушении объекта */
/*2) Необходимо заменить все Extended в развертках в виде приходящих аргументов -> 
  -> Изначальные строки закомментить, при реализации развертки для N > 2 пространства - переделать код
  -> Поиск по коду, где была замена по слову "ЗАМ" в комментах*/
/*3) --*/
 

// ********************************************************** //


#endif
