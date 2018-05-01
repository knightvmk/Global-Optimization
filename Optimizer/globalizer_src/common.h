/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      common.h                                                    //
//                                                                         //
//  Purpose:   Common Header file                                          //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __COMMON_H__
#define __COMMON_H__

#include "defines.h"

/* ============================================================================================= *\
**  Constants                                                                                    **
\* ============================================================================================= */

//const int MaxPathLength = 512
const int MaxNumOfTaskLevels = 5;
const int MaxNumOfFunc = 20;
const int MaxDim = MAX_DIM;
const int MaxNumOfGlobMinima = MAX_NUM_MIN;
const int MaxM = 20;
const int MaxL = 10;

const double MaxDouble = 1.7e308;
const double MinDouble = -1.7e308;
const double AccuracyDouble = 0.00000001;

const int DefaultQueueSize = 32767;//1023;//65535;//32767;//16383;//16777215;//8388607;//1048575;//8388607;//524287; //262143; // должно быть равно 2^k - 1
const int DefaultSearchDataSize = 100000;

const int TagChildStartSolve = 101;
const int TagChildSolved = 102;
const int ChildStopMsg = -101;

/* ============================================================================================= *\
**  Types                                                                                        **
\* ============================================================================================= */
class TProcess;
#ifdef WIN32
typedef void(__cdecl *tIterationHandler)(TProcess *pProcess);
typedef OBJECTIV_TYPE(__cdecl *tFunction)(const OBJECTIV_TYPE*);
#else
typedef void(*tIterationHandler)(TProcess *pProcess);
typedef OBJECTIV_TYPE(*tFunction)(const OBJECTIV_TYPE*);
#endif


/**
Тип итерации метода

При работе смешанного алгоритм в зависимости от параметра смешивания итерации могут выполняться
либо с использованием глобальных характеристик интервала, либо - локальных.
*/
enum IterationType
{
  ///Итерация выполняется по правилам глобального поиска
  Global,
  ///Итерация выполняется по правилам локального поиска
  Local
};


enum EParameterType
{
  Pbool,
  Pint,
  Pdouble,
  Pstring,
  PETypeMethod,
  PETypeCalculation,
  PELocalMethodScheme,
  PESeparableMethodType,
  PEStopCondition,
  PETypeProcess,
  PETypeSolver,
  PEMapType,
  Pints,
  Pdoubles,
  Pflag
};

enum ETypeMethod
{
  StandartMethod,
  HybridMethod,
  ManyNumPointMethod,
  UniformSearchMethod,
  GlobalSearchMethod,
  MCO_Method,
  MultievolventsMethod,
  ParallelMultievolventsMethod,
  IntegerMethod
};

enum ESeparableMethodType
{
  Off,
  GridSearch,
  GlobalMethod
};

enum ELocalMethodScheme
{
  None,
  FinalStart,
  IntegratedOnePoint,
  IntegratedManyPoints
};

enum EStopCondition
{
  Accuracy,
  OptimumVicinity,
  OptimumVicinity2,
  OptimumValue,
  AccuracyWithCheck
};

enum ETypeCalculation
{
  OMP,
  CUDA,
  PHI
};

enum ETypeProcess
{
  SynchronousProcess,
  SynchronousProcessNew,
  AsynchronousProcess
};

enum ETypeSolver
{
  SeparableSearch,
  SequentialSearch,
  SingleSearch
};

///Тип развертки
enum EMapType
{
  ///Вращаемая развертка
  mpRotated,
  ///Сдвиговая развертка
  mpShifted,
  ///Базовый, одиночный вариант
  mpBase
};

#endif
// - end of file ----------------------------------------------------------------------------------
