#ifndef __DEFINES_H__
#define __DEFINES_H__


/* ======================================================================== *\
**  Constants                                                               **
\* ======================================================================== */

/// Максимальная размерность
#define MAX_DIM 50
// Максимальное количество глобальных минимумов
#define MAX_NUM_MIN 20
// Если константа определена то на GPU будет использоваться double иначе float
#define CUDA_VALUE_DOUBLE_PRECISION
/// Если константа определена то на СPU будет использоваться double иначе float

#define CPU_VALUE_DOUBLE_PRECISION

/**
Данный параметр определяет уровень, меньше которого константа Липшица считается нулевой.

Проблема возникновения близкой к нулю оценки константы может возникнуть, если функция принимает
одинаковые значения во многих точках (например, функция Растригина).
*/
#define _M_ZERO_LEVEL 1e-12


/* ======================================================================== *\
**  Types                                                                   **
\* ======================================================================== */
/// Используемый тип данных для вычислений на GPU
#ifdef CUDA_VALUE_DOUBLE_PRECISION
#define CUDA_VALUE double
#else
#define CUDA_VALUE float
#endif
/// Используемый тип данных для вычислений на центральном процессоре(не доделано),
#ifdef CPU_VALUE_DOUBLE_PRECISION
#define CPU_VALUE double
#else
#define CPU_VALUE float
#endif
/// Массив целых чисел
#define ints int*
/// Массив действительных чисел
#define doubles double*
/// Тип флага
#define FLAG bool

/* ======================================================================== *\
**  Defines                                                                 **
\* ======================================================================== */

/* ======================================================================== *\
**  Параметры спецификаторы и прочее,  для разных версий компиляторов       **
\* ======================================================================== */

#ifdef _GPU_CUDA_ ///Для компиляции под CUDA

#define perm 1
#define SPECIFIER __shared__
#define ARRAY_SPECIFIER __constant__
#define concatenation cuda
#define F_DEVICE __device__
#define parameter_const
#define OBJECTIV_TYPE CUDA_VALUE
#define GET_FUNCTION_PARAMETERS OBJECTIV_TYPE* x, OBJECTIV_TYPE* f
#define FUNCTION_CALCOLATION_PREF (x)

#define GKLS_VARIABLES_SPECIFIER __shared__
//Формирование констант правильной точности
#ifdef CUDA_VALUE_DOUBLE_PRECISION
#define PRECISION(x) x##0
#else
#define PRECISION(x) x##f
#endif

#else ///Для компиляции под VS

#define perm 0
#define SPECIFIER extern
#define ARRAY_SPECIFIER extern
#define concatenation
#define F_DEVICE inline
#define parameter_const const
#define OBJECTIV_TYPE CPU_VALUE
#define GET_FUNCTION_PARAMETERS tFunction* f
#define FUNCTION_CALCOLATION_PREF

#define GKLS_VARIABLES_SPECIFIER extern
//Формирование констант правильной точности
#ifdef CPU_VALUE_DOUBLE_PRECISION
#define PRECISION(x) x##0
#else
#define PRECISION(x) x##f
#endif

#endif

/* ======================================================================== *\
**  Служебные                                                               **
\* ======================================================================== */
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? a : b)
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? a : b)
#endif

/// Объединить два слова
#define CAT(x, y) x##y
/// Объединить четыре слова
#define CAT4(a, b, c, d) a##b##c##d
/// Объединяет два слова в обратномпорядке
#define CONCATENATION2(name, console) CAT(console,name)
/** Добавляет префикс concatenation
(см. Параметры спецификаторы и прочее,  для разных версий компиляторов  )
*/
#define CONCATENATION(name) CONCATENATION2(name, concatenation)
/// Превращает слово в комментарий
#define CAT_COM(x) //##x

/// Добавляет префикс P к имени типа, что соответствует перечислению EParameterType
#define ParType(type) P##type
/// Добавляет префикс link к имени
#define LinkParameter(name) link##name
/// Добавляет префикс com к имени
#define ComParameter(name) com##name
/// Добавляет префикс help к имени
#define HelpParameter(name) help##name
/// Добавляет префикс inc к имени
#define IncParameter(name) inc##name
/// Переводит слово bar в строку типа char*
#define make_str(bar) # bar
/// Добавляет префикс IS_ к имени
#define IsChange(name) IS_##name

/// Инициализирует параметр из класса TParameters
#define InitParameter(type, name, defVal, com, help, sizeVal)                            \
IsChange(name) = false;                                                                  \
Inc(ParType(name), ParType(type), LinkParameter(name), com,                              \
  HelpParameter(name), help, (void*)(&name), make_str(defVal),                           \
  make_str(name), &IsChange(name), IncParameter(name), make_str(name), sizeVal);

#define OWNER_NAME Owner

/// Инициализирует параметр из класса Параметров
#define InitOption(name, defVal, com, help, sizeVal)                           \
  InitializationOption((TBaseProperty<OWNER_NAME>*)(&name), make_str(name), make_str(defVal), com, help, sizeVal);

/* ======================================================================== *\
**  Директиы для объявления переменных                                      **
\* ======================================================================== */

/** простое объявление переменной
name - имя переменной
type - тип переменной
specifier - спецификатор
*/
#define VARIABLES(name, type, specifier) specifier type name
/** Объявление переменной со спецификатором GKLS_VARIABLES_SPECIFIER
(см. Параметры спецификаторы и прочее,  для разных версий компиляторов  )
*/
#define GKLS_VARIABLES(name, type) VARIABLES(name, type, GKLS_VARIABLES_SPECIFIER)
/** Объявление переменной со спецификатором SPECIFIER
(см. Параметры спецификаторы и прочее,  для разных версий компиляторов  )
*/
#define CONSTANT_VARIABLES(name, type) SPECIFIER type name
/// Объявляет функцию для инициализации значения по умолчанию
#define NEW_FUNC_DEF(name) CONCATENATION(name##_func_def())
/// Создает статический массив
#define STATIC_ARRAY(type, name, count) type name[count]
/** Создает статический массив со спецификатором ARRAY_SPECIFIER
(см. Параметры спецификаторы и прочее,  для разных версий компиляторов  )
*/
#define NEW_ARRAY(type, name, count) ARRAY_SPECIFIER STATIC_ARRAY(type, name, count)
/** Создает статический массив размером MAX_DIM типа OBJECTIV_TYPE со спецификатором ARRAY_SPECIFIER
(см. Параметры спецификаторы и прочее,  для разных версий компиляторов  )
*/
#define NEW_ARRAY_MAX_SIZE(name) NEW_ARRAY(OBJECTIV_TYPE, name, MAX_DIM)
/** Создает переменнут типа OBJECTIV_TYPE со спецификатором SPECIFIER
(см.Параметры спецификаторы и прочее, для разных версий компиляторов)
*/
#define FLOAT_VARIABLES(name) SPECIFIER OBJECTIV_TYPE name


/**Глобальная переменная объявляемая в файле с задачей, могут быть переопеределены в параметрах,
!!!В problems.cpp требуется повторное объявление: type name;!!!
type - тип
name - имя
val - значение по умолчанию
*/
#define NEW_VARIABLES(type, name, val) CONSTANT_VARIABLES(name, type);\
 F_DEVICE type NEW_FUNC_DEF(name)\
  {\
  name = val;\
return name;\
}\

/// Объявляет переменные для создания параметра класса TParameters
#define CreateParameter(type, name)       \
public: type name;                        \
public: FLAG IsChange(name);              \
protected: int IncParameter(name);        \
protected: EParameterType ParType(name);  \
protected: string LinkParameter(name);    \
protected: string HelpParameter(name);

/// Базовые переопределения для классов типов данные
#define BasicMethods(ClassType, Type)                                             \
  virtual void operator =(Type data) {TTypedProperty<Type, Owner>::operator=(data);}                \
  virtual void operator=(ClassType<Owner>& data) {TParameterProperty<Type, Owner>::operator=(data);}     \
  virtual void operator = (ClassType<Owner>* data) {TParameterProperty<Type, Owner>::operator=(data);}  \
  virtual string ToString() {return operator string();}                           \
  virtual void FromString(string val) {operator = (val);}                           \
  virtual void operator = (TBaseProperty<Owner>& data) { TParameterProperty<Type, Owner>::operator=((TParameterProperty<Type, Owner>*)(& data)); } \
  virtual void operator=(char* data) {operator = (string(data));}

/* ======================================================================== *\
**  Директиы для объявления функций                                         **
\* ======================================================================== */

/** Создает функцию с имененм concatenation##name, возвращаемое значение типа type
(значение concatenation см. Параметры спецификаторы и прочее,  для разных версий компиляторов  )
со спецификатором console,
Параметром parameter типа массив type и модификатором parameterconst
*/
#define OBJECTIV_FUNCTION2(name, type, parameter, parameterconst, console) console type CONCATENATION(name)(parameterconst type parameter[])
/** Создает функцию с имененм concatenation##name, возвращаемое значение типа type
со спецификатором F_DEVICE,
Параметром parameter типа массив type и модификатором parameter_const
(см.Параметры спецификаторы и прочее, для разных версий компиляторов)
*/
#define OBJECTIV_FUNCTION(name, type, parameter) OBJECTIV_FUNCTION2(name, type, parameter, parameter_const, F_DEVICE)
/** Создает функцию с имененм concatenation##name, возвращаемое значение типа type
со спецификатором console,
Двумя параметрами parameter1, parameter2 типа массив type и модификатором parameterconst
*/
#define P2_FUNCTION2(name, type, parameter1, parameter2, parameterconst, console) console type CONCATENATION(name)(parameterconst type* parameter1, parameterconst type* parameter2)
/** Создает функцию с имененм concatenation##namen, возвращаемое значение типа type
со спецификатором F_DEVICE,
Двумя параметрами parameter1, parameter2 типа массив type и модификатором parameter_const
*/
#define P2_FUNCTION(name, type, parameter1, parameter2) P2_FUNCTION2(name, type, parameter1, parameter2, parameter_const, F_DEVICE)
//////////////////////////////////////////////////////////////////////////

/** Функция  с имененм concatenation##name
возвращает значение типа OBJECTIV_TYPE, со спецификатором F_DEVICE,
parameter - имя параметра, типа OBJECTIV_TYPE*  модификатором parameterconst
(см.Параметры спецификаторы и прочее, для разных версий компиляторов)
*/
#define HYBRID_OBJECTIV_FUNCTION(name, parameter) OBJECTIV_FUNCTION(name, OBJECTIV_TYPE, parameter)
/** Функция  с имененм concatenation##objfn
возвращает значение типа OBJECTIV_TYPE, со спецификатором F_DEVICE,
С параметром x , типа OBJECTIV_TYPE*  модификатором parameterconst
(см.Параметры спецификаторы и прочее, для разных версий компиляторов)
*/
#define OBJECTIV HYBRID_OBJECTIV_FUNCTION(objfn, x)
/** Функция ограничения, меет имя concatenation##constraint##num,
возвращает значение типа OBJECTIV_TYPE и со спецификатором F_DEVICE,
С параметром x , типа OBJECTIV_TYPE*  модификатором parameterconst
(см.Параметры спецификаторы и прочее, для разных версий компиляторов)
num - номер функции
*/
#define CONSTRAINT(num) HYBRID_OBJECTIV_FUNCTION(CAT(constraint,num), x)

/** Функция от дву параметров возвращает значение типа OBJECTIV_TYPE,
name - имя
parameter1 и parameter2 - имена параметров, типа OBJECTIV_TYPE*
возвращает значение типа OBJECTIV_TYPE
*/
#define HYBRID_P2_FUNCTION2(name, parameter1, parameter2) P2_FUNCTION(name, OBJECTIV_TYPE, parameter1, parameter2)
/** Изменяет функцию для использования в CUDA,
Принимает:
type - тип возвращаемого значения
name - имя функции
*/
#define NEW_FUNCTION(type, name) F_DEVICE type CONCATENATION(name)
/// Задает список функций задачи
#define GET_FUNCTION NEW_FUNCTION(void, GetFunction)(GET_FUNCTION_PARAMETERS)
/// целевая функция
#define OBJ CONCATENATION(objfn) FUNCTION_CALCOLATION_PREF
/// Функция ограничения
#define CON(num) CONCATENATION(CAT(constraint,num)) FUNCTION_CALCOLATION_PREF

/**Объявление Set и Get функции c именем
Get##N возвращает тип T
и Set##N параметр value типа T
*/
#define PROPERTY(T, N)     \
  T Get ## N() const;     \
  void Set ## N(T value);

/* ======================================================================== *\
**  Прочие константы                                                        **
\* ======================================================================== */

#define PI PRECISION(3.14159265359)







#endif
