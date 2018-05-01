/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      extended.h                                                  //
//                                                                         //
//  Purpose:   Class Extended is used to represent values with high        //
//             precision                                                   //
//                                                                         //
//  Author(s): Chigarev V.                                                 //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#ifndef NO_MPREAL
#include "mpreal.h"
#else
#include <math.h>
#endif
#include "extended.h"

//#define DEBUG_EXTENDED // defined in project properties for Debug configuration

#ifndef NO_MPREAL
ExtendedTypeID Extended::type = etMPFR; // etDouble, etMPFR
#else
ExtendedTypeID Extended::type = etDouble; // etDouble, etMPFR
#endif
double Extended::maxNumbers = 101;
int Extended::digitsNumber = 2;

int Extended::mpfrPrecision;
int Extended::mapmPrecision;
double Extended::precision = 0.01;
float Extended::precisionReserveFactor = 1.5;

int Extended::packedSize = -1;
bool Extended::precisionSet = false;

void Extended::Init()
{
#ifndef NO_MPREAL
  mpfrValue = NULL;
#endif
  if (!precisionSet)
    SetPrecision(0.01);
}

void Extended::Free()
{
#ifndef NO_MPREAL
  if (mpfrValue != NULL)
  {
    delete mpfrValue;
    mpfrValue = NULL;
  }
#endif
#ifdef DEBUG_EXTENDED
  dValue = 0.0;
#endif
}

Extended::Extended(double val)
{
  Init();
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    mpfrValue = new mpreal(val);
#ifdef DEBUG_EXTENDED
    dValue = val;
#endif
    break;
#endif
  case etDouble:
    dValue = val;
  }
}

Extended::Extended(const char* str)
{
  Init();
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    mpfrValue = new mpreal(str);
#ifdef DEBUG_EXTENDED
    dValue = mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    dValue = atof(str);
    break;
  }
}

Extended::Extended(const Extended& e)
{
  Init();
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    mpfrValue = new mpreal(*e.mpfrValue);
#ifdef DEBUG_EXTENDED
    dValue = mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    dValue = e.dValue;
  }
}

Extended::~Extended()
{
  Free();
}

float Extended::getPrecisionReserveFactor()
{
  return precisionReserveFactor;
}

double Extended::GetPrecision()
{
  return precision;
}

void Extended::SetPrecision(double prec)
{
  maxNumbers = ::ceil(1.0/prec+1);
  precisionSet = true;
  precision = prec;
  // transform human like precision to precision of mpfr and mapm
  mpfrPrecision = (int)(::ceil(precisionReserveFactor*log(maxNumbers)/log(2.0)));
  if (type == etDouble)
    digitsNumber = 20;
  else
    digitsNumber = (int)::ceil(::log10(maxNumbers));
  // mpfr precision is amount of bits, used to store numbers
  // the biggest number - 2 ^ (number of bits)
#ifndef NO_MPREAL
  mpreal::set_default_prec(mpfrPrecision);
#endif
  packedSize = -1;
}

int Extended::GetPackedSize()
{
  if (packedSize <= 0)
  {
    switch (type)
    {
#ifndef NO_MPREAL
  case etMPFR:
    // size = size of struct + precision rounded to 4 bytes (sizeof(mp_limb_t) = 4)
    packedSize = sizeof(__mpfr_struct) + ((mpfrPrecision-1)/sizeof(mp_limb_t)/8+1)*sizeof(mp_limb_t);
    break;
#endif
  case etDouble:
    packedSize = sizeof(double);
    break;
    }
  }
  return packedSize;
}

void Extended::Pack(void* result) const
{
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
  {
    char* offset = (char*)result;
    memcpy(offset, (mpfr_ptr)mpfrValue, sizeof(__mpfr_struct));
    offset += sizeof(__mpfr_struct);
#ifdef WIN32
    memcpy_s(offset, GetPackedSize() - sizeof(__mpfr_struct),
      ((mpfr_ptr)mpfrValue)->_mpfr_d, GetPackedSize() - sizeof(__mpfr_struct));
#else
    memcpy(offset,
      //GetPackedSize() - sizeof(__mpfr_struct),
      ((mpfr_ptr)mpfrValue)->_mpfr_d, GetPackedSize() - sizeof(__mpfr_struct));
#endif
    //TODO
  }
    break;
#endif
  case etDouble:
    memcpy(result, &dValue, sizeof(double));
    break;
  }
}

void Extended::Unpack(void* packed)
{
  //    unsigned char* intmem;

  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
  {
    char* offset = (char*)packed;
    // store pointer to data of current object
    void* tmp;
    if (mpfrValue == NULL)
      mpfrValue = new mpreal();
    // _mpfr_d after unpacking of a structure will be changed, so it must
    // be stored before unpacking
    tmp = (void*)((mpfr_ptr)mpfrValue)->_mpfr_d;
    memcpy(((mpfr_ptr)mpfrValue), offset, sizeof(__mpfr_struct));
    offset += sizeof(__mpfr_struct);
    // now we can assign _mpfr_d to the previous value
    ((mpfr_ptr)mpfrValue)->_mpfr_d = (mp_limb_t*)tmp;
    // and write data at that address
    memcpy(((mpfr_ptr)mpfrValue)->_mpfr_d, offset, GetPackedSize() - sizeof(__mpfr_struct));
#ifdef DEBUG_EXTENDED
    dValue = mpfrValue->toDouble();
#endif
  }
    break;
#endif
  case etDouble:
    memcpy(&dValue, packed, sizeof(double));
    break;
  }
}

ExtendedTypeID Extended::GetTypeID()
{
  return type;
}

void Extended::SetTypeID(ExtendedTypeID typeID)
{
  type = typeID;
  packedSize = -1;
}

int Extended::GetStringSize()
{
  int result=(digitsNumber)+10;
  //            |           |
  //            |           +- sign + . + e + e sign +
  //            |                   + exp + reserve
  //            +-- number of digits
  /*switch (type)
  {
  case etMPFR:
  result = ::log10(precision)+6;
  break;
  case etDouble:
  result = ::log10(precision)+6;
  break;
  }*/
  return result;
}

void Extended::toString(char* s) const
{
  switch (type)
  {

#ifndef NO_MPREAL
  case etMPFR:
    {
#ifdef WIN32
      std::string mpfrStr = mpfrValue->toString((size_t)digitsNumber);
      strcpy_s(s, GetStringSize(), mpfrStr.c_str());
#else
      std::string mpfrStr = mpfrValue->toString((size_t)digitsNumber);
      strcpy(s,mpfrStr.c_str());
#endif
    }
    break;
#endif
  case etDouble:
    //errno_t res = _gcvt_s(s, GetStringSize(), dValue, GetStringSize() - 2);
    //if (res != 0)
    //  s[0] = '\0';
#ifdef WIN32
    sprintf_s(s, GetStringSize() - 1, "%.15lf", dValue);
    _heapchk();
#else
    sprintf(s,"%.15lf", dValue);
#endif
    //TODO
    //heapchk();
    break;
  }

  // change '.' -> ',' don't know why, but do that more neatly
  int i = 0;
  while (s[i] != 0)
  {
    if (s[i] == '.')
    {
      s[i] = ','; // it doesn't depend on locale
      break;
    }
    i++;
  }
}

double Extended::toDouble() const
{
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    return mpfrValue->toDouble();
#endif
  case etDouble:
    return dValue;
    break;
  }
  return -1;
}

const Extended& Extended::operator=(const Extended& e)
{
  if (this == &e)
  {
    return e;
  }

  Free();
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    mpfrValue = new mpreal(*e.mpfrValue);
#ifdef DEBUG_EXTENDED
    dValue = mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    dValue = e.dValue;
  }

  return *this;
}

const Extended& Extended::operator=(const double& e)
{
  Free();
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    mpfrValue = new mpreal(e);
#ifdef DEBUG_EXTENDED
    dValue = mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    dValue = e;
  }

  return *this;
}

const Extended& Extended::operator=(const char* e)
{
  Free();
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    mpfrValue = new mpreal(e);
#ifdef DEBUG_EXTENDED
    dValue = mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    dValue = atof(e);
  }

  return *this;
}

Extended Extended::operator+(const Extended& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) += *(e.mpfrValue);
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue += e.dValue;
    break;
  }

  return tmp;
}

Extended Extended::operator+=(const Extended& e)
{
  *this = *this + e;

  return (*this);
}

Extended Extended::operator+(const double& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) += e;
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue += e;
    break;
  }

  return tmp;
}

void Extended::operator+=(const double& e)
{
  *this = this->operator +(e);
}

Extended Extended::operator-(const Extended& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) -= *(e.mpfrValue);
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue -= e.dValue;
    break;
  }

  return tmp;
}

Extended Extended::operator-=(const Extended& e)
{
  *this = *this - e;

  return (*this);
}

Extended Extended::operator-(const double& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) -= e;
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue -= e;
    break;
  }

  return tmp;
}

void Extended::operator-=(const double& e)
{
  *this = this->operator -(e);
}

Extended Extended::operator*(const Extended& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) *= *(e.mpfrValue);
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue *= e.dValue;
    break;
  }

  return tmp;
}

Extended Extended::operator*=(const Extended& e)
{
  *this = *this * e;

  return (*this);
}

Extended Extended::operator*(const double& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) *= e;
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue *= e;
    break;
  }

  return tmp;
}

void Extended::operator*=(const double& e)
{
  *this = this->operator *(e);
}

Extended Extended::operator/(const Extended& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) /= *(e.mpfrValue);
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue /= e.dValue;
    break;
  }

  return tmp;
}

Extended Extended::operator/(const double& e) const
{
  Extended tmp(*this);
  switch (type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *(tmp.mpfrValue) /= e;
#ifdef DEBUG_EXTENDED
    tmp.dValue = tmp.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    tmp.dValue /= e;
    break;
  }

  return tmp;
}

//------------------------------------------
// friend operations (+,-,*,/,floor,pow,root)
//------------------------------------------

Extended operator+(const double& e1, const Extended& e2)
{
  Extended tmp(e1);
  tmp += e2;

  return tmp;
}

Extended operator-(const double& e1, const Extended& e2)
{
  Extended tmp(e1);
  tmp -= e2;

  return tmp;
}

Extended operator*(const double& e1, const Extended& e2)
{
  Extended tmp(e1);
  tmp *= e2;

  return tmp;
}

Extended operator/(const double& e1, const Extended& e2)
{
  Extended tmp(e1);
  tmp = tmp / e2;

  return tmp;
}

Extended fabs(const Extended& e)
{
  Extended res;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *res.mpfrValue = mpfr::fabs(*e.mpfrValue);
    break;
#endif
  case etDouble:
    res.dValue = ::fabs(e.dValue);
    break;
  }

  return res;
}

double floor(const Extended& e)
{
  double res = 0.0;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    res = ::floor(e.mpfrValue->toDouble());
    break;
#endif
  case etDouble:
    res = floor(e.dValue);
    break;
  }

  return res;
}

Extended pow(const Extended& e, double degree)
{
  Extended res = 0.0;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    *res.mpfrValue = mpfr::pow(*e.mpfrValue, degree);
    // ::pow((double)*e.mpfrValue, degree);
#ifdef DEBUG_EXTENDED
    res.dValue = res.mpfrValue->toDouble();
#endif
    break;
#endif
  case etDouble:
    res.dValue = ::pow(e.dValue, degree);
    break;
  }

  return res;
}

double root(const Extended& e, int k)
{
  double res = 0.0;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    //*(tmp.mpfrValue) = mpfr::root(*(e.mpfrValue), k);
    //*(tmp.mpfrValue) = ::pow((double)*(e.mpfrValue), 1./k);
    //res = ::pow((double)*(e.mpfrValue), 1./k);
    res = ::pow(e.mpfrValue->toDouble(), 1./k);
    break;
#endif
  case etDouble:
    //tmp.dValue = ::pow(e.dValue, 1./k);
    res = ::pow(e.dValue, 1./k);
    break;
  }

  return res;
}

//---------------------------------------------
// friend comparison operators (==, !=, <, <=)
//---------------------------------------------

bool operator==(const Extended& e1, const Extended& e2)
{
  bool result = false;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    result = (*e1.mpfrValue == *e2.mpfrValue);
    break;
#endif
  case etDouble:
    result = (e1.dValue == e2.dValue);
    break;
  }

  return result;
}

bool operator==(const Extended& e1, double d)
{
  bool result = false;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    result = (e1.mpfrValue->toDouble() == d);
    break;
#endif
  case etDouble:
    result = (e1.dValue == d);
    break;
  }

  return result;
}

bool operator==(double d, const Extended& e1)
{
  bool result = false;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    result = (e1.mpfrValue->toDouble() == d);
    break;
#endif
  case etDouble:
    result = (e1.dValue == d);
    break;
  }

  return result;
}

bool operator!=(const Extended& e1, const Extended& e2)
{
  return !(e1 == e2);
}

bool operator!=(const Extended& e1, double d)
{
  return !(e1 == d);
}

bool operator!=(double d, const Extended& e1)
{
  return !(e1 == d);
}

bool operator<(const Extended& e1, const Extended& e2)
{
  bool result = false;
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    result = (*e1.mpfrValue < *e2.mpfrValue);
    break;
#endif
  case etDouble:
    result = (e1.dValue < e2.dValue);
    break;
  }

  return result;
}

bool operator<(const Extended& e, const double d)
{
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    return (e.mpfrValue->toDouble() < d);
#endif
  case etDouble:
    return (e.dValue < d);
  }

  return false;
}

bool operator<(const double d, const Extended& e)
{
  switch (Extended::type)
  {
#ifndef NO_MPREAL
  case etMPFR:
    return (d < e.mpfrValue->toDouble() );
#endif
  case etDouble:
    return (d < e.dValue);
    break;
  }

  return false;
}

bool operator<=(const Extended& e1, const Extended& e2)
{
  return (e1 < e2 || e1 == e2);
}

bool operator<=(const Extended& e, const double d)
{
  return (e < d || e == d);
}

bool operator<=(const double d, const Extended& e)
{
  return (d < e || d == e);
}

bool operator>(const Extended& e1, const Extended& e2)
{
  return !(e1 <= e2);
}

bool operator>(const Extended& e, const double d)
{
  return !(e <= d);
}

bool operator>(const double d, const Extended& e)
{
  return !(d <= e);
}

bool operator>=(const Extended& e1, const Extended& e2)
{
  return !(e1 < e2);
}

bool operator>=(const Extended& e, const double d)
{
  return !(e < d);
}

bool operator>=(const double d, const Extended& e)
{
  return !(d < e);
}
