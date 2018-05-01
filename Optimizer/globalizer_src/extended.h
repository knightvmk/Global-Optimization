/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      extended.h                                                  //
//                                                                         //
//  Purpose:   Header file for class, that represents numbers with         //
//             arbitrary accuracy                                          //
//                                                                         //
//  Author(s): Chigarev V.                                                 //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __EXTENDED_H__
#define __EXTENDED_H__

#define  _CRT_SECURE_NO_WARNINGS

#ifndef NO_MPREAL

namespace mpfr
{
  class mpreal;
}
using namespace mpfr;

#endif

/**
* Extended type
*/
enum ExtendedTypeID { etDouble = 0, etMPFR = 2 };

// ----------------------------------------------------------------------------
// class represents number that is used where arbitrary accuracy
// may be necessary
// MPFR is used by default

class Extended
{
protected:
  double          dValue;
#ifndef NO_MPREAL
  mpreal*         mpfrValue;
#endif

  static bool precisionSet;             // indicates whether mpfr_precision has been set or not
  static int mpfrPrecision;             // this precision is calculated from precision
                                        // it stores number of bits for mantissa
  static int mapmPrecision;             // this precision is calculated from precision
                                        // it stores number of bytes for presenting
                                        // mantissa, 2 digits per byte
  static double maxNumbers;             // amount of various numbers to be stored (e.g. 10^10)
                                        // by mantissa disregarding exponent
  static double precision;              // minimum distance between two nearby numbers
  static float precisionReserveFactor;  // coefficient of precision reserve
  static int digitsNumber;              // number of significant digits
  static int packedSize;                // number of bytes needed for representing packed number

public:
  static ExtendedTypeID type; // etMPFR etDouble

  void Init();
  void Free();

  Extended(double val = 0.0);
  Extended(const char* str);
  Extended(const Extended& e);
  ~Extended();

  static float getPrecisionReserveFactor();
  static void SetPrecision(double prec);
  static double GetPrecision();
  static int GetPackedSize();
  void Pack(void* result) const;
  void Unpack(void* packed);

  static ExtendedTypeID GetTypeID();
  static void SetTypeID(ExtendedTypeID typeID);

  static int GetStringSize();
  void toString(char* s) const;

  double toDouble() const;

  const Extended& operator=(const Extended& e);
  const Extended& operator=(const double& e);
  const Extended& operator=(const char* e);

  Extended operator+(const Extended& e) const;
  Extended operator+=(const Extended& e);
  Extended operator+(const double& e) const;
  void operator+=(const double& e);

  //int GetRefCount() const;

  Extended operator-(const Extended& e) const;
  Extended operator-=(const Extended& e);
  Extended operator-(const double& e) const;
  void operator-=(const double& e);

  Extended operator*(const Extended& e) const;
  Extended operator*=(const Extended& e);
  Extended operator*(const double& e) const;
  void operator*=(const double& e);

  Extended operator/(const Extended& e) const;
  Extended operator/(const double& e) const;

  //------------------------------------------
  // friend operations (+,-,*,/,floor,pow,root)
  //------------------------------------------
  friend Extended operator+(const double& e1, const Extended& e2);
  friend Extended operator-(const double& e1, const Extended& e2);
  friend Extended operator*(const double& e1, const Extended& e2);
  friend Extended operator/(const double& e1, const Extended& e2);

  friend Extended fabs(const Extended& e);
  friend double floor(const Extended& e);
  friend Extended pow(const Extended& e, double degree);
  // this method is used for distance between two points
  // there is a Sergey's proof of sufficiency of double precision for it
  friend double root(const Extended& e, int k);

  //---------------------------------------------
  // friend comparison operators (==, !=, <, <=)
  //---------------------------------------------

  friend bool operator==(const Extended& e1, const Extended& e2);
  friend bool operator==(const Extended& e1, double d);
  friend bool operator==(double d, const Extended& e1);

  friend bool operator!=(const Extended& e1, const Extended& e2);
  friend bool operator!=(const Extended& e1, double d);
  friend bool operator!=(double d, const Extended& e1);

  friend bool operator<(const Extended& e1, const Extended& e2);
  friend bool operator<(const Extended& e, const double d);
  friend bool operator<(const double d, const Extended& e);

  friend bool operator>(const Extended& e1, const Extended& e2);
  friend bool operator>(const Extended& e, const double d);
  friend bool operator>(const double d, const Extended& e);

  friend bool operator<=(const Extended& e1, const Extended& e2);
  friend bool operator<=(const Extended& e, const double d);
  friend bool operator<=(const double d, const Extended& e);

  friend bool operator>=(const Extended& e1, const Extended& e2);
  friend bool operator>=(const Extended& e, const double d);
  friend bool operator>=(const double d, const Extended& e);
};

#endif // __ABSEXTENDED_H__
