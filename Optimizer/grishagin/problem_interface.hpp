#pragma once

template <class fptype>
class IGOProblem
{
public:
  ~IGOProblem() {}

  virtual fptype Calculate(const fptype* y, int fNumber) const = 0;
  virtual int GetConstraintsNumber() const = 0;
  virtual int GetDimension() const = 0;
  virtual void GetBounds(fptype* left, fptype* right) const = 0;
  virtual int GetOptimumPoint(fptype* y) const = 0;
  virtual fptype GetOptimumValue() const = 0 ;
};
