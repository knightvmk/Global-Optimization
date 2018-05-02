/* Two-dimensional multiextremal test function of Vladimir A. Grishagin

   See:  Grishagin, V.A., Operating Characteristics of Some Global Search Algorithms.
         in: Problems in Random Search, Riga: Zinatne, 1978, issue 7, pp. 198--206;
   also: Strongin, R.G., and Sergeyev, Ya.D., Global Optimization with Non-Convex
         Constraints: Sequential and Parallel Algorithms. Dordrecht: Kluwer, 2000.
*/

#ifndef __GRISHAGIN_FUNCTION_H__
#define __GRISHAGIN_FUNCTION_H__

#include "..\problem_interface.hpp"

#define PROPERTY(T, N)     \
  T Get ## N() const;     \
  void Set ## N(T value);

namespace vagrish
{
  class GrishaginFunction : public IGOProblem<double>
  {
  private:
    int mFunctionNumber;
    unsigned char icnf[45];
    double af[7][7], bf[7][7], cf[7][7], df[7][7];

    double rndm20(unsigned char k[]);
    void gen(unsigned char k[], unsigned char k1[], int kap1, int kap2);

  public:
    GrishaginFunction();
    ~GrishaginFunction();

    PROPERTY(int, FunctionNumber);

    double CalculateXDerivative(const double* y) const;
    double CalculateYDerivative(const double* y) const;

    double Calculate(const double* y, int fNumber = 0) const override;
    int GetOptimumPoint(double* y) const override;
    double GetOptimumValue() const override;
    void GetBounds(double* lb, double* ub) const override;
    int GetDimension() const override;
    int GetConstraintsNumber() const override;
  };
}
#endif
