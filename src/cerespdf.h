#ifndef cerespdf_h
#define cerespdf_h


#include "parton.h"

#include <ceres/ceres.h>
#include <complex>

using ceres::SizedCostFunction;
using ceres::DynamicCostFunction;
using ceres::Solver;

using namespace std;

namespace cerespdf
{
  const int P = 3;
  const int M = 3;  //degree of logx polynomial
  const int N = 15; //degree of Bernstein polynomial

  void init(double xmin);
  void setoptions();
  void update(double muf);
  void updatepdfs(double muf);
  void free();
  void bernstein(int n, double t, double *bern);

  extern int np;
  extern int num_observations;
  
  //internal data
  extern double *xp;
  extern double *logx;
  extern double *log1mx;
  extern double *bber;

  class ResAnalytic : public SizedCostFunction<1,P,M,N>
    {
    public:
    ResAnalytic(const double x, const double y, const int i) : x_(x), y_(y), i_(i) {}
      virtual ~ResAnalytic() {}
      bool Evaluate(const double* const* parameters,
		    double* residuals,
		    double** jacobians) const;
    private:
      const double x_;
      const double y_;
      const int i_;
    };
  
  //extern int N, M, P;
  //class ResAnalytic : public DynamicCostFunction
  //{
  //public:
  //ResAnalytic(const double x, const double y, const int i) : x_(x), y_(y), i_(i) {}
  //    virtual ~ResAnalytic() {}
  //    virtual bool Evaluate(double const* const* parameters,
  //			    double* residuals,
  //			    double** jacobians) const;
  //  private:
  //    const double x_;
  //    const double y_;
  //    const int i_;
  //};

  //Same options for all minimisations
  extern Solver::Options options;
  
  class approxpdf
  {
  public:
    double p[P] = {1,1,1};
    double c[N] = {0};
    double l[M] = {0};

    double *f;

    void allocate() {f = new double[num_observations];};
    void free()     {delete [] f;};
    
    void copypar(approxpdf pdfin)
    {
      copy (pdfin.p, pdfin.p+P, p);
      copy (pdfin.c, pdfin.c+N, c);
      copy (pdfin.l, pdfin.l+M, l);
    };
    void resetpar()
    {
      fill (p, p+P, 0);
      fill (c, c+N, 0);
      fill (l, l+M, 0);
    };
    void init(double xmin, parton::pdgid p);
    
    void fit();
    complex <double> mellin(complex <double> s);
    complex <double> mellin_num(complex <double> s);
    double approx(double x);

    Solver::Summary summary;
  };

  extern approxpdf U,UB,D,DB,S,SB,C,CB,B,BB,G;
};
  


#endif
