// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the GNU Lesser General Public
// License (GNU LGPL) version 3.
// ====================================================================

#pragma once
#include <complex>

namespace polylogarithm {

/// real polylogarithm with n=2 (dilogarithm)
double Li2(double) noexcept;

/// real polylogarithm with n=2 (dilogarithm) with long double precision
long double Li2(long double) noexcept;

/// complex polylogarithm with n=2 (dilogarithm)
std::complex<double> Li2(const std::complex<double>&) noexcept;

/// complex polylogarithm with n=2 (dilogarithm) with long double precision
std::complex<long double> Li2(const std::complex<long double>&) noexcept;

} // namespace polylogarithm

extern "C" {
  double dli2_(double& x);
  //double li2_(double& x);
  //double ddilog_(double& x);
}
