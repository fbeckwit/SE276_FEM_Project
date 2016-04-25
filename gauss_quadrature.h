//
//  gauss_quadrature.h
//
//
//  Created by Jacob Koester on 1/5/15.
//
//

#ifndef GUARD_GAUSS_QUADRATURE_H
#define GUARD_GAUSS_QUADRATURE_H

#include <array>
#include <vector>

namespace util {
  std::vector<double> get_gauss_pts( int order );
  std::vector<double> get_gauss_pts(
      int order,
      const std::array<double, 2> & interval_ends
      );

  std::vector<double> get_gauss_wts( int order );
  std::vector<double> get_gauss_wts(
      int order,
      const std::array<double, 2> & interval_ends
      );

  double integrate( const std::vector<double> & values,
      const std::vector<double> & weights );

  template <typename Func>
  double integrate( const Func & f, int int_order )
  {
    // Get integration points and weights;
    std::vector<double> points  = get_gauss_pts( int_order );
    std::vector<double> weights = get_gauss_wts( int_order );

    // Perform integration;
    double ret{ 0.0 };
    for( int pt{ 0 }; pt != int_order; ++pt )
      ret += f( points[pt] ) * weights[pt];
    return ret;
  }

  void test( int num_tests );
  double test_order( int poly_order, int num_intervals );
}

#endif
