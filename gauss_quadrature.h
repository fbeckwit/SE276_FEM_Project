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
  std::vector<double> gauss_pts( int order );
  std::vector<double> gauss_pts(
      int order,
      const std::array<double, 2> & interval_ends
      );

  std::vector<double> gauss_wts( int order );
  std::vector<double> gauss_wts(
      int order,
      const std::array<double, 2> & interval_ends
      );

  //MatrixXd pts_and_wts( int order, VectorXd interval_ends );
  //MatrixXd pts_and_wts( int order );

  double integrate( const std::vector<double> & values,
      const std::vector<double> & weights );

  void test( int num_tests );
  double test_order( int poly_order, int num_intervals );
}

#endif
