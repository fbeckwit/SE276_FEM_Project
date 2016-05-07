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
#include <Eigen/Dense>

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
  double integrate( const Func & f, int int_order );

  void test( int num_tests );
  double test_order( int poly_order, int num_intervals );

  /* Given the dimensions of the matrix, a & b, and the function object used to
   * integrate each term, func, integrate the coefficients and return the
   * matrix. */
  template <typename Func>
  Eigen::MatrixXd integrate_matrix( Func & func, int int_order );

}

/* *************************  TEMPLATED FUNCTIONS  ************************** */

template <typename Func>
double util::integrate( const Func & f, int int_order )
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

/* -------------------------------------------------------------------------- */

/* Given the dimensions of the matrix, a & b, and the function object used to
 * integrate each term, func, integrate the coefficients and return the
 * matrix. */
template <typename Matrix_Func>
Eigen::MatrixXd util::integrate_matrix( Matrix_Func &func, int int_order )
{
  // Get the size of the matrix;
  const std::size_t num_rows = func.get_rows( );
  const std::size_t num_cols = func.get_cols( );

  // Calculate the matrix matrix;
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero( num_rows, num_cols );
  for( std::size_t a{ 0 }; a != num_rows; ++a ) {
    for( std::size_t b{ a }; b != num_cols; ++b ) {

      func.a = a;
      func.b = b;
      matrix( a, b ) = util::integrate( func, int_order );

      // If we're calculating off-diagonal terms, copy to the lower triangle;
      if ( a != b )
        matrix( b, a ) = matrix( a, b );
    }
  }
  return matrix;
}

#endif
