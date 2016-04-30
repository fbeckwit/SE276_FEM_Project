/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the `Linear' abstraction, which      *
 * inherits publicly from `Element.' Class definition given in Linear.h.      *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Linear.h"
#include "gauss_quadrature.h"

// System headers;
#include <vector>

/* *****************************  COPY CONTROL  ***************************** */

/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the shape function.  */
double Linear::shape_func( double xi, std::size_t a ) const
{
  // Determine the appropriate value of xi_a;
  int xi_a = ( a == 0 ) ? -1 : 1;

  // Calculate shape function and return;
  return 0.5 * ( 1 + xi_a * xi );
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the gradient matrix, B. */
Eigen::VectorXd Linear::get_gradient_matrix( double xi, std::size_t a ) const
{
  // Determine the appropriate value of xi_a;
  int xi_a = ( a == 0 ) ? -1 : 1;

  // Calculate coefficients used in gradient matrix;
  double radius = interp_coord( xi );
  double N_a = shape_func( xi, a );

  // Build the gradient matrix;
  Eigen::VectorXd B_a = Eigen::VectorXd::Zero( 2 );
  B_a[0] = xi_a / length;
  B_a[1] = N_a / radius;

  return B_a;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */

