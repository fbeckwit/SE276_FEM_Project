/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the Element abstraction.             *
 * Class definition given in Element.h.                                       *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Element.h"
#include "gauss_quadrature.h"

// System headers;
#include <vector>

/* *****************************  COPY CONTROL  ***************************** */


/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* get_stiffness( )
 * Returns the stiffness matrix for the given element using the current
 * consistent tangent.
 */
Eigen::MatrixXd Element::get_stiffness( std::size_t int_order )
{
  // Get elastic modulii tensor;
  // TODO:  Obtain this from a call from material model, hard-coded for now.
  double E{ 1000 };
  double nu{ 0.4999 };
  double lambda = nu * E / (( 1 + nu ) * ( 1 - 2 * nu ));
  double G = E / 2.0 / ( 1 + nu );

  Eigen::Matrix<double, 2, 2> elastic_mod;
  elastic_mod( 0, 0 ) = elastic_mod( 1, 1 ) = lambda + 2 * G;
  elastic_mod( 0, 1 ) = elastic_mod( 1, 0 ) = lambda;

  // Get Gauss points and weights, and the values of radius at the points;
  std::vector<double> gauss_pts = util::gauss_pts( int_order );
  std::vector<double> gauss_wts = util::gauss_wts( int_order );
  std::vector<double> func_eval( int_order );
  std::vector<double> radius( int_order );
  for( std::size_t pt{ 0 }; pt != int_order; ++pt )
    radius[pt] = interp_coord( gauss_pts[pt] );

  // Calculate the stiffness matrix;
  Eigen::MatrixXd stiffness( NEN, NEN );
  for( std::size_t a{ 0 }; a != NEN; ++a ) {
    for( std::size_t b{ a }; b != NEN; ++b) {

      // Calculate k_{ab}.  Get function evaluations;
      for( std::size_t pt{ 0 }; pt != int_order; ++pt) {
        Eigen::Vector2d B_a = get_gradient_matrix( gauss_pts[pt], a );
        Eigen::Vector2d B_b = get_gradient_matrix( gauss_pts[pt], b );

        func_eval[pt] = ( B_a.transpose( ) * elastic_mod * B_b ).value( );
        func_eval[pt] *= radius[pt] * length / 2.0;
      }
      stiffness( a, b ) = util::integrate( func_eval, gauss_wts );

      // If we're calculating off-diagonal terms, copy to the lower triangle;
      if ( a != b )
        stiffness( b, a ) = stiffness( a, b );
    }
  }

  return stiffness;
}

/* -------------------------------------------------------------------------- */

/* get_force_ext( )
 * Returns the external force acting on the element from tractions and body
 * forces.
 */
Eigen::MatrixXd Element::get_force_ext( )
{
  return Eigen::MatrixXd( );
}

/* -------------------------------------------------------------------------- */

/* get_force_int( )
 * Returns the internal force acting on the element due to strain energy.
 */
Eigen::MatrixXd Element::get_force_int( )
{
  return Eigen::MatrixXd( );
}

/* -------------------------------------------------------------------------- */

/* print_nodes( )
 * Given an output stream, print the node locations.
 */
void Element::print_nodes( std::ostream &out )
{
  for ( int a{ 0 }; a != NEN; ++a )
    out << nodes[ a ]->get_coord( ) << '\n';
}

/* -------------------------------------------------------------------------- */

/* calc_coord( )
 * Given the parametric coordinate, xi, interpolate the coordinate within the
 * element.
 */
double Element::interp_coord( double xi )
{
  // Sum N_a * x_a;
  double coord{0};
  for( int a{0}; a != NEN; ++a )
    coord += shape_func( xi, a ) * nodes[a]->get_coord( );
  return coord;
}

/* -------------------------------------------------------------------------- */

/* shape_func( )
 * Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the shape function.
 */
double Element::shape_func( double xi, std::size_t a )
{
  // Determine the appropriate value of xi_a;
  int xi_a = ( a == 0 ) ? -1 : 1;

  // Calculate shape function and return;
  return 0.5 * ( 1 + xi_a * xi );
}

/* -------------------------------------------------------------------------- */

/* get_gradient_matrix( )
 * Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the gradient matrix, B.
 */
Eigen::Vector2d Element::get_gradient_matrix( double xi, std::size_t a )
{
  // Determine the appropriate value of xi_a;
  int xi_a = ( a == 0 ) ? -1 : 1;

  // Calculate coefficients used in gradient matrix;
  double radius = interp_coord( xi );
  double N_a = shape_func( xi, a );

  // Build the gradient matrix;
  Eigen::Vector2d B_a;
  B_a[0] = xi_a / length;
  B_a[1] = N_a / radius;

  return B_a;
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
