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

// System headers;

/* *****************************  COPY CONTROL  ***************************** */


/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* get_stiffness( )
 * Returns the stiffness matrix for the given element using the current
 * consistent tangent.
 */
Eigen::MatrixXd Element::get_stiffness( )
{
  return Eigen::MatrixXd( );
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
double Element::shape_func( double xi, unsigned int a )
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
Eigen::Vector2d Element::get_gradient_matrix( double xi, unsigned int a )
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
