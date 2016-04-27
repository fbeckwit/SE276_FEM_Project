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

const std::size_t Element::NEN;

/* *****************************  COPY CONTROL  ***************************** */


/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Returns the stiffness matrix for the given element using the current
 * consistent tangent. */
Eigen::MatrixXd Element::get_stiffness( std::size_t int_order ) const
{
  // Get elastic modulii tensor;
  Eigen::Matrix2d elastic_mod = material->get_tangent( );

  // Get Gauss points and weights, and the values of radius at the points;
  // TODO:  Make Gauss quadrature work on function objects and convert this to
  // utilize that (would remove need to grab Gauss points & weights and others);
  std::vector<double> gauss_pts = util::get_gauss_pts( int_order );
  std::vector<double> gauss_wts = util::get_gauss_wts( int_order );
  std::vector<double> func_eval( int_order );
  std::vector<double> radius( int_order );
  for( std::size_t pt{ 0 }; pt != int_order; ++pt )
    radius[pt] = interp_coord( gauss_pts[pt] );

  // Calculate the stiffness matrix;
  Eigen::MatrixXd stiffness( NEN, NEN );
  stiffness.setZero( );
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

/* Returns the external force acting on the element from tractions and body
 * forces. */
Eigen::MatrixXd Element::get_force_ext( ) const
{
  Eigen::VectorXd force( 2 );
  force.setZero( );

  // Check if node is on the natural boundary (Node::NBC) and calc the force;
  for( std::size_t a{ 0 }; a != NEN; ++a ) {
    if( nodes[a]->type == Node::NBC )
      force[a] = nodes[a]->bound_cond * nodes[a]->coord;
    else
      force[a] = 0;
  }
  return force;
}

/* -------------------------------------------------------------------------- */

/* Returns the internal force acting on the element due to strain energy. */
Eigen::MatrixXd Element::get_force_int( ) const
{
  return Eigen::MatrixXd( );
}

/* -------------------------------------------------------------------------- */

/* Given an output stream, print the node locations. */
void Element::print_nodes( std::ostream &out ) const
{
  for ( int a{ 0 }; a != NEN; ++a )
    out << nodes[ a ]->get_coord( ) << '\n';
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, interpolate the coordinate within the
 * element. */
double Element::interp_coord( double xi ) const
{
  // Sum N_a * x_a;
  double coord{0};
  for( int a{0}; a != NEN; ++a )
    coord += shape_func( xi, a ) * nodes[a]->coord;
  return coord;
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the shape function.  */
double Element::shape_func( double xi, std::size_t a )
{
  // Determine the appropriate value of xi_a;
  int xi_a = ( a == 0 ) ? -1 : 1;

  // Calculate shape function and return;
  return 0.5 * ( 1 + xi_a * xi );
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the gradient matrix, B. */
Eigen::Vector2d Element::get_gradient_matrix( double xi, std::size_t a ) const
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

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, return the stresses from the resulting
 * displacement.
 * PRECONDITION:  Nodes must have updated displacements. */
Eigen::Vector3d Element::get_stress( double xi ) const
{
  // Calculate the strain components and pass to the material to get the stress;
  Eigen::Vector2d strain;
  strain.setZero( );
  for( std::size_t a{ 0 }; a != NEN; ++a )
    strain += get_gradient_matrix( xi, a ) * nodes[a]->disp;

  return material->get_stress( strain );
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
