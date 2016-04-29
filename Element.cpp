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
Eigen::MatrixXd Element::get_stiffness( std::size_t int_order )
{
  // Calculate the stiffness matrix;
  Eigen::MatrixXd stiffness = Eigen::MatrixXd::Zero( NEN, NEN );
  for( std::size_t a{ 0 }; a != NEN; ++a ) {
    for( std::size_t b{ a }; b != NEN; ++b) {

      stiff_eval.set_indeces( a, b );
      stiffness( a, b ) = util::integrate( stiff_eval, int_order );

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

/* Given the number of points to print, interpolate the coordinates within the
 * element. */
std::vector<double> Element::interp_coord( std::size_t num_pts) const
{
  // Get the points over the parametric domain;
  std::vector<double> xi = get_points( num_pts );

  // Loop the points, interpolate the displacement, and return;
  std::vector<double> coord( num_pts );
  for( std::size_t i{ 0 }; i != num_pts; ++i )
    coord[i] = interp_coord( xi[i] );
  return coord;
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, interpolate the displacement within the
 * element. */
double Element::interp_disp( double xi ) const
{
  // Sum N_a * d_a;
  double disp{0};
  for( int a{0}; a != NEN; ++a )
    disp += shape_func( xi, a ) * nodes[a]->disp;
  return disp;
}

/* -------------------------------------------------------------------------- */

/* Given the number of points to print, interpolate the displacements within the
 * element.
 * PRECONDITION:  Nodes must have updated displacements. */
std::vector<double> Element::interp_disp( std::size_t num_pts) const
{
  // Get the points over the parametric domain;
  std::vector<double> xi = get_points( num_pts );

  // Loop the points, interpolate the displacement, and return;
  std::vector<double> disp( num_pts );
  for( std::size_t i{ 0 }; i != num_pts; ++i )
    disp[i] = interp_disp( xi[i] );
  return disp;
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, return the stresses from the resulting
 * displacement.
 * PRECONDITION:  Nodes must have updated displacements. */
Eigen::Vector3d Element::interp_stress( double xi ) const
{
  // Calculate the strain components and pass to the material to get the stress;
  Eigen::Vector2d strain;
  strain.setZero( );
  for( std::size_t a{ 0 }; a != NEN; ++a )
    strain += get_gradient_matrix( xi, a ) * nodes[a]->disp;

  return material->get_stress( strain );
}

/* -------------------------------------------------------------------------- */

/* Given the number of points to print, interpolate the stresses from the
 * resulting displacement.
 * PRECONDITION:  Nodes must have updated displacements. */
std::vector<Eigen::Vector3d> Element::interp_stress( std::size_t num_pts) const
{
  // Get the points over the parametric domain;
  std::vector<double> xi = get_points( num_pts );

  // Loop the points, interpolate the displacement, and return;
  std::vector<Eigen::Vector3d> stress( num_pts );
  for( std::size_t i{ 0 }; i != num_pts; ++i )
    stress[i] = interp_stress( xi[i] );
  return stress;
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

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */

/* Given the number of intervals, return a set of equally spaced points over the
 * parametric domain. */
std::vector<double> Element::get_points( std::size_t num_pts )
{
  // Get interval size;
  double cell_size = 2.0 / ( num_pts - 1 );

  // Create vector of points in parametric domain and return;
  std::vector<double> xi( num_pts );
  for( std::vector<double>::size_type i{ 0 }; i != num_pts; ++i )
    xi[i] = -1.0 + i * cell_size;
  return xi;
}

double Element::Stiff_Eval::operator()( double xi ) const
{
  // Get required matrices and info;
  Eigen::Matrix2d elastic_mod = parent->material->get_tangent( );
  double radius = parent->interp_coord( xi );
  Eigen::Vector2d B_a = parent->get_gradient_matrix( xi, a );
  Eigen::Vector2d B_b = parent->get_gradient_matrix( xi, b );

  // Calculate value;
  double ret = ( B_a.transpose( ) * elastic_mod * B_b ).value( );
  return ret * radius * parent->length / 2.0;
}
