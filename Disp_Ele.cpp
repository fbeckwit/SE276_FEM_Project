/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the Disp_Ele abstraction.            *
 * Class definition given in Disp_Ele.h.                                      *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Disp_Ele.h"
#include "gauss_quadrature.h"

// System headers;
#include <vector>

/* *****************************  COPY CONTROL  ***************************** */


/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Returns the stiffness matrix for the given element using the current
 * consistent tangent. */
Eigen::MatrixXd Disp_Ele::get_stiffness( std::size_t int_order )
{
  // Calculate the stiffness matrix;
  std::vector<Node *>::size_type num_nodes = nodes.size( );
  Eigen::MatrixXd stiffness = Eigen::MatrixXd::Zero( num_nodes, num_nodes );
  for( std::vector<Node *>::size_type a{ 0 }; a != nodes.size( ); ++a ) {
    for( std::vector<Node *>::size_type b{ a }; b != nodes.size( ); ++b ) {

      stiff_eval.a = a;
      stiff_eval.b = b;
      stiffness( a, b ) = util::integrate( stiff_eval, int_order );

      // If we're calculating off-diagonal terms, copy to the lower triangle;
      if ( a != b )
        stiffness( b, a ) = stiffness( a, b );
    }
  }
  return stiffness;
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, return the stresses from the resulting
 * displacement.
 * PRECONDITION:  Nodes must have updated displacements. */
Eigen::Vector3d Disp_Ele::interp_stress( double xi ) const
{
  // Calculate the strain components and pass to the material to get the stress;
  Eigen::Vector2d strain;
  strain.setZero( );
  for( std::vector<Node *>::size_type a{ 0 }; a != nodes.size( ); ++a )
    strain += get_gradient_matrix( xi, a ) * nodes[a]->disp;

  return material->get_stress( strain );
}

/* -------------------------------------------------------------------------- */

double Disp_Ele::K_Func::operator()( double xi ) const
{
  // Get required matrices and info;
  Eigen::Matrix2d elastic_mod = parent->material->get_tangent( );
  double radius = parent->interp_coord( xi );
  double rad_deriv = parent->interp_coord_deriv( xi );
  Eigen::Vector2d B_a = parent->get_gradient_matrix( xi, a );
  Eigen::Vector2d B_b = parent->get_gradient_matrix( xi, b );

  // Calculate value;
  double ret = ( B_a.transpose( ) * elastic_mod * B_b ).value( );
  return ret * radius * rad_deriv;
}
