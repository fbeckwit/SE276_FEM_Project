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

// System headers;
#include <vector>

/* *****************************  COPY CONTROL  ***************************** */


/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

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

/* ************************  NESTED CLASS FUNCTIONS  ************************ */

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
