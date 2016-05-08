/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the UP_Ele abstraction.              *
 * Class definition given in UP_Ele.h.                                        *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "UP_Ele.h"
#include "gauss_quadrature.h"

// System headers;
#include <vector>

/* *****************************  COPY CONTROL  ***************************** */


/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Returns the stiffness matrix for the given element using the current
 * consistent tangent. */
Eigen::MatrixXd UP_Ele::get_stiffness( std::size_t int_order )
{
  return Eigen::MatrixXd::Zero( nodes.size( ), nodes.size( ) );
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, return the stresses from the resulting
 * displacement.
 * PRECONDITION:  Nodes must have updated displacements. */
Eigen::Vector3d UP_Ele::interp_stress( double xi ) const
{
  return Eigen::VectorXd::Zero( 3 );
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, and the node number, a, return the
 * divergence matrix, b^v. */
double UP_Ele::get_divergence_matrix( double xi, std::size_t a ) const
{
  // Get the gradient matrix and sum to get the divergence matrix at a;
  Eigen::Vector2d B_a = get_gradient_matrix( xi, a );
  return B_a.sum( );
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */

/* Returns the stiffness matrix for the given element from the mu-part of the
 * current consistent tangent of the material. */
Eigen::MatrixXd UP_Ele::get_K_mu( std::size_t int_order ) const
{
  return Eigen::MatrixXd::Zero( 1, 1 );
}

/* -------------------------------------------------------------------------- */

/* Returns the discrete gradient operator, G, which acts on pressures. */
Eigen::MatrixXd UP_Ele::get_G( std::size_t int_order ) const
{
  return Eigen::MatrixXd::Zero( 1, 1 );
}

/* -------------------------------------------------------------------------- */

/* Returns the discrete constraint operator, M. */
Eigen::MatrixXd UP_Ele::get_M( std::size_t int_order ) const
{
  return Eigen::MatrixXd::Zero( 1, 1 );
}

/* ************************  NESTED CLASS FUNCTIONS  ************************ */

double UP_Ele::K_Func::operator()( double xi ) const
{
  // Get required matrices and info;
  double mu = parent->material->get_mu( );
  double radius = parent->interp_coord( xi );
  double rad_deriv = parent->interp_coord_deriv( xi );
  Eigen::Vector2d B_a = parent->get_gradient_matrix( xi, a );
  Eigen::Vector2d B_b = parent->get_gradient_matrix( xi, b );

  // Calculate value and return;
  return ( B_a.transpose( ) * B_b ).value( ) * (2*mu) * radius * rad_deriv;
}

/* -------------------------------------------------------------------------- */

double UP_Ele::G_Func::operator()( double xi ) const
{
  // Get required matrices and info;
  double radius = parent->interp_coord( xi );
  double rad_deriv = parent->interp_coord_deriv( xi );
  double bv_a = parent->get_divergence_matrix( xi, a );
  double psi_b = parent->pressure_func( xi, b );

  // Calculate the value and return;
  return bv_a * psi_b * radius * rad_deriv;
}

/* -------------------------------------------------------------------------- */

double UP_Ele::M_Func::operator()( double xi ) const
{
  // Get required matrices and info;
  double radius = parent->interp_coord( xi );
  double rad_deriv = parent->interp_coord_deriv( xi );
  double psi_a = parent->pressure_func( xi, a );
  double psi_b = parent->pressure_func( xi, b );
  double bulk =
    parent->material->get_lambda( ) + 2.0/3.0 * parent->material->get_mu( );

  // Calculate the value and return;
  return -psi_a * psi_b * radius * rad_deriv / bulk;
}
