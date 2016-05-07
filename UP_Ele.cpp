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

double UP_Ele::K_Func::operator()( double xi ) const
{
  // Get required matrices and info;

  // Calculate value;
  return 0.0;
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
