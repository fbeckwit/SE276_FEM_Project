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
Eigen::MatrixXd fem::UP_Ele::get_stiffness( std::size_t int_order )
{
  // Calculate the various matrices;
  Eigen::MatrixXd stiff = quad::integrate_matrix( k_eval, int_order, true );
  Eigen::MatrixXd div_op = quad::integrate_matrix( g_eval, int_order );
  Eigen::MatrixXd constr_op = quad::integrate_matrix( m_eval, int_order, true );

  // Perform static condensation and return;
  stiff -= ( div_op * constr_op.inverse( ) * div_op.transpose( ) );
  return stiff;
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, return the stresses from the resulting
 * displacement.
 * PRECONDITION:  Nodes must have updated displacements. */
Eigen::Vector3d fem::UP_Ele::interp_stress( double xi ) const
{
  // Interpolate the strain and pressure at the coordinate, xi;
  Eigen::Vector2d strain = interp_strain( xi );
  double press = interp_pressure( xi );

  // Calculate the stress;
  Eigen::Vector3d stress = material->get_stress_mu( strain );
  stress += Eigen::Vector3d::Constant( press );
  return stress;
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, and the node number, a, return the
 * divergence matrix, b^v. */
double fem::UP_Ele::get_divergence_matrix( double xi, std::size_t a ) const
{
  // Get the gradient matrix and sum to get the divergence matrix at a;
  Eigen::Vector2d B_a = get_gradient_matrix( xi, a );
  return B_a.sum( );
}

/* -------------------------------------------------------------------------- */

/* Update the element info.
 * PRECONDITION:  Element nodes must be updated. */
void fem::UP_Ele::update( )
{
  // Update the pressures.  Get the G & M matrices;
  std::size_t int_order{ 2 }; // TODO:  hard-coded for now, remove later;
  Eigen::MatrixXd G = quad::integrate_matrix( g_eval, int_order );
  Eigen::MatrixXd M = quad::integrate_matrix( m_eval, int_order );

  // Perform matrix mult to get discrete pressure operator;
  Eigen::MatrixXd press_op = -( M.inverse( ) * G.transpose( ) );

  // Create a displacement vector for the element;
  Eigen::VectorXd disp = Eigen::VectorXd::Zero( nodes.size( ) );
  for( std::vector<Node *>::size_type a{ 0 }; a != nodes.size( ); ++a )
    disp(a) = nodes[a]->disp;

  // Calculate the pressures;
  for( std::vector<double>::size_type a{ 0 }; a != pressure.size( ); ++a )
    pressure[a] = ( press_op.row( a ) * disp ).value( );
}

/* -------------------------------------------------------------------------- */

/* Given the parametrix coordinate, xi, interpolate the pressure.
 * PRECONDITION:  Pressures must be updated after solving. */
double fem::UP_Ele::interp_pressure( double xi ) const
{
  // Perform summation over the pressure coefficients and their interpolation;
  double pres{ 0.0 };
  for( std::vector<double>::size_type a{ 0 }; a != pressure.size( ); ++a )
    pres += pressure_func( xi, a ) * pressure[a];
  return pres;
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */

/* ************************  NESTED CLASS FUNCTIONS  ************************ */

double fem::UP_Ele::K_Func::operator()( double xi ) const
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

double fem::UP_Ele::G_Func::operator()( double xi ) const
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

double fem::UP_Ele::M_Func::operator()( double xi ) const
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
