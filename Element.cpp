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

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
