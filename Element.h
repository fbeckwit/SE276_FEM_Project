/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** */

#ifndef GUARD_ELEMENT_H
#define GUARD_ELEMENT_H

// Project-specific headers;
#include "Node.h"

// System headers;
#include <array>
#include <cstddef>
#include <iostream>
#include <Eigen/LU>

class Element {

public:

  static const int NEN = 2;

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */
  Element( ) : nodes{ nullptr, nullptr }, length{0.0}
  { }

  Element( Node *n0, Node *n1 ) : nodes{ n0, n1 }
  {
    length = n1->get_coord( ) - n0->get_coord( );
  }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* get_stiffness( )
   * Returns the stiffness matrix for the given element using the current
   * consistent tangent.
   */
  Eigen::MatrixXd get_stiffness( std::size_t int_order = 2 );

  /* get_force_ext( )
   * Returns the external force acting on the element from tractions and body
   * forces.
   */
  Eigen::MatrixXd get_force_ext( );

  /* get_force_int( )
   * Returns the internal force acting on the element due to strain energy.
   */
  Eigen::MatrixXd get_force_int( );

  /* print_nodes( )
   * Given an output stream, print the node locations.
   */
  void print_nodes( std::ostream &out = std::cout );

  /* calc_coord( )
   * Given the parametric coordinate, xi, interpolate the coordinate within the
   * element.
   */
  double interp_coord( double xi );

  /* shape_func( )
   * Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function.
   */
  static double shape_func( double xi, std::size_t a );

  /* shape_deriv( )
   * Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative.
   */
  static inline double shape_deriv( std::size_t a ) {
    return ( a == 0 ) ? -0.5 : 0.5;
  }

  /* get_gradient_matrix( )
   * Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the gradient matrix, B.
   */
  Eigen::Vector2d get_gradient_matrix( double xi, std::size_t a );

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  std::array<Node *, NEN> nodes;
  double length;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

#endif
