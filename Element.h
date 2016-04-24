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
#include <iostream>
#include <vector>
#include <Eigen/LU>

class Element {

public:

  static const int NEN = 2;

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */
  Element( ) : nodes{ nullptr, nullptr }
  { }

  Element( Node *n1, Node *n2 ) : nodes{ n1, n2 }
  { }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* get_stiffness( )
   * Returns the stiffness matrix for the given element using the current
   * consistent tangent.
   */
  Eigen::MatrixXd get_stiffness( );

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

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  std::array<Node *, NEN> nodes;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

  /* calc_coord( )
   * Given the parametric coordinate, xi, interpolate the coordinate within the
   * element.
   */
  double interp_coord( double xi );

  /* shape_func( )
   * Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function.
   */
  static double shape_func( double xi, unsigned int a );

  /* shape_deriv( )
   * Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative.
   */
  static inline double shape_deriv( unsigned int a ) {
    return ( a == 0 ) ? -0.5 : 0.5;
  }

  /* get_gradient_matrix( )
   * Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the gradient matrix, B.
   */
  Eigen::MatrixXd get_gradient_matrix( double xi, unsigned int a );

};

#endif
