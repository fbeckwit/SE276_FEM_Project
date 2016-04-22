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
#include <Eigen/Dense>

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

};

#endif
