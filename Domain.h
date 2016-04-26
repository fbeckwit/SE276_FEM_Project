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

#ifndef GUARD_DOMAIN_H
#define GUARD_DOMAIN_H

// Project-specific headers;
#include "Element.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <Eigen/LU>
#include <iostream>
#include <vector>

class Domain {

public:

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Domain( ) : nodes{ }, elements{ }, materials{ }, stiff{ }, force{ }, disp{ }
  { }

  /* Domain should be unique, disallow copy and assignment operators */
  Domain( const Domain & other ) = delete;
  Domain( Domain && other ) = delete;
  Domain & operator=( const Domain & rhs ) = delete;
  Domain && operator=( Domain && rhs ) = delete;

  /* Destructor */
  ~Domain( );

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* Given material properties, Young's modulus and Poisson's ratio, create an
   * elastic material and store in `mats.' */
  void create_material( double E, double nu );

  /* Given a coordinate, create a node and store in `nodes.' */
  void create_node( double coord );

  /* Given a coordinate, node type, and BC, create a node and store in
   * `nodes.' */
  void create_node( double coord, Node::node_type type, double bc );

  /* Given the node ids and a material id, create an element and store in
   * `elements.'
   * PRECONDITION:  Nodes `n0' and `n1' and material `mat_id' must be created */
  void create_element( std::size_t n0, std::size_t n1, std::size_t mat_id );

  /* Builds the stiffness matrix by looping elements and assembling.
   * PRECONDITION:  `elements' must be properly initialized. */
  void build_stiffness( std::size_t int_order = 2 );

  /* Builds the force vector by looping elements and assembling.
   * PRECONDITION:  `elements' must be properly initialized. */
  void build_force( );

  /* Builds the system of equations and then solves.
   * PRECONDITION:  `elements' must be properly initialized. */
  Eigen::VectorXd solve( std::size_t int_order = 2 );

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  /* Domain elements */
  std::vector<Node *> nodes;
  std::vector<Element *> elements;
  std::vector<Material *> materials;

  /* Matrix elements and solvers */
  Eigen::MatrixXd stiff;
  Eigen::VectorXd force;
  Eigen::VectorXd disp;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

#endif
