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
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <iostream>
#include <Eigen/LU>
#include <utility>

class Element {

public:

  static const std::size_t NEN = 2;

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */
  Element( ) :
    nodes{ nullptr, nullptr }, material{ nullptr }, length{0.0}, ele_ID{ 0 }
  { }

  Element( std::size_t id, Node *n0, Node *n1, const Material *mat ) :
    nodes{ n0, n1 }, material{ mat->clone( ) }, length{ 0.0 }, ele_ID{ id }
  {
    length = n1->get_coord( ) - n0->get_coord( );
  }

  /* Copy Constructor */
  Element( const Element & other ) :
    nodes{ other.nodes }, material{ other.material->clone( ) },
    length{ other.length }, ele_ID{ other.ele_ID }
  { }

  /* Move Constructor */
  Element( Element && other ) :
    nodes{ std::move( other.nodes ) }, material{ other.material },
    length{ other.length }, ele_ID{ other.ele_ID }
  {
    other.nodes = { nullptr, nullptr };
    other.material = nullptr;
    other.length = 0.0;
    other.ele_ID = 0;
  }

  /* Assignment operators (Deleted) */
  Element & operator=( const Element & other ) = delete;
  Element && operator=( Element && other ) = delete;

  ~Element( ) { delete material; }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* Returns the stiffness matrix for the given element using the current
   * consistent tangent. */
  Eigen::MatrixXd get_stiffness( std::size_t int_order = 2 ) const;

  /* Returns the external force acting on the element from tractions and body
   * forces. */
  Eigen::MatrixXd get_force_ext( ) const;

  /* Returns the internal force acting on the element due to strain energy. */
  Eigen::MatrixXd get_force_int( ) const;

  /* Given the local node number, return the node type. */
  inline Node::node_type get_node_type( std::size_t a ) const {
    return nodes[a]->get_type( );
  }

  /* Given an output stream, print the node locations. */
  void print_nodes( std::ostream &out = std::cout ) const;

  /* Given the parametric coordinate, xi, interpolate the coordinate within the
   * element. */
  double interp_coord( double xi ) const;

  /* Given the number of points to print, interpolate the coordinates within the
   * element. */
  std::vector<double> interp_coord( std::size_t num_pts = 11 ) const;

  /* Given the parametric coordinate, xi, interpolate the displacement within
   * the element.
   * PRECONDITION:  Nodes must have updated displacements. */
  double interp_disp( double xi ) const;

  /* Given the number of points to print, interpolate the displacements within
   * the element.
   * PRECONDITION:  Nodes must have updated displacements. */
  std::vector<double> interp_disp( std::size_t num_pts = 11 ) const;

  /* Given the parametric coordinate, xi, interpolate the stresses from the
   * resulting displacement.
   * PRECONDITION:  Nodes must have updated displacements. */
  Eigen::Vector3d interp_stress( double xi ) const;

  /* Given the number of points to print, interpolate the stresses from the
   * resulting displacement.
   * PRECONDITION:  Nodes must have updated displacements. */
  std::vector<Eigen::Vector3d> interp_stress( std::size_t num_pts = 11 ) const;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function. */
  static double shape_func( double xi, std::size_t a );

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative. */
  static inline double shape_deriv( std::size_t a ) {
    return ( a == 0 ) ? -0.5 : 0.5;
  }

  /* Given the local node number (and eventually DOF number), return the global
   * equation number using the `LM' array. */
  inline std::size_t location_matrix( std::size_t a ) const {
    return nodes[a]->node_ID;
  }

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the gradient matrix, B. */
  Eigen::Vector2d get_gradient_matrix( double xi, std::size_t a ) const;

  inline std::size_t get_id( ) const { return ele_ID; }

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  std::vector<Node *> nodes;
  Material *material;
  double length;
  std::size_t ele_ID;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

  /* Given the number of intervals, return a set of equally spaced points over
   * the parametric domain. */
  static std::vector<double> get_points( std::size_t num_pts = 11 );

};

#endif
