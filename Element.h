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
#include <iomanip>
#include <iostream>
#include <Eigen/LU>
#include <utility>

class Element {

public:

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */
  Element( ) :
    nodes{ }, stiff_eval( this ), material{ nullptr },
    ele_ID{ 0 }, length{0.0}
  { }

  Element( std::size_t id, std::vector<Node *> nodes, const Material *mat ) :
    nodes{ nodes }, stiff_eval( this ),  material{ mat->clone( ) },
    ele_ID{ id }, length{ 0.0 }
  {
    length = nodes.back( )->get_coord( ) - nodes.front( )->get_coord( );
  }

  /* Copy Constructor */
  Element( const Element & other ) :
    nodes{ other.nodes }, stiff_eval( this ),
    material{ other.material->clone( ) },
    ele_ID{ other.ele_ID }, length{ other.length }
  { }

  /* Move Constructor */
  Element( Element && other ) :
    nodes{ std::move( other.nodes ) }, stiff_eval( this ),
    material{ other.material },
    ele_ID{ other.ele_ID }, length{ other.length }
  {
    other.nodes = { nullptr, nullptr };
    other.material = nullptr;
    other.length = 0.0;
    other.ele_ID = 0;
  }

  /* Assignment operators (Deleted) */
  Element & operator=( const Element & other ) = delete;
  Element && operator=( Element && other ) = delete;

  virtual ~Element( ) { delete material; }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* Returns the stiffness matrix for the given element using the current
   * consistent tangent. */
  virtual Eigen::MatrixXd get_stiffness( std::size_t int_order ) = 0;

  /* Returns the external force acting on the element from tractions and body
   * forces. */
  Eigen::MatrixXd get_force_ext( ) const;

  /* Returns the internal force acting on the element due to strain energy. */
  Eigen::MatrixXd get_force_int( ) const;

  /* Return the number of element nodes. */
  std::size_t get_num_nodes( ) const {
    return nodes.size( );
  }

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

  /* Given the parametric coordinate, xi, interpolate the derivative of the
   * coordinate within the element. */
  double interp_coord_deriv( double xi ) const;

  /* Given the number of points to print, interpolate the derivatives of the
   * coordinates within the element. */
  std::vector<double> interp_coord_deriv( std::size_t num_pts) const;

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
  virtual double shape_func( double xi, std::size_t a ) const = 0;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative. */
  virtual double shape_deriv( double xi, std::size_t a ) const = 0;

  /* Given the local node number (and eventually DOF number), return the global
   * equation number using the `LM' array. */
  inline std::size_t location_matrix( std::size_t a ) const {
    return nodes[a]->node_ID;
  }

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the gradient matrix, B. */
  Eigen::VectorXd get_gradient_matrix( double xi, std::size_t a) const;

  inline std::size_t get_id( ) const { return ele_ID; }

  /* Given a function object representing the exact solution, the field width,
   * an output stream, and the number of points to print, calculate the
   * displacements (exact and FEM) and print to the output stream. */
  template <typename Func>
    void print_disp( const Func & exact, std::size_t width,
        std::ostream & out = std::cout, std::size_t num_pts = 11 ) const
    {
      // Get the radius and displacements over the element;
      std::vector<double> radius = interp_coord( num_pts );
      std::vector<double> disp = interp_disp( num_pts );

      // Get the exact solution and the error from the function object;
      std::vector<double> disp_exact( num_pts );
      std::vector<double> disp_error( num_pts );
      for( std::size_t i{ 0 }; i != num_pts; ++i ) {
        disp_exact[i] = exact( radius[i] );
        disp_error[i] = disp_exact[i] - disp[i];
      }

      // Output the information;
      for( std::size_t i{ 0 }; i != num_pts; ++i ) {
        out << std::setw( width ) << radius[i];
        out << std::setw( width ) << disp[i];
        out << std::setw( width ) << disp_exact[i];
        out << std::setw( width ) << disp_error[i];
        out << '\n';
      }
    }

  /* Given a function object representing the exact solution, the field width,
   * an output stream, and the number of points to print, calculate the stresses
   * (exact and FEM) and print to the output stream. */
  template <typename Func>
    void print_stress( const Func & exact, std::size_t width,
        std::ostream & out = std::cout, std::size_t num_pts = 11 ) const
    {
      // Typedef;
      using Vector_array = std::vector<Eigen::Vector3d>;

      // Get the radius and displacements over the element;
      std::vector<double> radius = interp_coord( num_pts );
      Vector_array stress = interp_stress( num_pts );

      // Get the exact solution and the error from the function object;
      Vector_array stress_exact( num_pts, Eigen::Vector3d::Zero( ) );
      Vector_array stress_error( num_pts, Eigen::Vector3d::Zero( ) );
      for( std::size_t i{ 0 }; i != num_pts; ++i ) {
        stress_exact[i] = exact( radius[i] );
        stress_error[i] = stress_exact[i] - stress[i];
      }

      // Output the information;
      for( std::size_t i{ 0 }; i != num_pts; ++i ) {
        out << std::setw( width ) << radius[i];
        for( std::size_t j{ 0 }; j != 3; ++j ) {
          out << std::setw( width ) << stress[i][j];
          out << std::setw( width ) << stress_exact[i][j];
          out << std::setw( width ) << stress_error[i][j];
        }
        out << '\n';
      }
    }

private:

  /* ***************************  NESTED CLASSES  *************************** */

  /* Function object used in the evaluation of the stiffness matrix.  operator()
   * overloaded to return the internal energy density at the evaluation point,
   * xi. */
  struct Stiff_Eval {

    /* Constructor */
    Stiff_Eval( const Element * p, std::size_t _a = 0, std::size_t _b = 0 ) :
      parent{ p }, a{ _a }, b{ _b }
    { }

    /* Set the indeces as given.  Called before using integrate so that the
     * correct submatrix k_{ab} is calculated. */
    inline void set_indeces( std::size_t _a, std::size_t _b ) {
      a = _a;
      b = _b;
    }

    /* Given a parametric coordinate, xi, calculate the internal energy density
     * of the stiffness. */
    double operator()( double xi ) const;

    const Element * parent;
    std::size_t a;
    std::size_t b;
  };

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  std::vector<Node *> nodes;
  Stiff_Eval stiff_eval;
  Material *material;
  std::size_t ele_ID;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

  /* Given the number of intervals, return a set of equally spaced points over
   * the parametric domain. */
  static std::vector<double> get_points( std::size_t num_pts = 11 );

protected:

  /* ***********************  PROTECTED DATA MEMBERS  *********************** */
  double length;

  /* *********************  PROTECTED MEMBER FUNCTIONS  ********************* */

  /* Returns the stiffness matrix for the given element using the current
   * consistent tangent. */
  Eigen::MatrixXd def_stiffness( std::size_t int_order );

};

#endif
