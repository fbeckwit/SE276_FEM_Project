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
#include "Linear.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <Eigen/Cholesky>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class Domain {

public:

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Domain( ) : nodes{ }, elements{ }, materials{ }, num_equations{ 0 }
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

  /* Count the number of equations corresponding to free DOFs and store. */
  std::size_t get_eqn_count( );

  /* Builds the stiffness matrix by looping elements and assembling.
   * PRECONDITION:  `elements' must be properly initialized. */
  Eigen::MatrixXd build_stiffness( std::size_t int_order = 2 );

  /* Builds the force vector by looping elements and assembling.
   * PRECONDITION:  `elements' must be properly initialized. */
  Eigen::VectorXd build_force( );

  /* Builds the system of equations and then solves.
   * PRECONDITION:  `elements' must be properly initialized. */
  Eigen::VectorXd solve( std::size_t int_order = 2 );

  /* Given an output stream and the number of displacement points to print for
   * each element, compute the displacement and print to the output. */
  void print_disp( std::ostream & out = std::cout,
                     std::size_t num_pts = 11 ) const;

  /* Given an output stream and the number of stress points to print for each
   * element, compute the stress and print to the output. */
  void print_stress( std::ostream & out = std::cout,
                     std::size_t num_pts = 11 ) const;

  /* *********************  TEMPLATE MEMBER FUNCTIONS  ********************** */

  /* Given a function object to calculate the exact displacement, an output
   * stream,  and the number of displacement points to print for each element,
   * compute the displacement and print to the output. */
  template <typename Func>
  void print_disp( const Func & exact,
      std::ostream & out = std::cout,
      std::size_t num_pts = 11 ) const
  {
    // Set output precision and format;
    std::streamsize prec = out.precision( 6 );
    std::streamsize width = 14;
    out << std::scientific;

    // Print header;
    out << '#' << std::string( width - 1, ' ' );
    print_centered( std::string( "Radial Displacement" ), width * 3, out, true);
    out << "\n#" << std::setw( width - 1 ) << "Radius:";
    out << std::setw( width ) << "FEM:" << std::setw( width ) << "Exact:" <<
      std::setw( width ) << "Error:";
    out << '\n';

    // Loop over the elements and plotting points, grab their disp, and output;
    for( const auto elem : elements ) {
      elem->print_disp( exact, width, out, num_pts );
    }
    // Reset precision and format;
    out << std::fixed << std::setprecision( prec );
  }

  /* Given a function object to calculate the exact stresses, an output stream,
   * and the number of stress points to print for each element, compute the
   * stress and print to the output. */
  template <typename Func>
  void print_stress( const Func & exact,
      std::ostream & out = std::cout,
      std::size_t num_pts = 11 ) const
  {
    // Set output precision and format;
    std::streamsize prec = out.precision( 6 );
    std::size_t width = 14;
    out << std::scientific;

    // Print header;
    out << '#' << std::string( width - 1, ' ' );
    print_centered( std::string( "Radial Stress" ), width * 3, out);
    print_centered( std::string( "Hoop Stress" ), width * 3, out);
    print_centered( std::string( "Axial Stress" ), width * 3, out, true);
    out << "\n#" << std::setw( width - 1 ) << "Radius:";
    for( std::size_t j{ 0 }; j != 3; ++j )
      out << std::setw( width ) << "FEM:" << std::setw( width ) << "Exact:" <<
        std::setw( width ) << "Error:";
    out << '\n';

    // Loop over the elements and plotting points, grab their disp, and output;
    for( const auto elem : elements ) {
      elem->print_stress( exact, width, out, num_pts );
    }
    // Reset precision and format;
    out << std::fixed << std::setprecision( prec );
  }

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  /* Domain elements */
  std::vector<Node *> nodes;
  std::vector<Element *> elements;
  std::vector<Material *> materials;
  std::size_t num_equations;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

  /* Given a vector of displacements, update the nodes. */
  void update_nodes( const Eigen::VectorXd & displacement );

  /* Given a string, a field width, and an output stream, center the string and
   * print to the output. */
  void print_centered( const std::string & str,
      std::size_t width, std::ostream & out, bool trim = false ) const;

};

#endif
