/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the Domain abstraction.              *
 * Class definition given in Domain.h.                                        *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Domain.h"

// System headers;
#include <iomanip>

/* *****************************  COPY CONTROL  ***************************** */

/* Destructor */
Domain::~Domain( )
{
  // Delete the nodes, elements, and materials;
  for( auto it : nodes )
    delete it;
  for( auto it : elements )
    delete it;
  for( auto it : materials )
    delete it;
}

/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Given material properties, Young's modulus and Poisson's ratio, create an
 * elastic material and store in `mats.' */
void Domain::create_material( double E, double nu )
{
  materials.push_back( new Material( E, nu ) );
}

/* -------------------------------------------------------------------------- */

/* Given a coordinate, create a node and store in `nodes.' */
void Domain::create_node( double coord )
{
  // Use current size of nodes as ID of new node;
  std::size_t node_ID = nodes.size( );
  nodes.push_back( new Node( node_ID, coord ) );
}

/* -------------------------------------------------------------------------- */

/* Given a coordinate, node type, and BC, create a node and store in `nodes.' */
void Domain::create_node( double coord, Node::node_type type, double bc )
{
  // Use current size of nodes as ID of new node;
  std::size_t node_ID = nodes.size( );
  nodes.push_back( new Node( node_ID, coord, type, bc ) );
}

/* -------------------------------------------------------------------------- */

/* Given the node ids and a material id, create an element and store in
 * `elements.'
 * PRECONDITION:  Nodes `n0' and `n1' and material `mat_id' must be created */
void Domain::create_element(
    std::size_t n0,
    std::size_t n1,
    std::size_t mat_id
    )
{
  // Use current size of elements as ID of new element;
  std::size_t ele_ID = elements.size( );

  // Get the element nodes;
  std::vector<Node *> node_ele = { nodes[n0], nodes[n1] };
  elements.push_back(
      new Linear( ele_ID, node_ele, materials[mat_id] )
      );
}

/* -------------------------------------------------------------------------- */

/* Count the number of equations corresponding to free DOFs and store. */
std::size_t Domain::get_eqn_count( )
{
  // Reset the DOF counter;
  num_equations = 0;

  // Loop over node DOFs and count those not on the EBC;
  for( auto node : nodes ) {
    if( node->get_type( ) != Node::EBC )
      ++num_equations;
  }
  return num_equations;
}

/* -------------------------------------------------------------------------- */

/* Builds the stiffness matrix by looping elements and assembling.
 * PRECONDITION:  `elements' must be properly initialized. */
Eigen::MatrixXd Domain::build_stiffness( std::size_t int_order )
{
  // Resize the stiffness matrix;
  Eigen::MatrixXd stiff( num_equations, num_equations );
  stiff.setZero( );

  // Loop over elements, get each stiffness and assemble to global stiffness;
  for( std::vector<Element *>::const_iterator elem_it = elements.begin( );
      elem_it != elements.end( ); ++elem_it ) {
    Element * elem = *elem_it;

    Eigen::MatrixXd stiff_elem = elem->get_stiffness( int_order );
    for( std::size_t a{ 0 }; a != elem->num_nodes( ); ++a ) {
      for( std::size_t b{ 0 }; b != elem->num_nodes( ); ++b ) {

        // Check if node is free or not;
        if( elem->get_node_type( a ) != Node::EBC &&
            elem->get_node_type( b ) != Node::EBC ) {

          // Get the global index numbers and assemble component to global;
          std::size_t A = elem->location_matrix( a );
          std::size_t B = elem->location_matrix( b );
          stiff( A, B ) += stiff_elem( a, b );
        }
      }
    }
  }
  return stiff;
}

/* -------------------------------------------------------------------------- */

/* Builds the force vector by looping elements and assembling.
 * PRECONDITION:  `elements' must be properly initialized and num_equations must
 * be valid. */
Eigen::VectorXd Domain::build_force( )
{
  // Resize the force vector;
  Eigen::VectorXd force( num_equations );
  force.setZero( );

  // Loop over the elements, get each external force and assemble to global;
  for( std::vector<Element *>::const_iterator elem_it = elements.begin( );
      elem_it != elements.end( ); ++elem_it ) {
    Element * elem = *elem_it;

    Eigen::VectorXd force_elem = elem->get_force_ext( );
    for( std::size_t a{ 0 }; a != elem->num_nodes( ); ++a ) {

      // Check if node is free or not;
      if( elem->get_node_type( a ) != Node::EBC ) {

        // Get the global index number and assemble component to global;
        std::size_t A = elem->location_matrix( a );
        force( A ) += force_elem( a );
      }
    }
  }
  return force;
}

/* -------------------------------------------------------------------------- */

/* Given the integration order, builds the system of equations, solves, and
 * returns displacement.
 * PRECONDITION:  `elements' must be properly initialized. */
Eigen::VectorXd Domain::solve( std::size_t int_order )
{
  // Get the number of equations of the system;
  get_eqn_count( );

  // Build the stiffness and force vectors;
  Eigen::MatrixXd stiff = build_stiffness( int_order );
  Eigen::VectorXd force = build_force( );

  // Solve system;
  Eigen::VectorXd disp( num_equations );
  disp.setZero( );
  disp = stiff.llt( ).solve( force );

  // Update the nodes;
  update_nodes( disp );
  return disp;
}

/* -------------------------------------------------------------------------- */

/* Given an output stream and the number of displacement points to print for
 * each element, compute the displacement and print to the output. */
void Domain::print_disp( std::ostream & out, std::size_t num_pts ) const
{
  // Set output precision and format;
  std::streamsize prec = out.precision( 6 );
  std::streamsize width = 14;
  out << std::scientific;

  // Loop over the elements and plotting points, grab their disp, and output;
  for( const auto elem : elements ) {
    // Get the radius and displacements over the element;
    std::vector<double> radius = elem->interp_coord( num_pts );
    std::vector<double> disp = elem->interp_disp( num_pts );

    for( std::vector<double>::size_type i{ 0 }; i != disp.size( ); ++i ) {
      // Output the information;
      out << std::setw( width ) << radius[i];
      out << std::setw( width ) << disp[i];
      out << '\n';
    }
  }
  // Reset precision and format;
  out << std::fixed << std::setprecision( prec );
}

/* -------------------------------------------------------------------------- */

/* Given an output stream and the number of stress points to print for each
 * element, compute the stress and print to the output. */
void Domain::print_stress( std::ostream & out, std::size_t num_pts ) const
{
  // Set output precision and format;
  std::streamsize prec = out.precision( 6 );
  std::streamsize width = 14;
  out << std::scientific;

  // Loop over the elements and plotting points, grab their stress, and output;
  for( const auto elem : elements ) {
    std::vector<double> radius = elem->interp_coord( num_pts );
    std::vector<Eigen::Vector3d> stress = elem->interp_stress( num_pts );

    for( std::vector<double>::size_type i{ 0 }; i != radius.size( ); ++i ) {
      // Output the information;
      out << std::setw( width ) << radius[i];
      for( std::size_t j{ 0 }; j != 3; ++j )
        out << std::setw( width ) << stress[i][j];
      out << '\n';
    }
  }
  // Reset precision and format;
  out << std::fixed << std::setprecision( prec );
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */

/* Given a vector of displacements, update the nodes. */
void Domain::update_nodes( const Eigen::VectorXd & displacement )
{
  // Loop over the nodes and update any necessary information;
  for( std::vector<Node *>::iterator node_iter = nodes.begin( );
      node_iter != nodes.end( ); ++node_iter ) {
    Node *node = *node_iter;

    // Check if the node is a DOF of the system;
    if( node->get_type( ) != Node::EBC ) {
      // Grab the global equation number and update the dipslacement;
      std::size_t A = node->get_eqn_num( );
      node->update_disp( displacement[A] );
    }
  }
}

/* -------------------------------------------------------------------------- */

/* Given a string, a field width, and an output stream, center the string and
 * print to the output. */
void Domain::print_centered( const std::string & str,
    std::size_t width, std::ostream & out, bool trim ) const
{
  // Create the header string;
  int num_blanks = ( width - str.size( ) ) / 2;
  std::string header = std::string( num_blanks, ' ' ) + str;
  if( !trim )
    header += std::string( width - header.size( ), ' ' );

  // Print to output;
  out << header;
}
