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

/* *****************************  COPY CONTROL  ***************************** */

/* Destructor */
Domain::~Domain( )
{
  // Delete the nodes;
  for( std::vector<Node *>::iterator it = nodes.begin( );
      it != nodes.end( ); ++it )
    delete *it;

  // Delete the elements;
  for( std::vector<Element *>::iterator it = elements.begin( );
      it != elements.end( ); ++it )
    delete *it;

  // Delete the materials;
  for( std::vector<Material *>::iterator it = materials.begin( );
      it != materials.end( ); ++it )
    delete *it;
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
  elements.push_back(
      new Element( ele_ID, nodes[n0], nodes[n1], materials[mat_id] )
      );
}

/* -------------------------------------------------------------------------- */

/* Builds the stiffness matrix by looping elements and assembling.
 * PRECONDITION:  `elements' must be properly initialized. */
void Domain::build_stiffness( std::size_t int_order )
{
  // Get number of equations;
  // TODO:  Generalize to problems with essential boundary conditions.
  std::size_t NEQ = nodes.size( );

  // Resize the stiffness matrix;
  stiff.resize( NEQ, NEQ );

  // Loop over elements, get each stiffness and assemble to global stiffness;
  for( std::vector<Element *>::const_iterator elem_it = elements.begin( );
      elem_it != elements.end( ); ++elem_it ) {
    Element * elem = *elem_it;

    Eigen::MatrixXd stiff_elem = elem->get_stiffness( int_order );
    for( std::size_t a{ 0 }; a != elem->NEN; ++a ) {
      for( std::size_t b{ 0 }; b != elem->NEN; ++b ) {

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
}

/* -------------------------------------------------------------------------- */

/* Builds the force vector by looping elements and assembling.
 * PRECONDITION:  `elements' must be properly initialized. */
void Domain::build_force( )
{
  // Get number of equations;
  // TODO:  Generalize to problems with essential boundary conditions.  Also
  // make a private function (needed in the above as well).
  std::size_t NEQ = nodes.size( );

  // Resize the force vector;
  force.resize( NEQ );

  // Loop over the elements, get each external force and assemble to global;
  for( std::vector<Element *>::const_iterator elem_it = elements.begin( );
      elem_it != elements.end( ); ++elem_it ) {
    Element * elem = *elem_it;

    Eigen::VectorXd force_elem = elem->get_force_ext( );
    for( std::size_t a{ 0 }; a != elem->NEN; ++a ) {

      // Check if node is free or not;
      if( elem->get_node_type( a ) != Node::EBC ) {

        // Get the global index number and assemble component to global;
        std::size_t A = elem->location_matrix( a );
        force( A ) += force_elem( a );
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

/* Builds the system of equations and then solves.
 * PRECONDITION:  `elements' must be properly initialized. */
Eigen::VectorXd Domain::solve( std::size_t int_order )
{
  // Build the stiffness and force vectors;
  build_stiffness( int_order );
  build_force( );

  disp = stiff.inverse( ) * force;
  return disp;
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
