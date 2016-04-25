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
void Domain::create_mat( double E, double nu )
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

/* Given the node ids and a material id, create an element and store in
 * `elements.' */
void Domain::create_elem( std::size_t n0, std::size_t n1, std::size_t mat_id )
{
  // Use current size of elements as ID of new element;
  // std::size_t ele_ID = elements.size( );
  elements.push_back( new Element( nodes[n0], nodes[n1], materials[mat_id] ) );
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
