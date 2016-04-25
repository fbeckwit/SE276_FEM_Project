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
#include <vector>

class Domain {

public:

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Domain( ) : nodes{ }, elements{ }, materials{ }
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
  void create_mat( double E, double nu );

  /* Given a coordinate, create a node and store in `nodes.' */
  void create_node( double coord );

  /* Given the node ids and a material id, create an element and store in
   * `elements.' */
  void create_elem( std::size_t n0, std::size_t n1, std::size_t mat_id );

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  std::vector<Node *> nodes;
  std::vector<Element *> elements;
  std::vector<Material *> materials;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

#endif
