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
  Domain( ) : nodes{ }, elements{ }, mats{ }
  { }

  /* Domain should be unique, disallow copy and assignment operators */
  Domain( const Domain & other ) = delete;
  Domain( Domain && other ) = delete;
  Domain & operator=( const Domain & rhs ) = delete;
  Domain && operator=( Domain && rhs ) = delete;

  /* Destructor */
  ~Domain( );

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  std::vector<Node *> nodes;
  std::vector<Element *> elements;
  std::vector<Material> mats;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

#endif
