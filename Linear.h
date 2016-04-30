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

#ifndef GUARD_LINEAR_H
#define GUARD_LINEAR_H

// Project-specific headers;
#include "Element.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <Eigen/LU>
#include <utility>

class Linear : public Element {

public:

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */

  Linear( ) : Element( )
  { }

  Linear( std::size_t id, Node *n0, Node *n1, const Material *mat ) :
    Element( id, n0, n1, mat )
  { }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

};

#endif
