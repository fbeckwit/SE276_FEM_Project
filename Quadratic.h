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

#ifndef GUARD_QUADRATIC_H
#define GUARD_QUADRATIC_H

// Project-specific headers;
#include "Element.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <iostream>
#include <Eigen/LU>
#include <vector>

class Quadratic : public Element {

public:

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Quadratic( ) : Element( )
  { }

  Quadratic( std::size_t id, std::vector<Node *> nodes, const Material *mat ) :
    Element( id, nodes, mat )
  {
    if( get_num_nodes( ) != 3 )
      ; // TODO:  Put an actual exception here (not sure which to use);
  }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function. */
  virtual double shape_func( double xi, std::size_t a ) const;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative. */
  virtual double shape_deriv( double xi, std::size_t a ) const;

};

#endif
