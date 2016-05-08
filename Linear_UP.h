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

#ifndef GUARD_LINEAR_UP_H
#define GUARD_LINEAR_UP_H

// Project-specific headers;
#include "UP_Ele.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <vector>

class Linear_UP : public UP_Ele {

public:

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Linear_UP( ) : UP_Ele( )
  { }

  Linear_UP( std::size_t id, std::vector<Node *> nodes, const Material *mat ) :
    UP_Ele( id, nodes, 1, mat )
  {
    if( get_num_nodes( ) != 2 )
      ; // TODO:  Put an actual exception here (not sure which to use);
  }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function. */
  virtual double shape_func( double xi, std::size_t a ) const;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative. */
  virtual double shape_deriv( double, std::size_t a ) const
  {
    return ( a == 0 ) ? -0.5 : 0.5;
  }

  /* Given the parametric coordinate, xi, and the local index of the pressure
   * shape function, a, return the value of the pressure function. */
  virtual double pressure_func( double xi, std::size_t ) const
  {
    return 1.0;
  }

};

#endif
