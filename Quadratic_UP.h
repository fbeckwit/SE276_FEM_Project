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

#ifndef GUARD_QUADRATIC_UP_H
#define GUARD_QUADRATIC_UP_H

// Project-specific headers;
#include "UP_Ele.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <vector>

namespace fem {

class Quadratic_UP : public fem::UP_Ele {

public:

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Quadratic_UP( ) : UP_Ele( )
  { }

  Quadratic_UP( std::size_t id, std::vector<Node *> nodes, const Material *mat ) :
    UP_Ele( id, nodes, 2, mat )
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

  /* Given the parametric coordinate, xi, and the local index of the pressure
   * shape function, a, return the value of the pressure function. */
  virtual double pressure_func( double xi, std::size_t a ) const;

};

} // namespace fem;

#endif
