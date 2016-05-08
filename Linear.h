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
#include "Disp_Ele.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <vector>

namespace fem {

class Linear : public fem::Disp_Ele {

public:

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Linear( ) : Disp_Ele( )
  { }

  Linear( std::size_t id, std::vector<Node *> nodes, const Material *mat ) :
    Disp_Ele( id, nodes, mat )
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

};

} // namespace fem;

#endif
