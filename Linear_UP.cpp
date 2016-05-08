/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the `Linear_UP' abstraction, which   *
 * inherits publicly from `Element.' Class definition given in Linear_UP.h.   *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Linear_UP.h"

// System headers;

/* *****************************  COPY CONTROL  ***************************** */

/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the shape function.  */
double fem::Linear_UP::shape_func( double xi, std::size_t a ) const
{
  // Determine the appropriate value of xi_a;
  int xi_a = ( a == 0 ) ? -1 : 1;

  // Calculate shape function and return;
  return 0.5 * ( 1 + xi_a * xi );
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */

