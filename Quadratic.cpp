/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the `Quadratic' abstraction, which      *
 * inherits publicly from `Element.' Class definition given in Quadratic.h.      *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Quadratic.h"
#include "gauss_quadrature.h"

// System headers;

/* *****************************  COPY CONTROL  ***************************** */

/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the shape function.  */
double Quadratic::shape_func( double xi, std::size_t a ) const
{
  if( a == 0 )
    return 0.5 * xi * ( xi - 1 );
  else if (a == 1 )
    return 1 - xi * xi;
  else if (a == 2 )
    return 0.5 * xi * ( xi + 1 );
  else
    return -99999;
}

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the shape function derivative. */
double Quadratic::shape_deriv( double xi, std::size_t a ) const
{
  if( a == 0 )
    return xi - 0.5;
  else if ( a == 1 )
    return -2 * xi;
  else if ( a == 2 )
    return xi + 0.5;
  else
    return -99999;
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
