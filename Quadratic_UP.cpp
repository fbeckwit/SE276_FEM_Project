/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the `Quadratic_UP' abstraction,      *
 * which inherits publicly from `Element.' Class definition given in          *
 * Quadratic_UP.h.                                                            *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Quadratic_UP.h"

// System headers;
#include <cmath>

/* *****************************  COPY CONTROL  ***************************** */

/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Given the parametric coordinate, xi, and the local index of the shape
 * function, a, return the value of the shape function.  */
double Quadratic_UP::shape_func( double xi, std::size_t a ) const
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
double Quadratic_UP::shape_deriv( double xi, std::size_t a ) const
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

/* -------------------------------------------------------------------------- */

/* Given the parametric coordinate, xi, and the local index of the pressure
 * shape function, a, return the value of the pressure function. */
double Quadratic_UP::pressure_func( double xi, std::size_t a ) const
{
  // Determine the integration location;
  double xi_a = ( a == 0 ? -1 : 1 ) / std::sqrt( 3.0 );

  // Calculate the pressure shape function;
  return 0.5 * ( 1 + 3 * xi_a * xi );
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
