/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the Material abstraction.            *
 * Class definition given in Material.h.                                      *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Material.h"

// System headers;

/* *****************************  COPY CONTROL  ***************************** */

/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* Return the tangent elastic modulii tensor for the material. */
Eigen::Matrix2d Material::get_tangent( ) const
{
  // Create the tangent matrix, load the entries and return;
  Eigen::Matrix2d elastic_mod;
  elastic_mod( 0, 0 ) = elastic_mod( 1, 1 ) = lambda + 2 * mu;
  elastic_mod( 0, 1 ) = elastic_mod( 1, 0 ) = lambda;
  return elastic_mod;
}

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
