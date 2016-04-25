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

#ifndef GUARD_MATERIAL_H
#define GUARD_MATERIAL_H

// Project-specific headers;

// System headers;
#include <Eigen/LU>
#include <iostream>

class Material {

public:

  /* ****************************  ENUMERATIONS  **************************** */

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Material( double E = 0, double nu = 0 )
  {
    lambda = nu * E / (( 1 + nu ) * ( 1 - 2 * nu ));
    mu = E / 2 / ( 1 + nu );
  }

  Material * clone( ) const { return new Material( *this ); }

  /* Destructor */

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* get_tangent( )
   * Return the tangent elastic modulii tensor for the material.
   */
  Eigen::Matrix2d get_tangent( ) const;

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  double lambda;
  double mu;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

#endif

