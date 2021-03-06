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

namespace fem {

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

  /* Return the tangent elastic modulii tensor for the material. */
  Eigen::Matrix2d get_tangent( ) const;

  /* Given the strain, return the resulting stress. */
  Eigen::Vector3d get_stress( const Eigen::Vector2d & strain ) const;

  /* Given the strain, return the resulting stress. */
  Eigen::Vector3d get_stress_mu( const Eigen::Vector2d & strain ) const;

  /* Return the Lamé constants, lambda or mu. */
  double get_lambda( ) const { return lambda; }
  double get_mu( ) const { return mu; }

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  double lambda;
  double mu;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

} // namespace fem;

#endif
