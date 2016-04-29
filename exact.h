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

#ifndef GUARD_EXACT_H
#define GUARD_EXACT_H

// Project-specific headers;

// System headers;

/* Function object to pass the exact displacement solution. */
struct Exact_Disp {
  Exact_Disp( double E, double nu, double P, double a, double b ) :
    E{ E }, nu{ nu }, P{ P }, a{ a }, b{ b }
  { }

  double operator()( double r ) const {
    return P * a*a * r / E / (b*b - a*a) *
      ((1 + nu)*(1 - 2*nu) + b*b/(r*r) * (1 + nu));
  }

  double E;
  double nu;
  double P;
  double a;
  double b;
};

/* Function object to pass the exact stress solution. */
struct Exact_Stress {
  Exact_Stress( double nu, double P, double a, double b ) :
    nu{ nu }, P{ P }, a{ a }, b{ b }
  { }

  Eigen::Vector3d operator()( double r ) const {
    Eigen::Vector3d ret = Eigen::Vector3d::Zero();
    ret[0] = P * a*a / (b*b - a*a) * ( 1 - b*b / (r*r) );
    ret[1] = P * a*a / (b*b - a*a) * ( 1 + b*b / (r*r) );
    ret[2] = 2 * nu * P * a*a / (b*b - a*a);
    return ret;
  }

  double nu;
  double P;
  double a;
  double b;
};

#endif
