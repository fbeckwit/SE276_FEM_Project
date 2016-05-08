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

#ifndef GUARD_DISP_ELE_H
#define GUARD_DISP_ELE_H

// Project-specific headers;
#include "Element.h"
#include "gauss_quadrature.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <Eigen/LU>
#include <utility>

class Disp_Ele : public Element {

public:

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */
  Disp_Ele( ) : Element( ), material( nullptr ), stiff_eval( this ) { }

  Disp_Ele( std::size_t id, std::vector<Node *> nodes, const Material *mat ) :
    Element( id, nodes ), material( mat->clone( ) ), stiff_eval( this )
  { }

  Disp_Ele( const Disp_Ele & other ) :
    Element( other ), material{ other.material->clone( ) },
    stiff_eval( this ) { }

  Disp_Ele( Disp_Ele && other ) :
    Element( std::move( other ) ), material{ other.material }, stiff_eval( this )
  {
    other.material = nullptr;
  }

  /* Deleted assignment operators */
  Disp_Ele & operator=( const Disp_Ele & ) = delete;
  Disp_Ele && operator=( Disp_Ele && ) = delete;

  virtual ~Disp_Ele( ) { delete material; }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* Returns the stiffness matrix for the given element using the current
   * consistent tangent. */
  Eigen::MatrixXd get_stiffness( std::size_t int_order )
  {
    // Call function from quadrature namespace with `sym = true';
    return util::integrate_matrix( stiff_eval, int_order, true );
  }

  /* Given the parametric coordinate, xi, interpolate the stresses from the
   * resulting displacement.
   * PRECONDITION:  Nodes must have updated displacements. */
  Eigen::Vector3d interp_stress( double xi ) const;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function. */
  virtual double shape_func( double xi, std::size_t a ) const = 0;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative. */
  virtual double shape_deriv( double xi, std::size_t a ) const = 0;

protected:

  /* ***********************  PROTECTED DATA MEMBERS  *********************** */

  Material *material;

private:

  /* ***************************  NESTED CLASSES  *************************** */

  /* Function object used in the evaluation of the stiffness matrix.  operator()
   * overloaded to return the internal energy density at the evaluation point,
   * xi. */
  struct K_Func {

    /* Constructor */
    K_Func( const Disp_Ele * p ) : parent{ p }, a{ 0 }, b{ 0 } { }

    /* Functions to query the size of the final matrix. */
    std::size_t get_rows( ) const { return parent->nodes.size( ); }
    std::size_t get_cols( ) const { return parent->nodes.size( ); }

    /* Given a parametric coordinate, xi, calculate the internal energy density
     * of the stiffness. */
    double operator()( double xi ) const;

    const Disp_Ele * parent;
    std::size_t a;
    std::size_t b;
  };

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  K_Func stiff_eval;

};

#endif
