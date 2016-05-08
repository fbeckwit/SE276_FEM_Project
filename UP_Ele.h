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

#ifndef GUARD_UP_ELE_H
#define GUARD_UP_ELE_H

// Project-specific headers;
#include "Element.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cstddef>
#include <Eigen/LU>
#include <utility>

namespace fem {

class UP_Ele : public fem::Element {

public:

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */
  UP_Ele( ) :
    Element( ), pressure( ), material( nullptr ),
    k_eval( this ), g_eval( this ), m_eval( this )
  { }

  UP_Ele( std::size_t id, std::vector<Node *> nodes,
          std::size_t num_pres, const Material *mat ) :
    Element( id, nodes ), pressure( num_pres, 0.0 ), material( mat->clone( ) ),
    k_eval( this ), g_eval( this ), m_eval( this )
  { }

  UP_Ele( const UP_Ele & other ) :
    Element( other ), pressure{ other.pressure },
    material{ other.material->clone( ) },
    k_eval( this ), g_eval( this ), m_eval( this )
  { }

  UP_Ele( UP_Ele && other ) :
    Element( std::move( other ) ), pressure{ std::move( other.pressure ) },
    material{ other.material },
    k_eval( this ), g_eval( this ), m_eval( this )
  {
    other.material = nullptr;
  }

  /* Deleted assignment operators */
  UP_Ele & operator=( const UP_Ele & ) = delete;
  UP_Ele && operator=( UP_Ele && ) = delete;

  virtual ~UP_Ele( ) { delete material; }

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  /* Returns the stiffness matrix for the given element using the current
   * consistent tangent. */
  Eigen::MatrixXd get_stiffness( std::size_t int_order );

  /* Given the parametric coordinate, xi, interpolate the stresses from the
   * resulting displacement.
   * PRECONDITION:  Nodes must have updated displacements. */
  Eigen::Vector3d interp_stress( double xi ) const;

  /* Given the parametric coordinate, xi, and the node number, a, return the
   * divergence matrix, b^v. */
  double get_divergence_matrix( double xi, std::size_t a ) const;

  /* Update the element info.
   * PRECONDITION:  Element nodes must be updated. */
  void update( );

  /* Given the parametrix coordinate, xi, interpolate the pressure.
   * PRECONDITION:  Pressures must be updated after solving. */
  double interp_pressure( double xi ) const;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function. */
  virtual double shape_func( double xi, std::size_t a ) const = 0;

  /* Given the parametric coordinate, xi, and the local index of the shape
   * function, a, return the value of the shape function derivative. */
  virtual double shape_deriv( double xi, std::size_t a ) const = 0;

  /* Given the parametric coordinate, xi, and the local index of the pressure
   * shape function, a, return the value of the pressure function. */
  virtual double pressure_func( double xi, std::size_t a ) const = 0;

protected:

  /* ***********************  PROTECTED DATA MEMBERS  *********************** */

  std::vector<double> pressure;
  Material *material;

private:

  /* ***************************  NESTED CLASSES  *************************** */

  /* Function object used in the evaluation of the stiffness matrix.  operator()
   * overloaded to return the internal energy density (from the shear modulus)
   * at the evaluation point, xi. */
  struct K_Func {

    /* Constructor */
    K_Func( const UP_Ele * p ) : parent{ p }, a{ 0 }, b{ 0 } { }

    /* Functions to query the size of the final matrix. */
    std::size_t get_rows( ) const { return parent->nodes.size( ); }
    std::size_t get_cols( ) const { return parent->nodes.size( ); }

    /* Calculate the internal energy density of the stiffness. */
    double operator()( double xi ) const;

    const UP_Ele * parent;
    std::size_t a;
    std::size_t b;
  };

  /* Function object used in the evaluation of the stiffness matrix.  operator()
   * overloaded to return the internal energy density due to dilation at the
   * evaluation point, xi. */
  struct G_Func {
    /* Constructor */
    G_Func( const UP_Ele * p ) : parent{ p }, a{ 0 }, b{ 0 } { }

    /* Functions to query the size of the final matrix. */
    std::size_t get_rows( ) const { return parent->nodes.size( ); }
    std::size_t get_cols( ) const { return parent->pressure.size( ); }

    /* Calculate the internal energy density of the dilation. */
    double operator()( double xi ) const;

    const UP_Ele * parent;
    std::size_t a;
    std::size_t b;
  };

  /* Function object used in the evaluation of the stiffness matrix.  operator()
   * overloaded to return the internal energy density due to the penalty
   * constraint at the evaluation point, xi. */
  struct M_Func {

    /* Constructor */
    M_Func( const UP_Ele * p ) : parent{ p }, a{ 0 }, b{ 0 } { }

    /* Functions to query the size of the final matrix. */
    std::size_t get_rows( ) const { return parent->pressure.size( ); }
    std::size_t get_cols( ) const { return parent->pressure.size( ); }

    /* Calculate the internal energy density of the penalty constraint. */
    double operator()( double xi ) const;

    const UP_Ele * parent;
    std::size_t a;
    std::size_t b;
  };

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  K_Func k_eval;
  G_Func g_eval;
  M_Func m_eval;

  /* *********************  PRIVATE MEMBERS FUNCTIONS  ********************** */

};

} // namespace fem;

#endif
