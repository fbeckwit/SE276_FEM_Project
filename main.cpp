
#include "Domain.h"
#include "Element.h"
#include "gauss_quadrature.h"
#include "Material.h"
#include "Node.h"

#include <Eigen/Dense>
#include <iostream>
#include <cmath>

double calc_exact( double r, double P, double E, double nu )
{
  double a = 6.0;
  double b = 9.0;
  return P * a*a * r / E / (b*b - a*a) *
    ((1 + nu)*(1 - 2*nu) + b*b/(r*r) * (1 + nu));
}

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

int main( int argc, char *argv[] )
{

  // Create some nodes;
  // Problem inputs (hard-coded for now);
  double a = 6.0;
  double b = 9.0;
  double P = 10.0;
  double E = 1000.0;
  double nu = 0.25;
  std::size_t num_elem = 10;

  // Load nodes;
  Domain domain;
  double elem_size = ( b - a ) / num_elem;
  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i ) {
    double coord = a + node_i * elem_size;
    if( node_i == 0 )
      domain.create_node( coord, Node::NBC, P );
    else
      domain.create_node( coord );
  }

  // Create material;
  domain.create_material( E, nu );

  for( std::size_t ele_i{ 0 }; ele_i != num_elem; ++ele_i ) {
    domain.create_element( ele_i, ele_i + 1, 0 );
  }

  Eigen::VectorXd disp = domain.solve( 2 );

  Eigen::VectorXd exact( num_elem + 1 );
  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i ) {
    double coord = a + node_i * elem_size;
    exact( node_i ) = calc_exact( coord, P, E, nu );
  }


  Exact_Disp disp_func{ E, nu, P, a, b };
  Exact_Stress stress_func{ nu, P, a, b };
  std::cout << '\n';
  domain.print_disp( disp_func );
  domain.print_stress( stress_func );

  return 0;
}
