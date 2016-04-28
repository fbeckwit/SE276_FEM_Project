
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

int main( int argc, char *argv[] )
{

  // Create some nodes;
  // Problem inputs (hard-coded for now);
  double a = 6.0;
  double b = 9.0;
  double P = 10.0;
  double E = 1000.0;
  double nu = 0.4999;
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

  Eigen::VectorXd disp = domain.solve( 1 );

  Eigen::VectorXd exact( num_elem + 1 );
  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i ) {
    double coord = a + node_i * elem_size;
    exact( node_i ) = calc_exact( coord, P, E, nu );
  }

  Eigen::VectorXd error = exact - disp;

  std::cout << "disp = [\n" << disp << "\n]\n";
  std::cout << "exact = [\n" << exact << "\n]\n";
  std::cout << "error = [\n" << error << "\n]\n";

  std::cout << '\n';
  domain.print_disp( std::cout );
  domain.print_stress( std::cout );

  return 0;
}
