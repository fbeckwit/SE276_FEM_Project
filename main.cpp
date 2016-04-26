
#include "Domain.h"
#include "Element.h"
#include "gauss_quadrature.h"
#include "Material.h"
#include "Node.h"

#include <Eigen/LU>
#include <iostream>
#include <cmath>

void print_stiffness(
    Element & e,
    std::size_t int_order,
    std::ostream &out = std::cout
    )
{
  out << "Stiffness(" << int_order << "-pt) = [\n";
  out << e.get_stiffness( int_order ) << "\n]\n";
}

int main( int argc, char *argv[] )
{

  std::size_t num_elem = 1;
  // Create some nodes;
  std::vector<Node *> nodes;
  std::vector<Element> elems;

  // Problem inputs (hard-coded for now);
  double a = 6.0;
  double b = 9.0;
  double P = 10.0;
  double E = 1000.0;
  double nu = 0.49999;

  double elem_size = ( b - a ) / num_elem;
  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i ) {
    if( node_i == 0 )
      nodes.push_back( new Node( node_i,  a + node_i * elem_size, Node::NBC, 10.0 ) );
    else
      nodes.push_back( new Node( node_i,  a + node_i * elem_size ));
  }

  Material mat1( E, nu );
  std::cout << "mat1(" << &mat1 << ") = [\n" << mat1.get_tangent( ) << "\n]\n";

  for( std::size_t ele{ 0 }; ele != num_elem; ++ele )
    elems.push_back( Element( ele, nodes[ ele ], nodes[ ele + 1 ], &mat1 ) );

  Eigen::MatrixXd stiff = elems[0].get_stiffness( 1 );
  Eigen::VectorXd force = elems[0].get_force_ext( );

  Eigen::VectorXd disp = stiff.inverse( ) * force;
  Eigen::VectorXd exact( 2 );
  exact( 0 ) = P * a*a * a / E / (b*b - a*a) * ((1 + nu)*(1 - 2*nu) + b*b/(a*a)*(1 + nu));
  exact( 1 ) = P * a*a * b / E / (b*b - a*a) * ((1 + nu)*(1 - 2*nu) + b*b/(b*b)*(1 + nu));

  Eigen::VectorXd error = exact - disp;

  std::cout << "disp = [\n" << disp << "\n]\n";
  std::cout << "exact = [\n" << exact << "\n]\n";
  std::cout << "error = [\n" << error << "\n]\n";

  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i )
    delete nodes[node_i];

  return 0;
}
