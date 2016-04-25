
#include "Element.h"
#include "gauss_quadrature.h"
#include "Node.h"

#include <iostream>
#include <cmath>

void print_stiffness( Element & e, std::size_t int_order, std::ostream &out = std::cout )
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

  double a = 6.0;
  double b = 9.0;
  double elem_size = ( b - a ) / num_elem;
  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i ) {
    nodes.push_back( new Node( a + node_i * elem_size ));

    if( node_i > 0 )
      elems.push_back( Element( nodes[ node_i - 1 ], nodes[ node_i ] ));
  }


  double xi_0 = 0;
  double xi_1 = -xi_0;

  print_stiffness( elems[0], 1 );
  print_stiffness( elems[0], 2 );
  print_stiffness( elems[0], 3 );
  print_stiffness( elems[0], 4 );

  std::cout << "\nForce(e1) = [\n" << elems[0].get_force_ext( ) << "\n]\n";

  std::cout << "\nB_0(xi_0) = [\n" << elems[0].get_gradient_matrix( xi_0, 0 ) << "\n]\n";
  std::cout << "\nB_1(xi_0) = [\n" << elems[0].get_gradient_matrix( xi_0, 1 ) << "\n]\n";
  std::cout << "\nB_0(xi_1) = [\n" << elems[0].get_gradient_matrix( xi_1, 0 ) << "\n]\n";
  std::cout << "\nB_1(xi_1) = [\n" << elems[0].get_gradient_matrix( xi_1, 1 ) << "\n]\n";

  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i )
    delete nodes[node_i];

  return 0;
}
