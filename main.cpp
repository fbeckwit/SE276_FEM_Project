
#include "Element.h"
#include "gauss_quadrature.h"
#include "Node.h"

#include <iostream>
#include <cmath>

int main( int argc, char *argv[] )
{

  // Create some nodes;
  Node *n1 = new Node( 2 );
  Node *n2 = new Node( 4 );
  Node *n3 = new Node( );

  Element e1( n1, n2 );
  e1.print_nodes( );
  double xi_0 = -1.0 / std::sqrt( 3 );
  double xi_1 = -xi_0;

  std::cout << "\nB_0(xi_0) = [\n" << e1.get_gradient_matrix( xi_0, 0 ) << "\n]\n";
  std::cout << "\nB_1(xi_0) = [\n" << e1.get_gradient_matrix( xi_0, 1 ) << "\n]\n";
  std::cout << "\nB_0(xi_1) = [\n" << e1.get_gradient_matrix( xi_1, 0 ) << "\n]\n";
  std::cout << "\nB_1(xi_1) = [\n" << e1.get_gradient_matrix( xi_1, 1 ) << "\n]\n";

  delete n1;
  delete n2;
  delete n3;

  //util::test( 128 );

  return 0;
}
