
#include "Element.h"
#include "gauss_quadrature.h"
#include "Node.h"

#include <iostream>

int main( int argc, char *argv[] )
{

  // Create some nodes;
  Node *n1 = new Node( 3 );
  Node *n2 = new Node( 2 );
  Node *n3 = new Node( );

  Element e1( n1, n2 );
  e1.print_nodes( );

  delete n1;
  delete n2;
  delete n3;

  util::test( 128 );

  return 0;
}
