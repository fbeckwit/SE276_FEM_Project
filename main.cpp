
#include "Element.h"
#include "gauss_quadrature.h"
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

  double a = 6.0;
  double b = 9.0;
  double elem_size = ( b - a ) / num_elem;
  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i ) {
    if( node_i == 0 )
      nodes.push_back( new Node( a + node_i * elem_size, Node::NBC, 10.0 ) );
    else
      nodes.push_back( new Node( a + node_i * elem_size ));
  }

  for( std::size_t ele{ 0 }; ele != num_elem; ++ele )
    elems.push_back( Element( nodes[ ele ], nodes[ ele + 1 ] ) );

  print_stiffness( elems[0], 1 );
  print_stiffness( elems[0], 2 );
  print_stiffness( elems[0], 3 );
  print_stiffness( elems[0], 4 );

  Eigen::MatrixXd stiff = elems[0].get_stiffness( 1 );
  Eigen::VectorXd force = elems[0].get_force_ext( );

  Eigen::VectorXd disp = stiff.inverse( ) * force;
  std::cout << disp << "\n";

  for( std::size_t node_i{ 0 }; node_i != num_elem + 1; ++node_i )
    delete nodes[node_i];

  return 0;
}
