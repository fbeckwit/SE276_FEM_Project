/* ************************************************************************** *
 *                           Frank Nathan Beckwith                            *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 * ************************************************************************** *
 *                                                                            *
 * Source file for the implementation of the Domain abstraction.              *
 * Class definition given in Domain.h.                                        *
 *                                                                            *
 * ************************************************************************** */

// Project-specific headers;
#include "Domain.h"
#include "Element.h"
#include "exact.h"
#include "gauss_quadrature.h"
#include "Material.h"
#include "Node.h"

// System headers;
#include <cmath>
#include <cstdlib>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

/* ****************************  BEGIN PROGRAM  ***************************** */
int main( int argc, char *argv[] )
{
  // Make usage string;
  std::string usage( "Usage: ");
  usage += argv[0];
  usage += " poisson_ratio num_elements [disp_outfile stress_outfile]\n";

  // Check for appropriate number of inputs;
  if( argc < 3 ) {
    std::cerr << "ERROR:  Incorrect number of inputs.\n";
    std::cerr << usage;
    return -1;
  }

  // Check optional arguments, require user to provide two output file names;
  if( argc > 3 && argc < 5 ) {
    std::cerr << "Error:  Please give second output file.\n";
    std::cerr << usage;
    return -1;
  }

  // Grab inputs;
  double nu = std::atof( argv[1] );
  std::size_t num_elem = atoi( argv[2] );

  // Problem inputs (hard-coded for now);
  double a = 6.0;
  double b = 9.0;
  double P = 10.0;
  double E = 1000.0;

  // Print message;
  std::cout << "Problem Inputs:\n";
  std::cout << "    a = " << a << " in\n";
  std::cout << "    b = " << b << " in\n";
  std::cout << "    P = " << P << " psi\n";
  std::cout << "    E = " << E << " psi\n";
  std::cout << "    nu = " << nu << "\n";
  std::cout << "    No. Elements = " << num_elem << "\n";

  // Create domain nodes;
  std::cout << "\nCreating domain:\n    Nodes ...\n";
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
  std::cout << "    Materials ...\n";
  domain.create_material( E, nu );

  // Create elements;
  std::cout << "    Elements ...\n";
  for( std::size_t ele_i{ 0 }; ele_i != num_elem; ++ele_i ) {
    domain.create_element( ele_i, ele_i + 1, 0 );
  }

  // Solve system of equations;
  std::cout << "\nSolving system of equations:\n";
  Eigen::VectorXd disp = domain.solve( );

  // Output results with comparison to anayltical;
  Exact_Disp disp_func{ E, nu, P, a, b };
  Exact_Stress stress_func{ nu, P, a, b };

  // If optional arguments given, print to files.  Else, print to console.
  if( argc > 3 ) {
    std::ofstream out_disp( argv[3] );
    std::ofstream out_stress( argv[4] );

    domain.print_disp( disp_func, out_disp);
    domain.print_stress( stress_func, out_stress );

    out_disp.close( );
    out_stress.close( );
  }
  else {
    std::cout << '\n';
    domain.print_disp( disp_func );
    std::cout << '\n';
    domain.print_stress( stress_func );
  }

  return 0;
}
