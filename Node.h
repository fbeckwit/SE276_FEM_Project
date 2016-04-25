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

#ifndef GUARD_NODE_H
#define GUARD_NODE_H

// Project-specific headers;

// System headers;
#include <array>
#include <vector>

class Node {

public:

  /* ****************************  ENUMERATIONS  **************************** */

  /* Enumeration to indicate if node is interior, or if on natural or essential
   * boundary.
   */
  enum node_type { INT, EBC, NBC };

  /* ****************************  COPY CONTROL  **************************** */

  /* Default constructor */
  Node( double coord = 0, node_type type = INT, double bc = 0 ) :
    coord{ coord }, type{ type }, disp{ 0 }, bound_cond{ bc }
  { }

  /* Destructor */

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  inline double get_coord( ) { return coord; }

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  double coord;       // Coordinate of the node;
  node_type type;     // Marker indicating interior or boundary;
  double disp;        // Current displacement of the node;
  double bound_cond;  // Boundary condition if node is not interior;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

#endif
