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

  /* ****************************  COPY CONTROL  **************************** */
  /* Default constructor */
  Node( ) : coord{ 0 }, disp{ 0 }
  { }

  Node( double coord ) : coord{ coord }, disp{ 0 }
  { }

  /* Destructor */

  /* **********************  PUBLIC MEMBER FUNCTIONS  *********************** */

  inline double get_coord( ) { return coord; }

private:

  /* ************************  PRIVATE DATA MEMBERS  ************************ */

  double coord;     // Coordinate of the node;
  double disp;      // Current displacement of the node;

  /* **********************  PRIVATE MEMBER FUNCTIONS  ********************** */

};

#endif
