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

// System headers;

/* *****************************  COPY CONTROL  ***************************** */

/* Destructor */
Domain::~Domain( )
{
  // Delete the nodes;
  for( std::vector<Node *>::iterator it = nodes.begin( );
      it != nodes.end( ); ++it )
    delete *it;

  // Delete the elements;
  for( std::vector<Element *>::iterator it = elements.begin( );
      it != elements.end( ); ++it )
    delete *it;
}

/* ***********************  PUBLIC MEMBER FUNCTIONS  ************************ */

/* -------------------------------------------------------------------------- */

/* ***********************  PRIVATE MEMBER FUNCTIONS  *********************** */
