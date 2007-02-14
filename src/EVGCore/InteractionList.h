//____________________________________________________________________________
/*!

\class   genie::InteractionList

\brief   A vector of Interaction objects.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 13, 2005

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _INTERACTION_LIST_H_
#define _INTERACTION_LIST_H_

#include <vector>
#include <ostream>

using std::vector;
using std::ostream;

namespace genie {

class Interaction;

class InteractionList : public vector<Interaction *> {

public :
  InteractionList();
  InteractionList(const InteractionList & intl);
  ~InteractionList();

  void Reset  (void);
  void Append (const InteractionList & intl);
  void Copy   (const InteractionList & intl);
  void Print  (ostream & stream) const;

  InteractionList & operator =  (const InteractionList & intl);
  friend ostream &  operator << (ostream & stream, const InteractionList & intl);
};

}      // genie namespace

#endif // _INTERACTION_LIST_H_
