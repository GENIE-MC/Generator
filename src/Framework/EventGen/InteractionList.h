//____________________________________________________________________________
/*!

\class   genie::InteractionList

\brief   A vector of Interaction objects.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 13, 2005

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
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
class InteractionList;

ostream & operator << (ostream & stream, const InteractionList & intl);

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
