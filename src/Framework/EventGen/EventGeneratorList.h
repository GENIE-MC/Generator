//____________________________________________________________________________
/*!

\class   genie::EventGeneratorList

\brief   A vector of EventGeneratorI objects

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created January 25, 2004

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org         
*/
//____________________________________________________________________________

#ifndef _EVENT_GENERATOR_LIST_H_
#define _EVENT_GENERATOR_LIST_H_

#include <vector>
#include <ostream>

using std::vector;
using std::ostream;

namespace genie {

class EventGeneratorList;
class EventGeneratorI;

ostream & operator << (ostream & stream, const EventGeneratorList & evgl);

class EventGeneratorList : public vector<const EventGeneratorI *> {

public :

  EventGeneratorList();
  ~EventGeneratorList();

  void Print(ostream & stream) const;

  friend ostream & operator << (ostream & stream, const EventGeneratorList & evgl);
};

}      // genie namespace

#endif // _EVENT_GENERATOR_LIST_H_
