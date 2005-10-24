//____________________________________________________________________________
/*!

\class   genie::HadronicSystemGenerator

\brief   Abstract class. Is used to pass some commonly recurring methods to
         all concrete implementations of the EventRecordVisitorI interface
         generating the hadronic system for a specific processes (QEL,DIS,
         RES,...)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created Juy 16, 2005

*/
//____________________________________________________________________________

#ifndef _HADRONIC_SYSTEM_GENERATOR_H_
#define _HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class HadronicSystemGenerator : public EventRecordVisitorI {

public :

  //-- Do not implement the EventRecordVisitorI interface.
  //   Leave it for its children that are EventRecordVisitors too

  //-- Common methods for all concrete PrimaryLeptonGenerator-type
  //   EventRecordVisitors

  void AddTargetNucleusRemnant (GHepRecord * event_rec) const;

protected:

  //-- Abstract class - Can only be instantiated by its children.
  HadronicSystemGenerator();
  HadronicSystemGenerator(string name);
  HadronicSystemGenerator(string name, string config);
  ~HadronicSystemGenerator();
};

}      // genie namespace

#endif // _HADRONIC_SYSTEM_GENERATOR_H_
