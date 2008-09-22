//____________________________________________________________________________
/*!

\class    genie::MECGenerator

\brief    

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep. 22, 2008

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_GENERATOR_H_
#define _MEC_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class MECGenerator : public EventRecordVisitorI {

public :
  MECGenerator();
  MECGenerator(string config);
 ~MECGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

private:

};

}      // genie namespace
#endif // _MEC_GENERATOR_H_
