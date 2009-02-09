//____________________________________________________________________________
/*!

\class    genie::DfrcGenerator

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Feb 15, 2008

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIFFRACTIVE_GENERATOR_H_
#define _DIFFRACTIVE_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class DfrcGenerator : public EventRecordVisitorI {

public :
  DfrcGenerator();
  DfrcGenerator(string config);
 ~DfrcGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

private:

};

}      // genie namespace
#endif // _DIFFRACTIVE_GENERATOR_H_
