//____________________________________________________________________________
/*!

\class    genie::NuETargetRemnantGenerator

\brief    Generates all the non-primary lepton final state particles in 
          neutrino-electron events.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  July 17, 2005

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUE_TARGET_REMNANT_GENERATOR_H_
#define _NUE_TARGET_REMNANT_GENERATOR_H_

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class NuETargetRemnantGenerator : public EventRecordVisitorI {

public :
  NuETargetRemnantGenerator();
  NuETargetRemnantGenerator(string config);
 ~NuETargetRemnantGenerator();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * evrec) const;

private:
  void AddElectronNeutrino     (GHepRecord * evrec) const;
  void AddTargetNucleusRemnant (GHepRecord * evrec) const;
};

}      // genie namespace
#endif // _NUE_TARGET_REMNANT_GENERATOR_H_
