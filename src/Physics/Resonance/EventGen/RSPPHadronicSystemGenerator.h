//____________________________________________________________________________
/*!

\class    genie::RSPPHadronicSystemGenerator

\brief    Generates the 'final state' hadronic system in v SPP interactions.
          It adds the remnant nucleus (if any) and the baryon resonance decay 
          products at the GHEP record. The resonance decay products are pre-
          determined since in this thread we generate exclusive SPP reactions.
          The module uses a simple phase space decay.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 23, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RSPP_HADRONIC_SYSTEM_GENERATOR_H_
#define _RSPP_HADRONIC_SYSTEM_GENERATOR_H_

#include <TGenPhaseSpace.h>

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

class RSPPHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  RSPPHadronicSystemGenerator();
  RSPPHadronicSystemGenerator(string config);
 ~RSPPHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

private:
  void AddResonanceDecayProducts (GHepRecord * event_rec) const;

  mutable TGenPhaseSpace fPhaseSpaceGenerator;
};

}      // genie namespace
#endif // _RSPP_HADRONIC_SYSTEM_GENERATOR_H_
