//____________________________________________________________________________
/*!

\class    genie::NHLPrimaryVtxGenerator

\brief    Neutral Heavy Lepton primary vertex generator

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  February 10, 2020

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _NEUTRAL_HEAVY_LEPTON_PRIMARY_VTX_GENERATOR_H_
#define _NEUTRAL_HEAVY_LEPTON_PRIMARY_VTX_GENERATOR_H_

#include <TGenPhaseSpace.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"

namespace genie {

class NHLPrimaryVtxGenerator: public EventRecordVisitorI {

public:
  NHLPrimaryVtxGenerator();
  NHLPrimaryVtxGenerator(string config);
 ~NHLPrimaryVtxGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

   void LoadConfig            (void);
   void AddInitialState       (GHepRecord * event) const;
   void GenerateDecayProducts (GHepRecord * event) const;

   mutable NHLDecayMode_t fCurrDecayMode;
   mutable TGenPhaseSpace fPhaseSpaceGenerator;
};

} // genie namespace

#endif // _NEUTRAL_HEAVY_LEPTON_PRIMARY_VTX_GENERATOR_H_
