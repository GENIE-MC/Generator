//____________________________________________________________________________
/*!

\class    genie::NeutronOscPrimaryVtxGenerator

\brief    Utilities for simulating neutron oscillation

\author   Jeremy Hewes, Georgia Karagiorgi
          University of Manchester

          Adapted from the NucleonDecay package (Author: Costas Andreopoulos).

\created  November, 2016

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NNBAR_OSC_PRIMARY_VTX_GENERATOR_H_
#define _NNBAR_OSC_PRIMARY_VTX_GENERATOR_H_

#include <TGenPhaseSpace.h>
#include <TFile.h>
#include <TH1.h>
#include <string>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/NNBarOscillation/NNBarOscMode.h"

namespace genie {

class NuclearModelI;

class NNBarOscPrimaryVtxGenerator: public EventRecordVisitorI {

public:
  NNBarOscPrimaryVtxGenerator();
  NNBarOscPrimaryVtxGenerator(string config);
 ~NNBarOscPrimaryVtxGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

   void LoadConfig                         (void);
   void AddInitialState                    (GHepRecord * event) const;
   void GenerateOscillatingNeutronPosition (GHepRecord * event) const;
   void GenerateFermiMomentum              (GHepRecord * event) const;
   void GenerateDecayProducts              (GHepRecord * event) const;

   mutable int                fCurrInitStatePdg;
   mutable NNBarOscMode_t     fCurrDecayMode;
   mutable bool               fNucleonIsBound;
   mutable TGenPhaseSpace     fPhaseSpaceGenerator;

   const NuclearModelI * fNuclModel;
};

} // genie namespace

#endif // _NNBAR_OSC_PRIMARY_VTX_GENERATOR_H_
