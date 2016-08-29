//____________________________________________________________________________
/*!

\class    genie::NeutronOscPrimaryVtxGenerator

\brief    Utilities for simulating neutron oscillation

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 03, 2011

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NEUTRON_OSC_PRIMARY_VTX_GENERATOR_H_
#define _NEUTRON_OSC_PRIMARY_VTX_GENERATOR_H_

#include <TGenPhaseSpace.h>
#include <TFile.h>
#include <TH1.h>
#include <string>

#include "EVGCore/EventRecordVisitorI.h"
#include "NeutronOsc/NeutronOscMode.h"

namespace genie {

class NuclearModelI;

class NeutronOscPrimaryVtxGenerator: public EventRecordVisitorI {

public:
  NeutronOscPrimaryVtxGenerator();
  NeutronOscPrimaryVtxGenerator(string config);
 ~NeutronOscPrimaryVtxGenerator();

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
   mutable NeutronOscMode_t   fCurrDecayMode;
   mutable bool               fNucleonIsBound;
   mutable TGenPhaseSpace     fPhaseSpaceGenerator;

   const NuclearModelI * fNuclModel;
};

} // genie namespace

#endif // _NEUTRON_OSC_PRIMARY_VTX_GENERATOR_H_
