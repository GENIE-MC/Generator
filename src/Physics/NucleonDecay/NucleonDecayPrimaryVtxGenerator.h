//____________________________________________________________________________
/*!

\class    genie::NucleonDecayPrimaryVtxGenerator

\brief    Utilities for simulating nucleon decay

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 03, 2011

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NUCLEON_DECAY_PRIMARY_VTX_GENERATOR_H_
#define _NUCLEON_DECAY_PRIMARY_VTX_GENERATOR_H_

#include <TGenPhaseSpace.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/NucleonDecay/NucleonDecayMode.h"

namespace genie {

class NuclearModelI;

class NucleonDecayPrimaryVtxGenerator: public EventRecordVisitorI {

public:
  NucleonDecayPrimaryVtxGenerator();
  NucleonDecayPrimaryVtxGenerator(string config);
 ~NucleonDecayPrimaryVtxGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

   void LoadConfig                     (void);
   void AddInitialState                (GHepRecord * event) const;
   void GenerateDecayedNucleonPosition (GHepRecord * event) const;
   void GenerateFermiMomentum          (GHepRecord * event) const;
   void GenerateDecayProducts          (GHepRecord * event) const;

   mutable int                fCurrInitStatePdg;
   mutable NucleonDecayMode_t fCurrDecayMode;
   mutable int                fCurrDecayedNucleon;
   mutable bool               fNucleonIsBound;
   mutable TGenPhaseSpace     fPhaseSpaceGenerator;

   const NuclearModelI * fNuclModel;
};

} // genie namespace

#endif // _NUCLEON_DECAY_PRIMARY_VTX_GENERATOR_H_
