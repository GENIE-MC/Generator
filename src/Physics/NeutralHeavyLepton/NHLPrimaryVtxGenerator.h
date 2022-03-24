//____________________________________________________________________________
/*!

\class    genie::NHLPrimaryVtxGenerator

\brief    Neutral Heavy Lepton primary vertex generator

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory
	  John Plows <komninos-john.plows \at physics.ox.ac.uk>

\created  February 10, 2020

\cpright  Copyright (c) 2003-2022, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org        
*/
//____________________________________________________________________________

#ifndef _NEUTRAL_HEAVY_LEPTON_PRIMARY_VTX_GENERATOR_H_
#define _NEUTRAL_HEAVY_LEPTON_PRIMARY_VTX_GENERATOR_H_

#include <TGenPhaseSpace.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Numerical/RandomGen.h"

#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"
#include "Physics/NeutralHeavyLepton/NHLEnums.h" // to be removed later
#include "Physics/NeutralHeavyLepton/SimpleNHL.h"

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
   // to be introduced later!
   //void GenerateDecayPosition (GHepRecord * event) const;
   void SetNHLCouplings          (double Ue42, double Um42, double Ut42) const;

   mutable int                        fCurrInitStatePdg;
   mutable genie::NHL::NHLDecayMode_t fCurrDecayMode;
   mutable TGenPhaseSpace             fPhaseSpaceGenerator;

   mutable double                     fEnergy;
   mutable double                     fUe42 = -1.0, fUm42 = -1.0, fUt42 = -1.0;
   // to be implemented!
   //mutable TH3D *                     fProdVtxHist = 0;
};

} // genie namespace

#endif // _NEUTRAL_HEAVY_LEPTON_PRIMARY_VTX_GENERATOR_H_
