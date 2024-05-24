//____________________________________________________________________________
/*!

\class    genie::NeutronOscPrimaryVtxGenerator

\brief    Utilities for simulating neutron oscillation

\author   Jeremy Hewes, Georgia Karagiorgi
          University of Manchester

          Adapted from the NucleonDecay package (Author: Costas Andreopoulos).

\created  November, 2016

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _NNBAR_OSC_PRIMARY_VTX_GENERATOR_H_
#define _NNBAR_OSC_PRIMARY_VTX_GENERATOR_H_

#include <array>

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

   static constexpr int fNPModes = 12;
   static constexpr int fNNModes = 32;
   static constexpr std::array< double, fNPModes + fNNModes > fBR = {
     0.001, 0.007, 0.148, 0.014, 0.020, 0.170, 0.108, 0.301,
     0.055, 0.032, 0.020, 0.124, //pnbar
     0.001, 0.007, 0.003, 0.010, 0.001, 0.003, 0.016, 0.131, 0.112, 0.033,
     0.014, 0.060, 0.136, 0.157, 0.006, 0.022, 0.020, 0.018, 0.037, 0.035,
     0.024, 0.027, 0.071, 0.016, 0.017, 0.001, 0.002, 0.001, 0.003, 0.003,
     0.010, 0.003
   };

public:

  inline int GetNPModes() const { return fNPModes; }
  inline int GetNNModes() const { return fNNModes; }
  inline int GetNModes() const { return fNPModes + fNNModes; }
  inline std::array< double, fNPModes + fNNModes > GetBRs() const
    { return fBR; }

};

} // genie namespace

#endif // _NNBAR_OSC_PRIMARY_VTX_GENERATOR_H_
