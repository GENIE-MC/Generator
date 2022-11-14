//____________________________________________________________________________
/*!

\class    genie::HNL::HNLDecayer

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

#include <cassert>

#include <TGenPhaseSpace.h>
#include <TH3.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "Tools/Flux/GNuMIFlux.h"

#include "Physics/BeamHNL/HNLDecayMode.h"
#include "Physics/BeamHNL/HNLFluxReader.h"
#include "Physics/BeamHNL/SimpleHNL.h"

namespace genie {

namespace HNL {

class HNLDecayer: public EventRecordVisitorI {

public:
  HNLDecayer();
  HNLDecayer(string config);
 ~HNLDecayer();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

  // required public to reroll vertex positioning in main app 
  std::vector< double > * GenerateDecayPosition (GHepRecord * event) const;
  std::vector< double > * GenerateMomentum (GHepRecord * event) const;

  genie::HNL::SimpleHNL GetHNLInstance(string config) const;

  // get information about parent and polarisation from HNLFluxCreator
  void ReadCreationInfo( genie::flux::GNuMIFluxPassThroughInfo gnmf ) const;

private:

   void LoadConfig            (void);
   void AddInitialState       (GHepRecord * event) const;
   void GenerateDecayProducts (GHepRecord * event) const;
   void UpdateEventRecord     (GHepRecord * event) const;
   void SetHNLCouplings       (double Ue42, double Um42, double Ut42) const;
   void SetBeam2User          (std::vector<double> translation, std::vector<double> rotation) const; 
   void SetProdVtxPosition    (const TLorentzVector & v4) const; // in detector coordinates

   // PolMag is a bit of a misnomer. It is a polarisation modulus, i.e. not positive semidefinite.
   double CalcPolMag          (int parPdg, int lepPdg, double M) const;
   double CalcPolMod          (double polMag, int lepPdg, int hadPdg, double M) const;

   void UnpolarisedDecay      (TGenPhaseSpace & fPSG, PDGCodeList pdgv, double wm, bool failed) const;
   void PolarisedDecay        (TGenPhaseSpace & fPSG, PDGCodeList pdgv, double wm, TVector3 vPolDir, bool failed)const;

   mutable int                        fCurrInitStatePdg;
   mutable genie::HNL::HNLDecayMode_t fCurrDecayMode;

   mutable int                        fParentPdg, fProdLepPdg, fDecLepPdg, fDecHadPdg;
   mutable std::vector<double>        fPolDir;

   mutable bool                       fIsConfigLoaded = false;
   
   mutable double                     fMass; 
   mutable double                     fUe42 = -1.0, fUm42 = -1.0, fUt42 = -1.0;
   mutable bool                       fIsMajorana = false;
   mutable int                        fType = 2;

   mutable bool                       fDoPol = false;

   mutable std::vector< genie::HNL::HNLDecayMode_t > fIntChannels;

   mutable double                     fAngularDeviation = -1.0;
   mutable std::vector< double >      fB2UTranslation;
   mutable double                     fTx = -1.0, fTy = -1.0, fTz = -1.0;
   mutable std::vector< double >      fB2URotation;
   mutable double                     fR1 = -1.0, fR2 = -1.0, fR3 = -1.0;

   mutable TH3D *                     fProdVtxHist = 0;
   mutable TLorentzVector *           fProdVtx = 0;
   mutable TLorentzVector *           fISMom = 0;

   mutable bool                       fGetCMFrameInstead = false;
};

} // HNL namespace

} // genie namespace

#endif // _NEUTRAL_HEAVY_LEPTON_PRIMARY_VTX_GENERATOR_H_
