//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
         Queen Mary University of London

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <RVersion.h>
#include <TClonesArray.h>
// Avoid the inclusion of dlfcn.h by Pythia.h that CINT is not able to process
#ifdef __CINT__
#define _DLFCN_H_
#define _DLFCN_H
#endif

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Hadronization/Pythia8Hadronization.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Pythia8/Pythia.h"
#endif

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
Pythia8Hadronization::Pythia8Hadronization() :
PythiaHadronizationBase("genie::Pythia8Hadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Hadronization::Pythia8Hadronization(string config) :
PythiaHadronizationBase("genie::Pythia8Hadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Hadronization::~Pythia8Hadronization()
{

}
//____________________________________________________________________________
void Pythia8Hadronization::ProcessEventRecord(GHepRecord *
#ifdef __GENIE_PYTHIA8_ENABLED__
  event // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  PythiaHadronizationBase::ProcessEventRecord(event);
#else
  LOG("Pythia8Had", pFATAL)
    << "Calling GENIE/PYTHIA8 hadronization modules without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif
}
//____________________________________________________________________________
TClonesArray * Pythia8Hadronization::Hadronize(const Interaction *
#ifdef __GENIE_PYTHIA8_ENABLED__
  interaction // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  LOG("Pythia8Had", pNOTICE) << "Running PYTHIA8 hadronizer";

  const Kinematics & kinematics = interaction->Kine();
  double W = kinematics.W();

  LOG("Pythia8Had", pNOTICE)
    << "Fragmentation: "
    << "q = " << fLeadingQuark << ", qq = " << fRemnantDiquark
    << ", W = " << W;

  // Hadronize

  double mA    = fPythia->particleData.m0(fLeadingQuark);
  double mB    = fPythia->particleData.m0(fRemnantDiquark);
  double pzAcm = 0.5 * Pythia8::sqrtpos( (W + mA + mB) * (W - mA - mB) * (W - mA + mB) * (W + mA - mB) ) / W;
  double pzBcm = -pzAcm;
  double eA    = sqrt(mA*mA + pzAcm*pzAcm);
  double eB    = sqrt(mB*mB + pzBcm*pzBcm);

  fPythia->event.reset();

  // Pythia8 status code for outgoing particles of the hardest subprocesses is 23
  // anti/colour tags for these 2 particles must complement each other
  fPythia->event.append(fLeadingQuark,   23, 101, 0, 0., 0., pzAcm, eA, mA);
  fPythia->event.append(fRemnantDiquark, 23, 0, 101, 0., 0., pzBcm, eB, mB);
  fPythia->next();

  // List the event information
  fPythia->event.list();
  fPythia->stat();

  // Get LUJETS record
  Pythia8::Event &fEvent = fPythia->event;
  int np = fEvent.size();
  assert(np>0);

  // Offset the initial (system) particle
  int ioff = 0;
  if (fEvent[0].id() == 90) ioff = -1;

  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", np);
  particle_list->SetOwner(true);

  // Hadronic 4vec
  TLorentzVector p4Had = kinematics.HadSystP4();

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = p4Had.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4Had.P()/p4Had.Energy());

  for (int i = 1; i < np; ++i) {
     /*
      * Convert Pythia8 status code to Pythia6
      * Initial quark has a pythia6 status code of 12
      * The initial diquark and the fragmented particles have a pythia6 code
      * of 11 (kIStNucleonTarget)
      * Final state particles have a positive pythia8 code and a pythia6 code of
      * 1 (kIStStableFinalState)
      * The fragmentation products are generated in the hadronic CM frame
      * where the z>0 axis is the \vec{phad} direction. For each particle
      * returned by the hadronizer:
      * - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
      * - rotate its 3-momentum from LAB' to LAB
      */
     GHepStatus_t gStatus;
     if (i == 1) gStatus = kIStDISPreFragmHadronicState;
     else gStatus = (fEvent[i].status()>0) ? kIStStableFinalState : kIStNucleonTarget;

     LOG("Pythia8Had", pDEBUG)
         << "Adding final state particle pdgc = " << fEvent[i].id()
         << " with status = " << gStatus;

     if (fEvent[i].status() > 0){
       if( pdg::IsQuark  (fEvent[i].id()) ||
               pdg::IsDiQuark(fEvent[i].id()) ) {
         LOG("Pythia8Had", pERROR)
             << "Hadronization failed! Bare quark/di-quarks appear in final state!";
         particle_list->Delete();
         delete particle_list;
         return 0;
       }
     }

     TLorentzVector p4o(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
     p4o.Boost(beta);
     TVector3 p3 = p4o.Vect();
     p3.RotateUz(unitvq);
     TLorentzVector p4(p3,p4o.Energy());

     // insert the particle in the list
     new((*particle_list)[i]) GHepParticle(
             fEvent[i].id(),               // pdg
             gStatus,                      // status
             fEvent[i].mother1()   + ioff, // first parent
             fEvent[i].mother2()   + ioff, // second parent
             fEvent[i].daughter1() + ioff, // first daughter
             fEvent[i].daughter2() + ioff, // second daughter
             p4.Px(),                      // px [GeV/c]
             p4.Py(),                      // py [GeV/c]
             p4.Pz(),                      // pz [GeV/c]
             p4.Energy(),                  // e  [GeV]
             fEvent[i].xProd(),            // x  [mm]
             fEvent[i].yProd(),            // y  [mm]
             fEvent[i].zProd(),            // z  [mm]
             fEvent[i].tProd()             // t  [mm/c]
     );
  }
  return particle_list;

#else
  return 0;
#endif
}
//____________________________________________________________________________
void Pythia8Hadronization::CopyOriginalDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fOriDecayFlag_pi0 = fPythia->particleData.canDecay(kPdgPi0);
  fOriDecayFlag_K0  = fPythia->particleData.canDecay(kPdgK0);
  fOriDecayFlag_K0b = fPythia->particleData.canDecay(kPdgAntiK0);
  fOriDecayFlag_L0  = fPythia->particleData.canDecay(kPdgLambda);
  fOriDecayFlag_L0b = fPythia->particleData.canDecay(kPdgAntiLambda);
  fOriDecayFlag_Dm  = fPythia->particleData.canDecay(kPdgP33m1232_DeltaM);
  fOriDecayFlag_D0  = fPythia->particleData.canDecay(kPdgP33m1232_Delta0);
  fOriDecayFlag_Dp  = fPythia->particleData.canDecay(kPdgP33m1232_DeltaP);
  fOriDecayFlag_Dpp = fPythia->particleData.canDecay(kPdgP33m1232_DeltaPP);

  //#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("Pythia8Had", pDEBUG)
       << "Original PYTHIA6 decay flags:"
       << "\n  pi0           =  " << fOriDecayFlag_pi0
       << "\n  K0            =  " << fOriDecayFlag_K0
       << "\n  \bar{K0}      =  " << fOriDecayFlag_K0b
       << "\n  Lambda        =  " << fOriDecayFlag_L0
       << "\n  \bar{Lambda0} =  " << fOriDecayFlag_L0b
       << "\n  D-            =  " << fOriDecayFlag_Dm
       << "\n  D0            =  " << fOriDecayFlag_D0
       << "\n  D+            =  " << fOriDecayFlag_Dp
       << "\n  D++           =  " << fOriDecayFlag_Dpp;
  //#endif

#endif
}
//____________________________________________________________________________
void Pythia8Hadronization::SetDesiredDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia->particleData.mayDecay(kPdgPi0,              fReqDecayFlag_pi0 );
  fPythia->particleData.mayDecay(kPdgK0,               fReqDecayFlag_K0  );
  fPythia->particleData.mayDecay(kPdgAntiK0,           fReqDecayFlag_K0b );
  fPythia->particleData.mayDecay(kPdgLambda,           fReqDecayFlag_L0  );
  fPythia->particleData.mayDecay(kPdgAntiLambda,       fReqDecayFlag_L0b );
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaM,  fReqDecayFlag_Dm  );
  fPythia->particleData.mayDecay(kPdgP33m1232_Delta0,  fReqDecayFlag_D0  );
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaP,  fReqDecayFlag_Dp  );
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaPP, fReqDecayFlag_Dpp );
#endif
}
//____________________________________________________________________________
void Pythia8Hadronization::RestoreOriginalDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia->particleData.mayDecay(kPdgPi0,              fOriDecayFlag_pi0 );
  fPythia->particleData.mayDecay(kPdgK0,               fOriDecayFlag_K0  );
  fPythia->particleData.mayDecay(kPdgAntiK0,           fOriDecayFlag_K0b );
  fPythia->particleData.mayDecay(kPdgLambda,           fOriDecayFlag_L0  );
  fPythia->particleData.mayDecay(kPdgAntiLambda,       fOriDecayFlag_L0b );
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaM,  fOriDecayFlag_Dm  );
  fPythia->particleData.mayDecay(kPdgP33m1232_Delta0,  fOriDecayFlag_D0  );
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaP,  fOriDecayFlag_Dp  );
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaPP, fOriDecayFlag_Dpp );
#endif
}
//____________________________________________________________________________
void Pythia8Hadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia8Hadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia8Hadronization::LoadConfig(void)
{
  PythiaHadronizationBase::LoadConfig();

#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia->settings.parm("StringFlav:probStoUD",         fSSBarSuppression);
  fPythia->settings.parm("Diffraction:primKTwidth",      fGaussianPt2);
  fPythia->settings.parm("StringPT:enhancedFraction",    fNonGaussianPt2Tail);
  fPythia->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);
  fPythia->settings.parm("StringFlav:probQQtoQ",         fDiQuarkSuppression);
  fPythia->settings.parm("StringFlav:mesonUDvector",     fLightVMesonSuppression);
  fPythia->settings.parm("StringFlav:mesonSvector",      fSVMesonSuppression);
  fPythia->settings.parm("StringZ:aLund",                fLunda);
  fPythia->settings.parm("StringZ:bLund",                fLundb);
  fPythia->settings.parm("StringZ:aExtraDiquark",        fLundaDiq);
#endif

  LOG("Pythia8Had", pDEBUG) << this->GetConfig();
}
//____________________________________________________________________________
void Pythia8Hadronization::Initialize(void)
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = new Pythia8::Pythia();

  // sync GENIE/PYTHIA8 seed number
  RandomGen::Instance(); // <---- NOT syncing yet - To do! Rephaps this need revisiting
#endif
}
//____________________________________________________________________________
