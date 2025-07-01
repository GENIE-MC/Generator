//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

 Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
 Queen Mary University of London
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
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Hadronization/Pythia8Hadro2019.h"

#ifdef __GENIE_PYTHIA8_ENABLED__
#include "Pythia8/Pythia.h"
#endif

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
Pythia8Hadro2019::Pythia8Hadro2019() :
PythiaBaseHadro2019("genie::Pythia8Hadro2019")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Hadro2019::Pythia8Hadro2019(string config) :
PythiaBaseHadro2019("genie::Pythia8Hadro2019", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia8Hadro2019::~Pythia8Hadro2019()
{

}
//____________________________________________________________________________
void Pythia8Hadro2019::ProcessEventRecord(GHepRecord *
#ifdef __GENIE_PYTHIA8_ENABLED__
  event // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  PythiaBaseHadro2019::ProcessEventRecord(event);
#else
  LOG("Pythia8Had", pFATAL)
    << "Calling GENIE/PYTHIA8 hadronization modules without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif
}
//____________________________________________________________________________
bool Pythia8Hadro2019::Hadronize(GHepRecord *
#ifdef __GENIE_PYTHIA8_ENABLED__
  event // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  LOG("Pythia8Had", pNOTICE) << "Running PYTHIA8 hadronizer";

  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  const Interaction * interaction = event->Summary();
  const Kinematics & kinematics = interaction->Kine();
  double W = kinematics.W();

  LOG("Pythia8Had", pNOTICE)
    << "Fragmentation: "
    << "q = " << fLeadingQuark << ", qq = " << fRemnantDiquark
    << ", W = " << W << " GeV";

  // Hadronize

  LOG("Pythia8Had", pDEBUG) << "Reseting PYTHIA8 event";
  gPythia->event.reset();

  // Get quark/diquark masses
  double mA = gPythia->particleData.m0(fLeadingQuark);
  double mB = gPythia->particleData.m0(fRemnantDiquark);

  LOG("Pythia8Had", pINFO)
    << "Leading quark mass = " << mA
    << " GeV, remnant diqurak mass = " << mB << ", GeV";

  // Calculate quark/diquark energy/momentum
  double pzAcm = 0.5 * Pythia8::sqrtpos(
     (W + mA + mB) * (W - mA - mB) * (W - mA + mB) * (W + mA - mB) ) / W;
  double pzBcm = -pzAcm;
  double eA    = sqrt(mA*mA + pzAcm*pzAcm);
  double eB    = sqrt(mB*mB + pzBcm*pzBcm);

  LOG("Pythia8Had", pINFO)
   << "Quark: (pz = " << pzAcm << ", E = " << eA << ") GeV, "
   << "Diquark: (pz = " << pzBcm << ", E = " << eB << ") GeV";

  // Pythia8 status code for outgoing particles of the hardest subprocesses is 23
  // anti/colour tags for these 2 particles must complement each other
  LOG("Pythia8Had", pDEBUG) << "Appending quark/diquark into the PYTHIA8 event";
  gPythia->event.append(fLeadingQuark,   23, 101, 0, 0., 0., pzAcm, eA, mA);
  gPythia->event.append(fRemnantDiquark, 23, 0, 101, 0., 0., pzBcm, eB, mB);
  gPythia->event.list();

  LOG("Pythia8Had", pDEBUG) << "Generating next PYTHIA8 event";
  gPythia->next();

  // List the event information
  //gPythia->event.list();
  //gPythia->stat();

  // Get LUJETS record
  LOG("Pythia8Had", pDEBUG) << "Copying PYTHIA8 event record into GENIE's";
  Pythia8::Event &fEvent = gPythia->event;
  int np = fEvent.size();
  assert(np>0);

  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", np);
  particle_list->SetOwner(true);

  // Hadronic 4vec
  TLorentzVector p4Had = kinematics.HadSystP4();

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = p4Had.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4Had.P()/p4Had.Energy());

  // Check target and decide appropriate status code for f/s particles
  bool is_nucleus = interaction->InitState().Tgt().IsNucleus();
  GHepStatus_t istfin = is_nucleus ? kIStHadronInTheNucleus : kIStStableFinalState ;

  // Get the index of the mother of the hadronic system
  int mom = event->FinalStateHadronicSystemPosition();
  assert(mom!=-1);

  // Get the neutrino vertex position (all hadrons positions set to this point)
  GHepParticle * neutrino  = event->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  // Loop over PYTHIA8 event particles and copy relevant entries
  for (int i = 0; i < np; i++) {

     if (fEvent[i].id() == 90) continue; // ignore (system) pseudoparticle

     // Get PYTHIA8 particle ID and status codes
     int particle_pdg_code      = fEvent[i].id();
     int pythia_particle_status = fEvent[i].status();

     // Sanity check
     if(pythia_particle_status > 0) {
       if( pdg::IsQuark  (particle_pdg_code) ||
           pdg::IsDiQuark(particle_pdg_code) )
       {
          LOG("Pythia8Had", pERROR)
             << "Hadronization failed! Bare quarks appear in final state!";
          return false;
       }
     }

     // Copy the initial q, qq (status = -23) and all undecayed particles
     // Ignore particles we asked PYTHIA to decay (eg omega, Delta) but
     // record their decay products
     bool copy = (pythia_particle_status==-23) || (pythia_particle_status > 0);
     if(copy) {
        // The fragmentation products are generated in the hadronic CM frame
        // where the z>0 axis is the \vec{phad} direction. For each particle
        // returned by the hadronizer:
        // - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
        // - rotate its 3-momentum from LAB' to LAB
        TLorentzVector p4o(
          fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
        p4o.Boost(beta);
        TVector3 p3 = p4o.Vect();
        p3.RotateUz(unitvq);
        TLorentzVector p4(p3,p4o.Energy());

        // Set the proper GENIE status according to a number of things:
        // interaction on a nucleus or nucleon, particle type
        GHepStatus_t ist = (pythia_particle_status > 0) ?
           istfin : kIStDISPreFragmHadronicState;
        // Handle gammas, and leptons that might come from internal pythia decays
        // mark them as final state particles
        bool is_gamma = (particle_pdg_code == kPdgGamma);
        bool is_nu    = pdg::IsNeutralLepton(particle_pdg_code);
        bool is_lchg  = pdg::IsChargedLepton(particle_pdg_code);
        bool not_hadr = is_gamma || is_nu || is_lchg;
        if(not_hadr)  { ist = kIStStableFinalState; }

        // Set mother/daugher indices
        // int mother1 = mom+ fEvent[i].mother1();
        // int mother2 = (pythia_particle_status > 0) ? mom + fEvent[i].mother2() : -1;
        int mother1 = mom; // fEvent[i].mother1();
        int mother2 = -1; //(pythia_particle_status > 0) ? mom + fEvent[i].mother2() : -1;
        if(pythia_particle_status > 0) {
          mother1 = mom+1;
          mother2 = mom+2;
        }
        int daughter1 = -1;//(fEvent[i].daughter1() <= 0 ) ? -1 : mom  + fEvent[i].daughter1();
        int daughter2 = -1;//(fEvent[i].daughter1() <= 0 ) ? -1 : mom  + fEvent[i].daughter2();

        // Create GHepParticle
        GHepParticle particle = GHepParticle(
            particle_pdg_code, // pdg
            ist,               // status
            mother1,           // first parent
            mother2,           // second parent
            daughter1,         // first daughter
            daughter2,         // second daughter
            p4.Px(),           // px
            p4.Py(),           // py
            p4.Pz(),           // pz
            p4.Energy(),       // e
            vtx.X(),           // x
            vtx.Y(),           // y
            vtx.Z(),           // z
            vtx.T()            // t
        );

        LOG("Pythia8Had", pDEBUG)
             << "Adding final state particle pdgc = " << particle.Pdg()
             << " with status = " << particle.Status();

        // Insert the particle in the list
        event->AddParticle(particle);
     }// copy?
  }// loop over particles

  return true;

#else
  return false;
#endif
}
//____________________________________________________________________________
void Pythia8Hadro2019::CopyOriginalDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  fOriDecayFlag_pi0 = gPythia->particleData.canDecay(kPdgPi0);
  fOriDecayFlag_K0  = gPythia->particleData.canDecay(kPdgK0);
  fOriDecayFlag_K0b = gPythia->particleData.canDecay(kPdgAntiK0);
  fOriDecayFlag_L0  = gPythia->particleData.canDecay(kPdgLambda);
  fOriDecayFlag_L0b = gPythia->particleData.canDecay(kPdgAntiLambda);
  fOriDecayFlag_Dm  = gPythia->particleData.canDecay(kPdgP33m1232_DeltaM);
  fOriDecayFlag_D0  = gPythia->particleData.canDecay(kPdgP33m1232_Delta0);
  fOriDecayFlag_Dp  = gPythia->particleData.canDecay(kPdgP33m1232_DeltaP);
  fOriDecayFlag_Dpp = gPythia->particleData.canDecay(kPdgP33m1232_DeltaPP);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
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
#endif

#endif
}
//____________________________________________________________________________
void Pythia8Hadro2019::SetDesiredDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  gPythia->particleData.mayDecay(kPdgPi0,              fReqDecayFlag_pi0 );
  gPythia->particleData.mayDecay(kPdgK0,               fReqDecayFlag_K0  );
  gPythia->particleData.mayDecay(kPdgAntiK0,           fReqDecayFlag_K0b );
  gPythia->particleData.mayDecay(kPdgLambda,           fReqDecayFlag_L0  );
  gPythia->particleData.mayDecay(kPdgAntiLambda,       fReqDecayFlag_L0b );
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaM,  fReqDecayFlag_Dm  );
  gPythia->particleData.mayDecay(kPdgP33m1232_Delta0,  fReqDecayFlag_D0  );
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaP,  fReqDecayFlag_Dp  );
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaPP, fReqDecayFlag_Dpp );
#endif
}
//____________________________________________________________________________
void Pythia8Hadro2019::RestoreOriginalDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  gPythia->particleData.mayDecay(kPdgPi0,              fOriDecayFlag_pi0 );
  gPythia->particleData.mayDecay(kPdgK0,               fOriDecayFlag_K0  );
  gPythia->particleData.mayDecay(kPdgAntiK0,           fOriDecayFlag_K0b );
  gPythia->particleData.mayDecay(kPdgLambda,           fOriDecayFlag_L0  );
  gPythia->particleData.mayDecay(kPdgAntiLambda,       fOriDecayFlag_L0b );
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaM,  fOriDecayFlag_Dm  );
  gPythia->particleData.mayDecay(kPdgP33m1232_Delta0,  fOriDecayFlag_D0  );
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaP,  fOriDecayFlag_Dp  );
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaPP, fOriDecayFlag_Dpp );
#endif
}
//____________________________________________________________________________
void Pythia8Hadro2019::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia8Hadro2019::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia8Hadro2019::LoadConfig(void)
{
  PythiaBaseHadro2019::LoadConfig();

#ifdef __GENIE_PYTHIA8_ENABLED__
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  gPythia->settings.parm("StringFlav:probStoUD",         fSSBarSuppression);
  gPythia->settings.parm("Diffraction:primKTwidth",      fGaussianPt2);
  gPythia->settings.parm("StringPT:enhancedFraction",    fNonGaussianPt2Tail);
  gPythia->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);
  gPythia->settings.parm("StringFlav:probQQtoQ",         fDiQuarkSuppression);
  gPythia->settings.parm("StringFlav:mesonUDvector",     fLightVMesonSuppression);
  gPythia->settings.parm("StringFlav:mesonSvector",      fSVMesonSuppression);
  gPythia->settings.parm("StringZ:aLund",                fLunda);
  gPythia->settings.parm("StringZ:bLund",                fLundb);
  gPythia->settings.parm("StringZ:aExtraDiquark",        fLundaDiq);

  gPythia->init(); // needed again to read the above?

#endif

  LOG("Pythia8Had", pDEBUG) << this->GetConfig();
}
//____________________________________________________________________________
void Pythia8Hadro2019::Initialize(void)
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  gPythia->readString("ProcessLevel:all = off");
  gPythia->readString("Print:quiet      = on");

  // sync GENIE and PYTHIA8 seeds
  RandomGen * rnd = RandomGen::Instance();
  long int seed = rnd->GetSeed();
  gPythia->readString("Random:setSeed = on");
  gPythia->settings.mode("Random:seed", seed);
  LOG("Pythia8Had", pINFO)
    << "PYTHIA8  seed = " << gPythia->settings.mode("Random:seed");

  gPythia->init();

#endif
}
//____________________________________________________________________________
