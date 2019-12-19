//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <RVersion.h>
#include <TClonesArray.h>

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
#include "Physics/Hadronization/Pythia6Hadro2019.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#endif // __GENIE_PYTHIA6_ENABLED__

using namespace genie;
using namespace genie::constants;

// the actual PYTHIA call
extern "C" void py2ent_(int *,  int *, int *, double *);

//____________________________________________________________________________
Pythia6Hadro2019::Pythia6Hadro2019() :
PythiaBaseHadro2019("genie::Pythia6Hadro2019")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia6Hadro2019::Pythia6Hadro2019(string config) :
PythiaBaseHadro2019("genie::Pythia6Hadro2019", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia6Hadro2019::~Pythia6Hadro2019()
{

}
//____________________________________________________________________________
void Pythia6Hadro2019::ProcessEventRecord(GHepRecord *
  #ifdef __GENIE_PYTHIA6_ENABLED__
    event // avoid unused variable warning if PYTHIA6 is not enabled
  #endif
) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__
  PythiaBaseHadro2019::ProcessEventRecord(event);
#else
  LOG("Pythia6Had", pFATAL)
    << "Calling GENIE/PYTHIA6 hadronization modules without enabling PYTHIA6";
  gAbortingInErr = true;
  std::exit(1);
#endif
}
//____________________________________________________________________________
bool Pythia6Hadro2019::Hadronize(GHepRecord *
#ifdef __GENIE_PYTHIA6_ENABLED__
  event // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__

  LOG("Pythia6Had", pNOTICE) << "Running PYTHIA6 hadronizer";

  const Interaction * interaction = event->Summary();
  const Kinematics & kinematics = interaction->Kine();
  double W = kinematics.W();

  LOG("Pythia6Had", pNOTICE)
    << "Fragmentation: "
    << "q = " << fLeadingQuark << ", qq = " << fRemnantDiquark
    << ", W = " << W;

  // Hadronize
  int ip = 0;
  py2ent_(&ip, &fLeadingQuark, &fRemnantDiquark, &W); // hadronizer

  // Get LUJETS record
  fPythia->GetPrimaries();
  TClonesArray * pythia_particles =
       (TClonesArray *) fPythia->ImportParticles("All");

  // Copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method

  int np = pythia_particles->GetEntries();
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
  unsigned int i = 0;
  TMCParticle * p = 0;
  TIter particle_iter(pythia_particles);
  while( (p = (TMCParticle *) particle_iter.Next()) ) {

     int particle_pdg_code      = p->GetKF();
     int pythia_particle_status = p->GetKS();

     // Sanity check
     if(pythia_particle_status == 1) {
       if( pdg::IsQuark  (particle_pdg_code) ||
           pdg::IsDiQuark(particle_pdg_code) )
       {
         LOG("Pythia6Had", pERROR)
            << "Hadronization failed! Bare quarks appear in final state!";
         return false;
       }
     }

     // The fragmentation products are generated in the hadronic CM frame
     // where the z>0 axis is the \vec{phad} direction. For each particle
     // returned by the hadronizer:
     // - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
     // - rotate its 3-momentum from LAB' to LAB
     TLorentzVector p4o(p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy());
     p4o.Boost(beta);
     TVector3 p3 = p4o.Vect();
     p3.RotateUz(unitvq);
     TLorentzVector p4(p3,p4o.Energy());

     // Set the proper GENIE status according to a number of things:
     // interaction on a nucleus or nucleon, particle type
     GHepStatus_t ist = (pythia_particle_status == 1) ?
         istfin : kIStDISPreFragmHadronicState;
     // Handle gammas, and leptons that might come from internal pythia decays
     // mark them as final state particles
     bool is_gamma = (particle_pdg_code == kPdgGamma);
     bool is_nu    = pdg::IsNeutralLepton(particle_pdg_code);
     bool is_lchg  = pdg::IsChargedLepton(particle_pdg_code);
     bool not_hadr = is_gamma || is_nu || is_lchg;
     if(not_hadr)  { ist = kIStStableFinalState; }

     // Set mother/daugher indices
     int mother1   = mom + p->GetParent();
     int mother2   = -1;
     int daughter1 = (p->GetFirstChild() <= 0 ) ? -1 : mom  + p->GetFirstChild();
     int daughter2 = (p->GetLastChild()  <= 0 ) ? -1 : mom  + p->GetLastChild();

     // Create GHepParticle
     GHepParticle particle = GHepParticle(
         particle_pdg_code,  // pdg
         ist,                // status
         mother1,            // first parent
         mother2,            // second parent
         daughter1,          // first daughter
         daughter2,          // second daughter
         p4.Px(),            // px
         p4.Py(),            // py
         p4.Pz(),            // pz
         p4.Energy(),        // e
         vtx.X(),            // x
         vtx.Y(),            // y
         vtx.Z(),            // z
         vtx.T()             // t
     );

     LOG("Pythia6Had", pDEBUG)
          << "Adding final state particle pdgc = " << particle.Pdg()
          << " with status = " << particle.Status();

     // Insert the particle in the list
     event->AddParticle(particle);
  }
  return true;

#else
  return false;
#endif // __GENIE_PYTHIA6_ENABLED__
}
//____________________________________________________________________________
void Pythia6Hadro2019::CopyOriginalDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__
  fOriDecayFlag_pi0 = (fPythia->GetMDCY(fPythia->Pycomp(kPdgPi0),              1) == 1);
  fOriDecayFlag_K0  = (fPythia->GetMDCY(fPythia->Pycomp(kPdgK0),               1) == 1);
  fOriDecayFlag_K0b = (fPythia->GetMDCY(fPythia->Pycomp(kPdgAntiK0),           1) == 1);
  fOriDecayFlag_L0  = (fPythia->GetMDCY(fPythia->Pycomp(kPdgLambda),           1) == 1);
  fOriDecayFlag_L0b = (fPythia->GetMDCY(fPythia->Pycomp(kPdgAntiLambda),       1) == 1);
  fOriDecayFlag_Dm  = (fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),  1) == 1);
  fOriDecayFlag_D0  = (fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),  1) == 1);
  fOriDecayFlag_Dp  = (fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),  1) == 1);
  fOriDecayFlag_Dpp = (fPythia->GetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP), 1) == 1);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("Pythia6Had", pDEBUG)
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

#endif // __GENIE_PYTHIA6_ENABLED__
}
//____________________________________________________________________________
void Pythia6Hadro2019::SetDesiredDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__

  fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),
      1, (fReqDecayFlag_pi0) ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgK0),
     1, (fReqDecayFlag_K0)  ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiK0),
     1, (fReqDecayFlag_K0b) ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgLambda),
     1, (fReqDecayFlag_L0)  ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiLambda),
     1, (fReqDecayFlag_L0b) ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),
     1, (fReqDecayFlag_Dm)  ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),
     1, (fReqDecayFlag_D0)  ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),
     1, (fReqDecayFlag_Dp)  ? 1 : 0);
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP),
     1, (fReqDecayFlag_Dpp) ? 1 : 0);

#endif // __GENIE_PYTHIA6_ENABLED__
}
//____________________________________________________________________________
void Pythia6Hadro2019::RestoreOriginalDecayFlags(void) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__

fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),
   1, (fOriDecayFlag_pi0) ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgK0),
   1, (fOriDecayFlag_K0)  ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiK0),
   1, (fOriDecayFlag_K0b) ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgLambda),
   1, (fOriDecayFlag_L0)  ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgAntiLambda),
   1, (fOriDecayFlag_L0b) ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),
   1, (fOriDecayFlag_Dm)  ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),
   1, (fOriDecayFlag_D0)  ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),
   1, (fOriDecayFlag_Dp)  ? 1 : 0);
fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP),
   1, (fOriDecayFlag_Dpp) ? 1 : 0);

#endif  // __GENIE_PYTHIA6_ENABLED__
}
//____________________________________________________________________________
void Pythia6Hadro2019::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia6Hadro2019::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia6Hadro2019::LoadConfig(void)
{
  PythiaBaseHadro2019::LoadConfig();

#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia->SetPARJ(2,  fSSBarSuppression       );
  fPythia->SetPARJ(21, fGaussianPt2            );
  fPythia->SetPARJ(23, fNonGaussianPt2Tail     );
  fPythia->SetPARJ(33, fRemainingECutoff       );
  fPythia->SetPARJ(1,  fDiQuarkSuppression     );
  fPythia->SetPARJ(11, fLightVMesonSuppression );
  fPythia->SetPARJ(12, fSVMesonSuppression     );
  fPythia->SetPARJ(41, fLunda                  );
  fPythia->SetPARJ(42, fLundb                  );
  fPythia->SetPARJ(45, fLundaDiq               );
#endif

  LOG("Pythia6Had", pDEBUG) << this->GetConfig() ;
}
//____________________________________________________________________________
void Pythia6Hadro2019::Initialize(void)
{
  PythiaBaseHadro2019::Initialize();
#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();
  // sync GENIE/PYTHIA6 seed number
  // PYTHIA6 is a singleton, so do this from RandomGen for all
  // GENIE algorithms that use PYTHIA6
  RandomGen::Instance();
#endif
}
//____________________________________________________________________________
