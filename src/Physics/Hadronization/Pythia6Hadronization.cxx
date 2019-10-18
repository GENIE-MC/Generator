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
#include "Physics/Hadronization/Pythia6Hadronization.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"

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
Pythia6Hadronization::Pythia6Hadronization() :
PythiaHadronizationBase("genie::Pythia6Hadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia6Hadronization::Pythia6Hadronization(string config) :
PythiaHadronizationBase("genie::Pythia6Hadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
Pythia6Hadronization::~Pythia6Hadronization()
{

}
//____________________________________________________________________________
void Pythia6Hadronization::ProcessEventRecord(GHepRecord *
  #ifdef __GENIE_PYTHIA6_ENABLED__
    event // avoid unused variable warning if PYTHIA6 is not enabled
  #endif
) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__
  PythiaHadronizationBase::ProcessEventRecord(event);
#else
  LOG("Pythia6Had", pFATAL)
    << "Calling GENIE/PYTHIA6 hadronization modules without enabling PYTHIA6";
  gAbortingInErr = true;
  std::exit(1);
#endif
}
//____________________________________________________________________________
TClonesArray * Pythia6Hadronization::Hadronize(const Interaction *
#ifdef __GENIE_PYTHIA6_ENABLED__
  interaction // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__

  LOG("Pythia6Had", pNOTICE) << "Running PYTHIA6 hadronizer";

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

  unsigned int i = 0;
  TMCParticle * p = 0;
  TIter particle_iter(pythia_particles);

  // Hadronic 4vec
  TLorentzVector p4Had = kinematics.HadSystP4();

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = p4Had.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4Had.P()/p4Had.Energy());

  while( (p = (TMCParticle *) particle_iter.Next()) ) {
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

     // Convert from TMCParticle to GHepParticle
     GHepParticle particle = GHepParticle(
         p->GetKF(),                // pdg
         GHepStatus_t(p->GetKS()),  // status
         p->GetParent(),            // first parent
         -1,                        // second parent
         p->GetFirstChild(),        // first daughter
         p->GetLastChild(),         // second daughter
         p4.Px(),                   // px
         p4.Py(),                   // py
         p4.Pz(),                   // pz
         p4.Energy(),               // e
         p->GetVx(),                // x
         p->GetVy(),                // y
         p->GetVz(),                // z
         p->GetTime()               // t
     );

     LOG("Pythia6Had", pDEBUG)
          << "Adding final state particle pdgc = " << particle.Pdg()
          << " with status = " << particle.Status();

     if(particle.Status() == 1) {
       if( pdg::IsQuark  (particle.Pdg()) ||
	         pdg::IsDiQuark(particle.Pdg()) ) {

          LOG("Pythia6Had", pERROR)
	           << "Hadronization failed! Bare quark/di-quarks appear in final state!";

	        particle_list->Delete();
	        delete particle_list;
	        return 0;
       }
     }

     // insert the particle in the list
     new ( (*particle_list)[i++] ) GHepParticle(particle);
  }
  return particle_list;

#else
  return 0;
#endif // __GENIE_PYTHIA6_ENABLED__
}
//____________________________________________________________________________
void Pythia6Hadronization::CopyOriginalDecayFlags(void) const
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

#endif // // __GENIE_PYTHIA6_ENABLED__
}
//____________________________________________________________________________
void Pythia6Hadronization::SetDesiredDecayFlags(void) const
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
void Pythia6Hadronization::RestoreOriginalDecayFlags(void) const
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
void Pythia6Hadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia6Hadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void Pythia6Hadronization::LoadConfig(void)
{
  PythiaHadronizationBase::LoadConfig();

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
void Pythia6Hadronization::Initialize(void)
{
  PythiaHadronizationBase::Initialize();
#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();
  RandomGen::Instance(); // sync GENIE/PYTHIA6 seed number
#endif
}
//____________________________________________________________________________
