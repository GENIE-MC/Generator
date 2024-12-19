//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

 Changes required to implement the GENIE Boosted Dark Matter module
 were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <RVersion.h>
#include <TVector3.h>
#include <TF1.h>
#include <TROOT.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/Hadronization/AGCharmPythia6Hadro2023.h"
#include "Physics/Hadronization/FragmentationFunctionI.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
  #if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
    #include <TMCParticle.h>
  #else
    #include <TMCParticle6.h>
  #endif

  #include <TPythia6.h>
#endif

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

#ifdef __GENIE_PYTHIA6_ENABLED__
extern "C" void py2ent_(int *,  int *, int *, double *);
#endif

//____________________________________________________________________________
AGCharmPythia6Hadro2023::AGCharmPythia6Hadro2023() :
AGCharmPythiaBaseHadro2023("genie::AGCharmPythia6Hadro2023")
{
  this->Initialize();
}
//____________________________________________________________________________
AGCharmPythia6Hadro2023::AGCharmPythia6Hadro2023(string config) :
AGCharmPythiaBaseHadro2023("genie::AGCharmPythia6Hadro2023", config)
{
  this->Initialize();
}
//____________________________________________________________________________
AGCharmPythia6Hadro2023::~AGCharmPythia6Hadro2023()
{

}
//____________________________________________________________________________
void AGCharmPythia6Hadro2023::Initialize(void) const
{
  AGCharmPythiaBaseHadro2023::Initialize();
#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();
  // sync GENIE/PyTHIA6 seed number
  // PYTHIA6 is a singleton, so do this from RandomGen for all
  // GENIE algorithms that use PYTHIA6
  RandomGen::Instance();
#else
  LOG("AGCharmPythia6Hadro2023", pFATAL)
    << "calling GENIE/PYTHIA6 charm hadronization without enabling PYTHIA6";
  gAbortingInErr = true;
  std::exit(1);
#endif

}
//____________________________________________________________________________

bool AGCharmPythia6Hadro2023::HadronizeRemnant (int qrkSyst1, int qrkSyst2,
                                                double WR, TLorentzVector p4R,
                       unsigned int& rpos, TClonesArray * particle_list) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__

  //
  // Run PYTHIA for the hadronization of remnant system
  //
  fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),              1,0); // don't decay pi0
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),  1,1); // decay Delta+
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),  1,1); // decay Delta++
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),  1,1); // decay Delta++
  fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP), 1,1); // decay Delta++
  //   fPythia->SetMDCY(fPythia->Pycomp(kPdgDeltaP),  1,1); // decay Delta+
  //   fPythia->SetMDCY(fPythia->Pycomp(kPdgDeltaPP), 1,1); // decay Delta++
  int ip = 0;
  py2ent_(&ip, &qrkSyst1, &qrkSyst2, &WR); // hadronize

  fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),1,1); // restore

  //-- Get PYTHIA's LUJETS event record
  TClonesArray * pythia_remnants = 0;
  fPythia->GetPrimaries();
  pythia_remnants = dynamic_cast<TClonesArray *>(fPythia->ImportParticles("All"));

  int np = pythia_remnants->GetEntries();
  assert(np>0);

  // PYTHIA performs the hadronization at the *remnant hadrons* centre of mass
  // frame  (not the hadronic centre of mass frame).
  // Boost all hadronic blob fragments to the HCM', fix their mother/daughter
  // assignments and add them to the fragmentation record.

  TVector3 rmnbeta = +1 * p4R.BoostVector(); // boost velocity

  TMCParticle * pythia_remn  = 0; // remnant
  GHepParticle * bremn = 0; // boosted remnant
  TIter remn_iter(pythia_remnants);
  while( (pythia_remn = (TMCParticle *) remn_iter.Next()) ) {

    // insert and get a pointer to inserted object for mods
    bremn = new ((*particle_list)[rpos++]) GHepParticle ( pythia_remn->GetKF(),                // pdg
                                                          GHepStatus_t(pythia_remn->GetKS()),  // status
                                                          pythia_remn->GetParent(),            // first parent
                                                          -1,                                  // second parent
                                                          pythia_remn->GetFirstChild(),        // first daughter
                                                          pythia_remn->GetLastChild(),         // second daughter
                                                          pythia_remn -> GetPx(),              // px
                                                          pythia_remn -> GetPy(),              // py
                                                          pythia_remn -> GetPz(),              // pz
                                                          pythia_remn -> GetEnergy(),          // e
                                                          pythia_remn->GetVx(),                // x
                                                          pythia_remn->GetVy(),                // y
                                                          pythia_remn->GetVz(),                // z
                                                          pythia_remn->GetTime()               // t
                                                          );

    // boost
    bremn -> P4() -> Boost( rmnbeta ) ;

    // handle insertion of charmed hadron
    int jp  = bremn->FirstMother();
    int ifc = bremn->FirstDaughter();
    int ilc = bremn->LastDaughter();

    bremn -> SetFirstMother( (jp  == 0 ?  1 : jp +1) );
    bremn -> SetFirstDaughter ( (ifc == 0 ? -1 : ifc+1) );
    bremn -> SetLastDaughter  ( (ilc == 0 ? -1 : ilc+1) );
  }
  return true;
#else
  LOG("AGCharmPythia6Hadro2023", pFATAL)
    << "calling GENIE/PYTHIA6 charm hadronization without enabling PYTHIA6"
    << " qrkSyst " << qrkSyst1 << "," << qrkSyst2 << " WR " << WR;
  gAbortingInErr = true;
  std::exit(1);
  return false;
#endif


}

//____________________________________________________________________________
