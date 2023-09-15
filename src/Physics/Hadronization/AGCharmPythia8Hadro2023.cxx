//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
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
#include "Physics/Hadronization/AGCharmPythia8Hadro2023.h"
#include "Physics/Hadronization/FragmentationFunctionI.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
AGCharmPythia8Hadro2023::AGCharmPythia8Hadro2023() :
AGCharmPythiaBaseHadro2023("genie::AGCharmPythia8Hadro2023")
{
  this->Initialize();
}
//____________________________________________________________________________
AGCharmPythia8Hadro2023::AGCharmPythia8Hadro2023(string config) :
AGCharmPythiaBaseHadro2023("genie::AGCharmPythia8Hadro2023", config)
{
  this->Initialize();
}
//____________________________________________________________________________
AGCharmPythia8Hadro2023::~AGCharmPythia8Hadro2023()
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  delete fPythia;
#endif
}
//____________________________________________________________________________
void AGCharmPythia8Hadro2023::Initialize(void) const
{
  AGCharmPythiaBaseHadro2023::Initialize();
#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = new Pythia8::Pythia();

  fPythia->readString("ProcessLevel:all = off");
  fPythia->readString("Print:quiet = on");

  // sync GENIE and PYTHIA8 seeds
  RandomGen * rnd = RandomGen::Instance();
  long int seed = rnd->GetSeed();
  fPythia->readString("Random:setSeed = on");
  fPythia->settings.mode("Random:seed", seed);
  LOG("AGCharmPythia8Hadro2023", pINFO)
    << "PYTHIA8 seed = " << fPythia->settings.mode("Random:seed");
  fPythia->init();

#else
  LOG("AGCharmPythia8Hadro2023", pFATAL)
    << "calling GENIE/PYTHIA8 charm hadronization without enabling PYTHIA8";
  gAbortingInErr = true;
  std::exit(1);
#endif

}
//____________________________________________________________________________

bool AGCharmPythia8Hadro2023::HadronizeRemnant (int qrkSyst1, int qrkSyst2,
                                                double WR, TLorentzVector p4R,
                       unsigned int& rpos, TClonesArray * particle_list) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  /*
  //
  // Run PYTHIA6 for the hadronization of remnant system
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
  return pythia_remnants;
  */
  fPythia->particleData.mayDecay(kPdgPi0,false); // don't decay pi0
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaM,true); // decay Delta+
  fPythia->particleData.mayDecay(kPdgP33m1232_Delta0,true); // decay Delta++
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaP,true); // decay Delta++
  fPythia->particleData.mayDecay(kPdgP33m1232_DeltaPP,true); // decay Delta++

  LOG("AGCharmPythia8Hadro2023", pFATAL)
    << "HadronizeRemnant not yet converted for PYTHIA8"
    << " qrkSyst " << qrkSyst1 << "," << qrkSyst2 << " WR " << WR;
  gAbortingInErr = true;
  std::exit(1);


  fPythia->particleData.mayDecay(kPdgPi0,true); // restore

#else
  LOG("AGCharmPythia8Hadro2023", pFATAL)
    << "calling GENIE/PYTHIA8 charm hadronization without enabling PYTHIA8"
    << " qrkSyst " << qrkSyst1 << "," << qrkSyst2 << " WR " << WR;
  gAbortingInErr = true;
  std::exit(1);
#endif


}

//____________________________________________________________________________
