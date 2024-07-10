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

#endif
}
//____________________________________________________________________________
void AGCharmPythia8Hadro2023::Initialize(void) const
{
  AGCharmPythiaBaseHadro2023::Initialize();
#ifdef __GENIE_PYTHIA8_ENABLED__
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  gPythia->readString("ProcessLevel:all = off");
  gPythia->readString("Print:quiet = on");

  // sync GENIE and PYTHIA8 seeds
  RandomGen * rnd = RandomGen::Instance();
  long int seed = rnd->GetSeed();
  gPythia->readString("Random:setSeed = on");
  gPythia->settings.mode("Random:seed", seed);
  LOG("AGCharmPythia8Hadro2023", pINFO)
    << "PYTHIA8 seed = " << gPythia->settings.mode("Random:seed");
  gPythia->init();

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
  Pythia8::Pythia* gPythia = Pythia8Singleton::Instance()->Pythia8();

  /*
  //
  // Run PYTHIA6 for the hadronization of remnant system
  //
  gPythia->SetMDCY(gPythia->Pycomp(kPdgPi0),              1,0); // don't decay pi0
  gPythia->SetMDCY(gPythia->Pycomp(kPdgP33m1232_DeltaM),  1,1); // decay Delta+
  gPythia->SetMDCY(gPythia->Pycomp(kPdgP33m1232_Delta0),  1,1); // decay Delta++
  gPythia->SetMDCY(gPythia->Pycomp(kPdgP33m1232_DeltaP),  1,1); // decay Delta++
  gPythia->SetMDCY(gPythia->Pycomp(kPdgP33m1232_DeltaPP), 1,1); // decay Delta++
  //   gPythia->SetMDCY(gPythia->Pycomp(kPdgDeltaP),  1,1); // decay Delta+
  //   gPythia->SetMDCY(gPythia->Pycomp(kPdgDeltaPP), 1,1); // decay Delta++
  int ip = 0;
  py2ent_(&ip, &qrkSyst1, &qrkSyst2, &WR); // hadronize

  gPythia->SetMDCY(gPythia->Pycomp(kPdgPi0),1,1); // restore

  //-- Get PYTHIA's LUJETS event record
  TClonesArray * pythia_remnants = 0;
  gPythia->GetPrimaries();
  pythia_remnants = dynamic_cast<TClonesArray *>(gPythia->ImportParticles("All"));
  return pythia_remnants;
  */
  gPythia->particleData.mayDecay(kPdgPi0,false); // don't decay pi0
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaM,true); // decay Delta+
  gPythia->particleData.mayDecay(kPdgP33m1232_Delta0,true); // decay Delta++
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaP,true); // decay Delta++
  gPythia->particleData.mayDecay(kPdgP33m1232_DeltaPP,true); // decay Delta++

  gPythia->event.reset();

  double m1 = gPythia->particleData.m0(qrkSyst1);
  double m2 = gPythia->particleData.m0(qrkSyst2);

  LOG("AGCharmPythia8Hadro2023", pINFO) // debug
    << " qrkSyst1 " << qrkSyst1 << " m1 " << m1
    << " qrkSyst2 " << qrkSyst2 << " m2 " << m2
    << " WR " << WR << " p4R " << genie::utils::print::P4AsString(&p4R);

  double pz1cm = 0.5 * Pythia8::sqrtpos(
    (WR + m1 + m2)*(WR - m1 -m2)*(WR - m1 + m2)*(WR + m1 - m2) ) / WR;

  double pz2cm = - pz1cm;
  double e1    = sqrt(m1*m1 + pz1cm*pz1cm);
  double e2    = sqrt(m2*m2 + pz2cm*pz2cm);

  LOG("AGCharmPythia8Hadro2023", pINFO) // debug
    << " qrkSyst1 pz=" << pz1cm << ", E=" << e1
    << " qrkSyst2 pz=" << pz2cm << ", E=" << e2;

  // status codes and anti/collor tags must complement each other
  // id, status, int col, int acol, px,py,pz,E,m,scale=0,pol=9
  gPythia->event.append(qrkSyst1,23,101,  0,0.,0.,pz1cm,e1,m1);
  gPythia->event.append(qrkSyst2,23,  0,101,0.,0.,pz2cm,e2,m2);
  //gPythia->event.list();

  // Generating next pythia8 event
  gPythia->next();

  //gPythia->event.list();
  //gPythia->stat();

  // get the LUJETss record
  Pythia8::Event &fEvent = gPythia->event;
  int np = fEvent.size();
  assert(np>0);

  // Vector defining rotation from LAB to LAB' (z:= \vec{p4R}J)
  TVector3 unitvq = p4R.Vect().Unit();

  // Boost velocity LAB' -> HCM
  TVector3 beta(0,0,p4R.P()/p4R.Energy());

  GHepStatus_t istfin = kIStHadronInTheNucleus; /// should be testing on is_nucleus
  // from interaction->InitState().Tgt().IsNucleus();

  int mom = -1; // event->FinalStateHadronicSystemPosition();
  const TLorentzVector vtx(0,0,0,0);

  for (int i=0; i<np; ++i) {
    if (fEvent[i].id() == 90) continue;   // ignore (system) pseudoparticle

    int particle_pdg_code      = fEvent[i].id();
    int pythia_particle_status = fEvent[i].status();

    GHepParticle* bremn = 0; // boosted remnant

    // sanity check
    if (pythia_particle_status > 0 ) {
      if ( pdg::IsQuark  (particle_pdg_code ) ||
           pdg::IsDiQuark(particle_pdg_code )    )
        {
          LOG("AGCharmPythia8Hadro2023", pERROR)
            << "Hadronization failed! Bare quarks appear in final state!";
          return false;
        }
    }
    // copy the initial quark systemcs (satus = -23) and all undecayed particles
    // ignore particle we asked PYTHIA to decay but record their decay products
    bool copy = (pythia_particle_status==-23) ||
      (pythia_particle_status > 0);
    if (copy) {
      // the fragmenation products are generated in the remnant hadronic CM frame
      // where the z>0 asix is the \vec{p4R} direction.
      // for each particle returned by the hadronizer:
      //   -- boost it back to LAB' frame {z:\vec{p4T}} / doesn't affect pT
      //   -- rotate its 3-momentum from LAB' to LAB
      TLorentzVector p4o(fEvent[i].px(),fEvent[i].py(),fEvent[i].pz(),fEvent[i].e());
      p4o.Boost(beta);
      TVector3 p3 = p4o.Vect();
      p3.RotateUz(unitvq);
      TLorentzVector p4(p3,p4o.Energy());


      //!!!!
      istfin = GHepStatus_t(1);

      // set the proper GENIE status according to a number of things:
      GHepStatus_t ist = (pythia_particle_status > 0) ?
        istfin : kIStDISPreFragmHadronicState;

      // handle gammas, and leptons that might come from internal pythia decay
      // mark them as final state particles
      bool is_gamma = (particle_pdg_code == kPdgGamma);
      bool is_nu    = pdg::IsNeutralLepton(particle_pdg_code);
      bool is_lchg  = pdg::IsChargedLepton(particle_pdg_code);
      bool not_hadr = is_gamma || is_nu || is_lchg;
      if (not_hadr) { ist = kIStStableFinalState; }

      int mother1 = mom;
      int mother2 = -1;

      int daughter1 = -1;
      int daughter2 = -1;

      // create GHepParticle

      bremn = new ((*particle_list)[rpos++]) GHepParticle(
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

      //std::cout << "bremn " << rpos-1 << "||" << *bremn << std::endl;

    } // if (copy)
  } // loop over particles in pythia event record

  gPythia->particleData.mayDecay(kPdgPi0,true); // restore

  return true;

#else
  LOG("AGCharmPythia8Hadro2023", pFATAL)
    << "calling GENIE/PYTHIA8 charm hadronization without enabling PYTHIA8"
    << " qrkSyst " << qrkSyst1 << "," << qrkSyst2 << " WR " << WR;
  gAbortingInErr = true;
  std::exit(1);
#endif


}

//____________________________________________________________________________
