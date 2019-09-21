//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2009 - CA
   Was first added in v2.5.1 
 @ Apr 24, 2010 - CA
   Add code to decay the off-the-mass-shell W- using PYTHIA8. 
   First complete version of the GLRES event thread.
*/
//____________________________________________________________________________

#include <cstring>

#include <RVersion.h>
#include <TClonesArray.h>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/GlashowResonance/EventGen/GLRESP8Generator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
GLRESP8Generator::GLRESP8Generator() :
EventRecordVisitorI("genie::GLRESP8Generator")
{

}
//___________________________________________________________________________
GLRESP8Generator::GLRESP8Generator(string config) :
EventRecordVisitorI("genie::GLRESP8Generator", config)
{

}
//___________________________________________________________________________
GLRESP8Generator::~GLRESP8Generator()
{

}
//___________________________________________________________________________
void GLRESP8Generator::ProcessEventRecord(GHepRecord * event) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  GHepParticle * nu = event -> Probe();
  GHepParticle * el = event -> HitElectron();
  assert(nu);
  assert(el);

  // Pick electron velocity
  //
  //... neglect for now

  //
  // Add remnant nucleus 
  //
  GHepParticle * target = event -> TargetNucleus();
  if(target) {
    int pdgc = target->Pdg();
    TLorentzVector p4 ( * target->P4() );
    TLorentzVector x4 ( * target->X4() );
    event->AddParticle(pdgc, kIStStableFinalState, 1,-1,-1,-1, p4, x4);
  }

  //
  // Add resonance
  //

  TLorentzVector p4_nu (* nu->P4());
  TLorentzVector p4_el (* el->P4());

  TLorentzVector x4(* el->X4());

  TLorentzVector p4_W = p4_nu + p4_el;

  event->AddParticle(kPdgWM, kIStDecayedState, 0,-1,-1,-1, p4_W, x4);

  //
  // Decay resonance and add decay products
  //

  double mass = p4_W.M();
  fPythia->Pythia8()->readString("Print:quiet           = on");
  fPythia->Pythia8()->readString("PDF:lepton            = off");
  fPythia->Pythia8()->readString("WeakBosonExchange:all = on");

  fPythia->Pythia8()->settings.mode("Beams:idA", -12); // nu_ebar
  fPythia->Pythia8()->settings.mode("Beams:idB", 11); // e-
  fPythia->Pythia8()->settings.mode("Beams:frameType", 1);
  fPythia->Pythia8()->settings.parm("Beams:eCM", mass);

  fPythia->Pythia8()->next();
  fPythia->Pythia8()->event.list();

  Pythia8::Event &fEvent = fPythia->Pythia8()->event;
  int numpart = fEvent.size();
  assert(numpart>0);

  // Vector defining rotation from LAB to LAB' (z:= \vec{resonance momentum})
  TVector3 unitvq = p4_W.Vect().Unit();

  // Boost velocity LAB' -> Resonance rest frame
  TVector3 beta(0,0,p4_W.P()/p4_W.Energy());

  for (int i = 1; i < numpart; ++i) {
    int pdgc = fEvent[i].id();
    int ist  = fEvent[i].status();
    if(ist > 0) {
        TLorentzVector p4o(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
        p4o.Boost(beta); 
        TVector3 p3 = p4o.Vect();
        p3.RotateUz(unitvq); 
        TLorentzVector p4(p3,p4o.Energy());
        event->AddParticle(pdgc, kIStStableFinalState, 4,-1,-1,-1, p4, x4);
     }
  }
#endif

}
//___________________________________________________________________________
void GLRESP8Generator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESP8Generator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESP8Generator::LoadConfig(void)
{
#ifdef __GENIE_PYTHIA8_ENABLED__
 fPythia = Pythia8Singleton::Instance();

 // sync GENIE/PYTHIA8 seed number
 RandomGen::Instance();
#endif
}
//____________________________________________________________________________

