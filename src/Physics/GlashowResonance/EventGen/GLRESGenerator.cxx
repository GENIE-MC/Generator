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
   Add code to decay the off-the-mass-shell W- using PYTHIA6. 
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
#include "Physics/GlashowResonance/EventGen/GLRESGenerator.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
GLRESGenerator::GLRESGenerator() :
EventRecordVisitorI("genie::GLRESGenerator")
{

}
//___________________________________________________________________________
GLRESGenerator::GLRESGenerator(string config) :
EventRecordVisitorI("genie::GLRESGenerator", config)
{

}
//___________________________________________________________________________
GLRESGenerator::~GLRESGenerator()
{

}
//___________________________________________________________________________
void GLRESGenerator::ProcessEventRecord(GHepRecord * event) const
{
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
  char p6frame[10], p6nu[10], p6tgt[10];
  strcpy(p6frame, "CMS"    );
  strcpy(p6nu,    "nu_ebar");
  strcpy(p6tgt,   "e-"     );
  fPythia->Pyinit(p6frame, p6nu, p6tgt, mass);
  fPythia->Pyevnt();
  fPythia->Pylist(1);

  // get LUJETS record
  fPythia->GetPrimaries();
  TClonesArray * pythia_particles =
       (TClonesArray *) fPythia->ImportParticles("All");

  // copy PYTHIA container to a new TClonesArray so as to transfer ownership
  // of the container and of its elements to the calling method
  int np = pythia_particles->GetEntries();
  assert(np>0);
  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", np);
  particle_list->SetOwner(true);

  // Vector defining rotation from LAB to LAB' (z:= \vec{resonance momentum})
  TVector3 unitvq = p4_W.Vect().Unit();

  // Boost velocity LAB' -> Resonance rest frame
  TVector3 beta(0,0,p4_W.P()/p4_W.Energy());

  GHepParticle * p = 0;
  TIter piter(pythia_particles);
  while( (p = (GHepParticle *) piter.Next()) ) {
     int pdgc = p->Pdg();
     int ist  = p->Status();
     if(ist == 1) {
        TLorentzVector p4o(p->Px(), p->Py(), p->Pz(), p->Energy());
        p4o.Boost(beta); 
        TVector3 p3 = p4o.Vect();
        p3.RotateUz(unitvq); 
        TLorentzVector p4(p3,p4o.Energy());
        event->AddParticle(pdgc, kIStStableFinalState, 4,-1,-1,-1, p4, x4);
     }
  }

}
//___________________________________________________________________________
void GLRESGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::LoadConfig(void)
{
 fPythia = TPythia6::Instance();

 // sync GENIE/PYTHIA6 seed number
 RandomGen::Instance();
}
//____________________________________________________________________________

