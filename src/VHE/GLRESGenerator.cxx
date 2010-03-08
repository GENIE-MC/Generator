//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 14, 2009 - CA
   Was first added in v2.5.1 

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "VHE/GLRESGenerator.h"

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
  this -> SelectElectronVelocity (event);
  this -> AddRemnantNucleus      (event); 
  this -> AddResonance           (event);
}
//___________________________________________________________________________
void GLRESGenerator::SelectElectronVelocity(GHepRecord * event) const
{

}
//___________________________________________________________________________
void GLRESGenerator::AddRemnantNucleus(GHepRecord * event) const
{
  GHepParticle * target = event -> TargetNucleus();

  int pdgc = target->Pdg();

  TLorentzVector p4 ( * target->P4() );
  TLorentzVector x4 ( * target->X4() );

  event->AddParticle(pdgc, kIStStableFinalState, 1,-1,-1,-1, p4, x4);
}
//___________________________________________________________________________
void GLRESGenerator::AddResonance(GHepRecord * event) const
{
  GHepParticle * nu = event -> Probe();
  GHepParticle * el = event -> HitElectron();

  assert(nu);
  assert(el);

  TLorentzVector p4_nu (* nu->P4());
  TLorentzVector p4_el (* el->P4());

  TLorentzVector p4_W = p4_nu + p4_el;
  TLorentzVector x4_W(* el->X4());

  event->AddParticle(kPdgWM, kIStStableFinalState, 0,-1,-1,-1, p4_W, x4_W);
}
//___________________________________________________________________________
