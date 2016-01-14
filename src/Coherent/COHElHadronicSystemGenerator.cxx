//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include "Coherent/COHElHadronicSystemGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
COHElHadronicSystemGenerator::COHElHadronicSystemGenerator() :
HadronicSystemGenerator("genie::COHElHadronicSystemGenerator")
{

}
//___________________________________________________________________________
COHElHadronicSystemGenerator::COHElHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::COHElHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
COHElHadronicSystemGenerator::~COHElHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void COHElHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system in COH elastic events
//
  int imom = evrec->TargetNucleusPosition();
  int pdgc = evrec->Particle(imom)->Pdg();

  TLorentzVector p4nuc = this->Hadronic4pLAB(evrec);
  TLorentzVector v4(0.,0.,0.,0.);

  evrec->AddParticle(pdgc,kIStStableFinalState, imom,-1,-1,-1, p4nuc, v4);
}
//___________________________________________________________________________

