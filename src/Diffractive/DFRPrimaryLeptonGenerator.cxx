//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - Feb 15, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 15, 2009 - CA
   This class was first added in version 2.5.1.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Diffractive/DFRPrimaryLeptonGenerator.h"
#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
DFRPrimaryLeptonGenerator::DFRPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::DFRPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
DFRPrimaryLeptonGenerator::DFRPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::DFRPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
DFRPrimaryLeptonGenerator::~DFRPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void DFRPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in DIS events

  // no modification is required to the std implementation
  PrimaryLeptonGenerator::ProcessEventRecord(evrec);

  if(evrec->FinalStatePrimaryLepton()->IsOffMassShell()) {
    LOG("LeptonicVertex", pERROR)
               << "*** Selected kinematics lead to off mass shell lepton!";
     evrec->EventFlags()->SetBitNumber(kLeptoGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("E<m for final state lepton");
     exception.SwitchOnFastForward();
     throw exception;
  }
}
//___________________________________________________________________________
