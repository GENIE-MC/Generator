//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new DIS package from its previous location (EVGModules).

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/DeepInelastic/EventGen/DISPrimaryLeptonGenerator.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
DISPrimaryLeptonGenerator::DISPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::DISPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
DISPrimaryLeptonGenerator::DISPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::DISPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
DISPrimaryLeptonGenerator::~DISPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void DISPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
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
