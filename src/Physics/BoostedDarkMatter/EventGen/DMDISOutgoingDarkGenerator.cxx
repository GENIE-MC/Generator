//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/BoostedDarkMatter/EventGen/DMDISOutgoingDarkGenerator.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
DMDISOutgoingDarkGenerator::DMDISOutgoingDarkGenerator() :
OutgoingDarkGenerator("genie::DMDISOutgoingDarkGenerator")
{

}
//___________________________________________________________________________
DMDISOutgoingDarkGenerator::DMDISOutgoingDarkGenerator(string config) :
OutgoingDarkGenerator("genie::DMDISOutgoingDarkGenerator", config)
{

}
//___________________________________________________________________________
DMDISOutgoingDarkGenerator::~DMDISOutgoingDarkGenerator()
{

}
//___________________________________________________________________________
void DMDISOutgoingDarkGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton in DIS events

  // no modification is required to the std implementation
  OutgoingDarkGenerator::ProcessEventRecord(evrec);

  if(evrec->FinalStatePrimaryLepton()->IsOffMassShell()) {
    LOG("LeptonicVertex", pERROR)
               << "*** Selected kinematics lead to off mass shell dark matter!";
     evrec->EventFlags()->SetBitNumber(kLeptoGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("E<m for final state dark matter");
     exception.SwitchOnFastForward();
     throw exception;
  }
}
//___________________________________________________________________________
