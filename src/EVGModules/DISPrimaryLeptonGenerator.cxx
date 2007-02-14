//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "EVGCore/EVGThreadException.h"
#include "EVGModules/DISPrimaryLeptonGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"

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
