//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - December 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/FLUKA.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
FLUKA::FLUKA() :
EventRecordVisitorI("genie::FLUKA")
{

}
//___________________________________________________________________________
FLUKA::FLUKA(string config) :
EventRecordVisitorI("genie::FLUKA",  config)
{

}
//___________________________________________________________________________
FLUKA::~FLUKA()
{

}
//___________________________________________________________________________
void FLUKA::ProcessEventRecord(GHepRecord * event) const
{
  //-- Check that we have an interaction with a nuclear target. If not skip...
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("FLUKA", pINFO)
       << "No nuclear target found - FLUKA will not be called....";
    return;
  }

  //-- Translate GENIE GHepRecord to whatever FLUKA needs
  LOG("FLUKA", pDEBUG) << "Translating: GENIE GHepRecord ---> FLUKA input";

  //
  // ...

  //-- Run FLUKA
  LOG("FLUKA", pNOTICE) << "************ Running FLUKA ************";

  //
  // ...

  //-- Get FLUKA output & add it to GENIE GHepRecord
  LOG("FLUKA", pDEBUG) << "Copying: FLUKA output ---> GENIE GHepRecord";

  //
  // ...

}
//___________________________________________________________________________
void FLUKA::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void FLUKA::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void FLUKA::LoadConfig (void)
{
// Access this module's configuration options from its designated Registry 
// and pass them to the actual FLUKA code


}
//___________________________________________________________________________

