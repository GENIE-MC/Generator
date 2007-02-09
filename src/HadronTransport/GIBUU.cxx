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
#include "HadronTransport/GIBUU.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
GIBUU::GIBUU() :
EventRecordVisitorI("genie::GIBUU")
{

}
//___________________________________________________________________________
GIBUU::GIBUU(string config) :
EventRecordVisitorI("genie::GIBUU", config)
{

}
//___________________________________________________________________________
GIBUU::~GIBUU()
{

}
//___________________________________________________________________________
void GIBUU::ProcessEventRecord(GHepRecord * event) const
{
  //-- Check that we have an interaction with a nuclear target. If not skip...
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("GIBUU", pINFO)
       << "No nuclear target found - GIBUU will not be called....";
    return;
  }

  //-- Translate GENIE GHepRecord to whatever FLUKA needs
  LOG("GIBUU", pDEBUG) << "Translating: GENIE GHepRecord ---> GIBUU input";

  //
  // ...

  //-- Run GIBUU
  LOG("GIBUU", pNOTICE) << "************ Running GIBUU ************";

  //
  // ...

  //-- Get FLUKA output & add it to GENIE GHepRecord
  LOG("GIBUU", pDEBUG) << "Copying: GIBUU output ---> GENIE GHepRecord";

  //
  // ...

}
//___________________________________________________________________________
void GIBUU::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void GIBUU::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void GIBUU::LoadConfig (void)
{
// Access this module's configuration options from its designated Registry
// and pass them to the actual GIBUU code
   
}
//___________________________________________________________________________


