//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - November 05, 2007

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 15, 2007 - CA
   The NEUT-cascade wrapper template added in 2.0.1

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "EVGCore/EVGThreadException.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/NEUTCascade.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
NEUTCascade::NEUTCascade() :
EventRecordVisitorI("genie::NEUTCascade")
{

}
//___________________________________________________________________________
NEUTCascade::NEUTCascade(string config) :
EventRecordVisitorI("genie::NEUTCascade",  config)
{

}
//___________________________________________________________________________
NEUTCascade::~NEUTCascade()
{

}
//___________________________________________________________________________
void NEUTCascade::ProcessEventRecord(GHepRecord * event) const
{
  //-- Check that we have an interaction with a nuclear target. If not skip...
  GHepParticle * nucltgt = event->TargetNucleus();
  if (!nucltgt) {
    LOG("NEUTc", pINFO)
       << "No nuclear target. Skipping....";
    return;
  }

#ifdef __GENIE_NEUT_CASCADE_ENABLED__

  //-- Translate the GENIE event record to whatever the NEUT cascade needs
  LOG("NEUTCascade", pDEBUG) 
        << "Translating: GENIE GHepRecord ---> NEUT cascade input";

  //
  // ...

  //-- Run the NEUT cascade code
  LOG("NEUTCascade", pNOTICE) 
        << "************ Running the NEUT Cascade ************";

  //
  // ...

  //-- Get NEUT cascade output & add it to the GENIE event record
  LOG("NEUTCascade", pDEBUG) 
	<< "Copying: NEUT output ---> GENIE GHepRecord";

  //
  // ...

#else
  LOG("NEUTCascade", pFATAL)
       << "\n"
       << "\n******************************************************"
       << "\n*** NEUT IS NOT INSTALLED AT YOUR SYSTEM!          ***"
       << "\n*** Please contact the NEUT team (Y.Hayato et al.) ***"
       << "\n******************************************************"
       << "\n";
  exit(1);
#endif
}
//___________________________________________________________________________
void NEUTCascade::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NEUTCascade::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NEUTCascade::LoadConfig (void)
{
// Access this module's configuration options from its designated Registry 
// and pass them to the actual NEUT cascade code


}
//___________________________________________________________________________

