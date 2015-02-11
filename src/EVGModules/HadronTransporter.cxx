//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 15, 2009 - CA
   IsFake() is no longer available in GHepParticle.Use pdg::IsPseudoParticle() 
 @ Nov 17, 2011 - CA
   Removed unused GenerateVertex() method. This is handled by another module.

*/
//____________________________________________________________________________

#include <cstdlib>

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "EVGModules/HadronTransporter.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
HadronTransporter::HadronTransporter() :
EventRecordVisitorI("genie::HadronTransporter")
{

}
//___________________________________________________________________________
HadronTransporter::HadronTransporter(string config) :
EventRecordVisitorI("genie::HadronTransporter", config)
{

}
//___________________________________________________________________________
HadronTransporter::~HadronTransporter()
{

}
//___________________________________________________________________________
void HadronTransporter::ProcessEventRecord(GHepRecord * evrec) const
{
  // Return if the neutrino was not scatterred off a nuclear target
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("HadTransp", pINFO) 
           << "No nuclear target found - Hadron transporter exits";
    return;
  }

  // If rescattering is turned off but the interaction was on a nuclear
  // target then simply transfer the hadrons outside the nucleus inhibiting
  // any rescattering
  if(!fEnabled) {
    LOG("HadTransp", pNOTICE) 
             << "*** Intranuclear rescattering has been turned off";
    this->TransportInTransparentNuc(evrec);
    return;
  }

  // Use the specified intrsnuclear rescattering model
  LOG("HadTransp", pINFO)  << "Calling the selected hadron transport MC";
  fHadTranspModel->ProcessEventRecord(evrec);
}
//___________________________________________________________________________
void HadronTransporter::TransportInTransparentNuc(GHepRecord * evrec) const
{
// Transport all hadrons assuming a transparent nucleus - used when the
// realistic hadron transport is tuned off.

  LOG("HadTransp", pNOTICE) 
     << "Getting the nucleons out of the nucleus as if it was transparent";

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr=-1;

  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;
    assert(p);

    // Check whether the particle needs rescattering, otherwise skip it

    bool had_in_nuc = (p->Status() == kIStHadronInTheNucleus);
    if(!had_in_nuc) continue;

    LOG("HadTransp", pINFO)
       << "Transporting " << p->Name() << " out of the nuclear target";

    // move it outside the nucleus
    GHepParticle * cp = new GHepParticle(*p); // create a clone

    cp->SetFirstMother(icurr);                // clone's mother
    cp->SetStatus(kIStStableFinalState);      // mark it & done with it

    evrec->AddParticle(*cp); // add it at the event record

    //LOG("HadronTransporter", pDEBUG) << *evrec;
  }
}
//___________________________________________________________________________
void HadronTransporter::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void HadronTransporter::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//___________________________________________________________________________
void HadronTransporter::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fHadTranspModel = 0; 
  fEnabled = fConfig->GetBoolDef("enable", gc->GetBool("HadronTransp-Enable")); 

  LOG("HadTransp", pDEBUG) 
       << "Hadron transport was " << ((fEnabled) ? "" : "not ") << " enabled";
  if(fEnabled) {
     RgAlg hadtransp_model = 
        fConfig->GetAlgDef("model", gc->GetAlg("HadronTransp-Model"));
     LOG("HadTransp", pDEBUG) 
         << "Loading the hadron transport model: " << hadtransp_model;

     // Note: this relies on the fact that if a missing registry entry was
     // found above then it would be filled up by the default choice so that
     // when the registry is looked-up again the same key would point to the 
     // correct value. Probably it would be better to convert the RgAlg to AlgId
     // and query the AlgFactory directly. 
     fHadTranspModel = 
         dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("model"));
     assert(fHadTranspModel);
  }
}
//___________________________________________________________________________


