//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 15, 2009 - CA
   IsFake() is no longer available in GHepParticle.Use pdg::IsPseudoParticle() 
 @ Nov 17, 2011 - CA
   Removed unused GenerateVertex() method. This is handled by another module.

*/
//____________________________________________________________________________

#include <cstdlib>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/HadronTransport/HadronTransporter.h"
#include "Physics/HadronTransport/CascadeReweight.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

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

  // Process Event record for the cascade reweight if required:
  if(fCascadeReweight) fCascadeReweight->ProcessEventRecord(evrec);

  return ;
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

  Registry * algos = AlgConfigPool::Instance() -> GlobalParameterList() ;
  Registry r( "HadronTransporter_specific", false ) ;

  r.Set("HadronTransp-Enable", algos -> GetBool("HadronTransp-Enable") ) ;
  r.Set("HadronTransp-Model",  algos -> GetAlg("HadronTransp-Model")   ) ;
  if( algos->Exists("CascadeReweightAlg") ) {
    r.Set("CascadeReweightAlg",  algos -> GetAlg("CascadeReweightAlg")   ) ;
  }

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//___________________________________________________________________________
void HadronTransporter::LoadConfig(void)
{
  bool good_config = true ;

  fHadTranspModel = 0;
  GetParam("HadronTransp-Enable", fEnabled ) ;

  LOG("HadTransp", pDEBUG) 
       << "Hadron transport was " << ((fEnabled) ? "" : "not ") << " enabled";
  if(fEnabled) {
     RgAlg hadtransp_model ;
     GetParam( "HadronTransp-Model", hadtransp_model ) ;
     LOG("HadTransp", pDEBUG) 
         << "Loading the hadron transport model: " << hadtransp_model;

     fHadTranspModel = 
         dynamic_cast<const EventRecordVisitorI *> ( this -> SubAlg("HadronTransp-Model") );

     if(!fHadTranspModel){
       good_config = false ;
       LOG("HadronTransporter", pERROR) << "The SubAlg HadronTransp-Model is not a XSecAlgorithmI";
     }

     // Read optional CascadeReweight EventRecordVisitorI
     fCascadeReweight = nullptr; 
     if( GetConfig().Exists("CascadeReweightAlg") ) {
       fCascadeReweight = dynamic_cast<const CascadeReweight *> ( this->SubAlg("CascadeReweightAlg") );
       if( !fCascadeReweight ) {
	 good_config = false ; 
	 LOG("HadronTransporter", pERROR) << "The required CascadeReweightAlg does not exist. AlgID is : " << SubAlg("CascadeReweightAlg")->Id() ;
       }
     }

     if( ! good_config ) {
       LOG("HadronTransporter", pERROR) << "Configuration has failed.";
       exit(78) ;
     }
  }
}
//___________________________________________________________________________


