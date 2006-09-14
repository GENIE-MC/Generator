//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, CCLRC, Rutherford Lab
         September 14, 2006

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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
  if(!fRescatON) {
    LOG("HadTransp", pNOTICE) 
                     << "*** Intranuclear rescattering has been turned off";
    this->GenerateVertex(evrec);
    this->TransportInTransparentNuc(evrec);
    return;
  }

  // Use the specified intrsnuclear rescattering model
  fHadTranspModel->ProcessEventRecord(evrec);
}
//___________________________________________________________________________
void HadronTransporter::GenerateVertex(GHepRecord * evrec) const
{
// generate a vtx and set it to all GHEP physical particles

  GHepParticle * nucltgt = evrec->TargetNucleus();
  assert(nucltgt);

  RandomGen * rnd = RandomGen::Instance();

  int    A     = nucltgt->A();
  double Ro    = nuclear::Radius(A); //fm
  double R     = Ro * rnd->RndFsi().Rndm();
  double cos9  = -1. + 2. * rnd->RndFsi().Rndm();    
  double sin9  = TMath::Sqrt(1.-cos9*cos9);   
  double fi    = 2 * kPi * rnd->RndFsi().Rndm();
  double cosfi = TMath::Cos(fi);
  double sinfi = TMath::Sin(fi);

  TVector3 vtx(R*sin9*cosfi,R*sin9*sinfi,R*cos9);

  LOG("HadTransp", pINFO) << "Vtx (in fm) = " << print::Vec3AsString(&vtx);

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(p->IsFake()) continue;
    p->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);
  }
}
//___________________________________________________________________________
void HadronTransporter::TransportInTransparentNuc(GHepRecord * evrec) const
{
// transport all hadrons assuming a transparent nucleus

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

  fRescatON = fConfig->GetBoolDef("rescatter",gc->GetBool("InuclRescat-On")); 

  if(fRescatON) {
    string name   = fConfig->GetStringDef(
        "rescattering-model-name",   gc->GetString("InuclRescat-ModelName")); 
    string config = fConfig->GetStringDef(
      "rescattering-model-config", gc->GetString("InuclRescat-ModelConfig"));

    AlgFactory * algf = AlgFactory::Instance();
    fHadTranspModel = dynamic_cast<const EventRecordVisitorI *>(
					    algf->GetAlgorithm(name,config));
  }
}
//___________________________________________________________________________


