//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Igor Kakorin <kakorin@jinr.ru>
 Joint Institute for Nuclear Research
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>


#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/Common/NormGenerator.h"

using namespace genie;

//___________________________________________________________________________
NormGenerator::NormGenerator():
EventRecordVisitorI("genie::NormGenerator")
{

}
//___________________________________________________________________________
NormGenerator::NormGenerator(string config):
EventRecordVisitorI("genie::NormGenerator")
{

}
//___________________________________________________________________________
NormGenerator::~NormGenerator()
{

}
//___________________________________________________________________________
void NormGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction -> InitState();

  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //bool isHeavyNucleus = tgt->A()>=3;
  GHepParticle * probe  = evrec->Probe();
  int pdg_probe = probe->Pdg();
  int iprobe = evrec->ProbePosition();
  TLorentzVector * p4v = probe->GetP4();
  TLorentzVector * vtx = probe->X4();
  
  GHepParticle * nucltgt = evrec->TargetNucleus();
  int pdg_tgt = nucltgt->Pdg();
  int inucltgt = evrec->TargetNucleusPosition();
  TLorentzVector * p4_tgt = nucltgt->GetP4();
  TLorentzVector * vtx_tgt = nucltgt->X4();

  
  evrec->AddParticle(pdg_probe, kIStStableFinalState, iprobe,-1,-1,-1, *p4v, *vtx);
  evrec->AddParticle(pdg_tgt, kIStStableFinalState, inucltgt,-1,-1,-1, *p4_tgt, *vtx_tgt);
  
   //-- Determine the status code
  //const Target & tgt = interaction->InitState().Tgt();


  // update the interaction summary
  evrec->Summary()->KinePtr()->SetFSLeptonP4(*p4v);
  
  double xsec = fXSecModel->XSec(interaction, kPSfE);
  
  evrec->SetDiffXSec(xsec, kPSfE);
  
  return;
}
//___________________________________________________________________________
void NormGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NormGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NormGenerator::LoadConfig(void)
{
}
//____________________________________________________________________________
