//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Algorithm/AlgConfigPool.h"
#include "EVGModules/CorrelatedNucleonGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
CorrelatedNucleonGenerator::CorrelatedNucleonGenerator() :
EventRecordVisitorI("genie::CorrelatedNucleonGenerator")
{

}
//___________________________________________________________________________
CorrelatedNucleonGenerator::CorrelatedNucleonGenerator(string config) :
EventRecordVisitorI("genie::CorrelatedNucleonGenerator", config)
{

}
//___________________________________________________________________________
CorrelatedNucleonGenerator::~CorrelatedNucleonGenerator()
{

}
//___________________________________________________________________________
void CorrelatedNucleonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(!fSimulateCorrelN) return;

  Interaction *  interaction = evrec       -> Summary();
  InitialState * init_state  = interaction -> InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  // do nothing for non-nuclear targets
  if(!tgt->IsNucleus()) return;



}
//___________________________________________________________________________
void CorrelatedNucleonGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CorrelatedNucleonGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CorrelatedNucleonGenerator::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fSimulateCorrelN = fConfig->GetBoolDef(
              "Enable", gc->GetBool("CorrelNN-Enable"));
  fMomentumThr  = fConfig->GetDoubleDef(
              "MomentumThreshold", gc->GetBool("CorrelNN-MomentumThreshold"));
}
//____________________________________________________________________________

