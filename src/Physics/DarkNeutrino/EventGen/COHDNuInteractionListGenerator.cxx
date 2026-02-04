//____________________________________________________________________________
/*
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
  University of Sussex

  Costas Andreopoulos <c.andreopoulos \at cern.ch>
  University of Liverpool
*/
//____________________________________________________________________________

#include "Physics/DarkNeutrino/EventGen/COHDNuInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
COHDNuInteractionListGenerator::COHDNuInteractionListGenerator() :
InteractionListGeneratorI("genie::COHDNuInteractionListGenerator")
{

}
//___________________________________________________________________________
COHDNuInteractionListGenerator::COHDNuInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::COHDNuInteractionListGenerator", config)
{

}
//___________________________________________________________________________
COHDNuInteractionListGenerator::~COHDNuInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * COHDNuInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
      << "InitialState = " << init_state.AsString();

  int probe_pdg = init_state.ProbePdg();
  bool isnu = pdg::IsNeutrino(probe_pdg) || pdg::IsAntiNeutrino(probe_pdg);
  if( !isnu) {
     // shouldn't happen... warn
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }
  const Target & target = init_state.Tgt();
  if(!target.IsNucleus()) {
     // happens as this code is also indiscriminately both for free-nucleon and
     // nuclear targets - don't warn
     LOG("IntLst", pINFO)
       << "Not a nuclear target! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  InteractionList * intlist = new InteractionList;

  ProcessInfo proc_info( kScCoherentElastic, kIntDarkNC );
  Interaction * interaction = new Interaction( init_state, proc_info);

  intlist->push_back(interaction);

  return intlist;
}
//___________________________________________________________________________
void COHDNuInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void COHDNuInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void COHDNuInteractionListGenerator::LoadConfigData(void)
{

}
//____________________________________________________________________________
