//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - Feb 15, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 15, 2009 - CA
   This class was first added in version 2.5.1.

*/
//____________________________________________________________________________

#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Diffractive/EventGen/DFRInteractionListGenerator.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
DFRInteractionListGenerator::DFRInteractionListGenerator() :
InteractionListGeneratorI("genie::DFRInteractionListGenerator")
{

}
//___________________________________________________________________________
DFRInteractionListGenerator::DFRInteractionListGenerator(string config):
InteractionListGeneratorI("genie::DFRInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DFRInteractionListGenerator::~DFRInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DFRInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  InteractionType_t inttype;
  if      (fIsCC) inttype = kIntWeakCC;
  else if (fIsNC) inttype = kIntWeakNC;
  else {
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  int nupdg = init_state.ProbePdg();
  if( !pdg::IsNeutrino(nupdg) && !pdg::IsAntiNeutrino(nupdg) ) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  // at least for now, we don't have any models
  // that provide diffractive scattering from nuclear targets--
  // just free protons
  const Target & target = init_state.Tgt();
  if( !target.IsProton() ) {
     LOG("IntLst", pINFO)
       << "Target is not free proton -- Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }


  // The code below does way more work than necessary
  // to assign the final state when scattering
  // is just for free protons (as ensured by the checks above).
  // However, it does the correct thing,
  // and we leave it here in the event that someday we
  // introduce a diffractive model that can handle other targets.
  InteractionList * intlist = new InteractionList;

  ProcessInfo proc_info(kScDiffractive, inttype);

  bool hasP = (init_state.Tgt().Z() > 0);
  bool hasN = (init_state.Tgt().N() > 0);

  int hit_nucleon[2] = {kPdgProton, kPdgNeutron};

  for(int i=0; i<2; i++) {

    Interaction * interaction = new Interaction(init_state, proc_info);

    int nuc = hit_nucleon[i];

    if(nuc == kPdgProton  && !hasP) continue;
    if(nuc == kPdgNeutron && !hasN) continue;

    if(fIsCC) {
      if(pdg::IsNeutrino(nupdg)) {
        //v N -> l- N pi+
        interaction->ExclTagPtr()->SetNPions(1,0,0);  
      } else {
        //vbar N -> l+ N pi-
        interaction->ExclTagPtr()->SetNPions(0,0,1); 
      }
    }
    else {
      //v N -> v N pi0
      interaction->ExclTagPtr()->SetNPions(0,1,0);
    }

    if (nuc == kPdgProton)
      interaction->ExclTagPtr()->SetNProtons(1);
    else if (nuc == kPdgNeutron)
      interaction->ExclTagPtr()->SetNNeutrons(1);

    interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(nuc);
    intlist->push_back(interaction);
  }

  return intlist;
}
//___________________________________________________________________________
void DFRInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DFRInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DFRInteractionListGenerator::LoadConfigData(void)
{
	GetParamDef( "is-CC",fIsCC, false  ) ;
	GetParamDef( "is-NC",fIsNC, false  ) ;
}
//____________________________________________________________________________


