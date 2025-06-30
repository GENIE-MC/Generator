//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Multinucleon/EventGen/MECInteractionListGenerator.h"

using namespace genie;

//___________________________________________________________________________
MECInteractionListGenerator::MECInteractionListGenerator() :
InteractionListGeneratorI("genie::MECInteractionListGenerator")
{

}
//___________________________________________________________________________
MECInteractionListGenerator::MECInteractionListGenerator(string config):
InteractionListGeneratorI("genie::MECInteractionListGenerator", config)
{

}
//___________________________________________________________________________
MECInteractionListGenerator::~MECInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList *
  MECInteractionListGenerator::CreateInteractionList(
      const InitialState & init_state) const
{
  LOG("IntLst", pINFO) << "InitialState = " << init_state.AsString();

  int nupdg  = init_state.ProbePdg();
  int tgtpdg = init_state.Tgt().Pdg();

  const Target & target = init_state.Tgt();

  if(target.A() < 4) return 0;

  InteractionList * intlist = new InteractionList;

  if(!fSetDiNucleonCode&&fIsCC) {
	 LOG("IntLst", pWARN) << "fIsCC(val) = " << fIsCC;
      Interaction * interaction = Interaction::MECCC(tgtpdg, nupdg, 0.0);
      intlist->push_back(interaction);
  }

  if(!fSetDiNucleonCodeEM&&fIsEM) {
   LOG("IntLst", pWARN) << "fIsCC(val) = " << fIsCC;
      Interaction * interaction = Interaction::MECEM(tgtpdg, nupdg, 0.0);
      intlist->push_back(interaction);
  }

  const int nc = 3;
  const int nucleon_cluster[nc] = {
    kPdgClusterNN, kPdgClusterNP, kPdgClusterPP };

  for(int ic = 0; ic < nc; ic++) {
     int ncpdg = nucleon_cluster[ic];
     if(fIsCC&&fSetDiNucleonCode) {
       bool allowed = false;
       LOG("IntLst", pWARN) << "fIsCC(emp) = " << fIsCC;
       if(pdg::IsNeutrino(nupdg)) {
         // neutrino CC => final state primary lepton is -1
         // therefore the nucleon-cluster charge needs to be incremented by +1.
	 if(ncpdg == kPdgClusterNN || ncpdg == kPdgClusterNP) {
	   allowed = true;
	 }
       }
       else
	 if(pdg::IsAntiNeutrino(nupdg)) {
	   // anti-neutrino CC => final state primary lepton is +1
	   // therefore the nucleon-cluster charge needs to be incremented by -1.
	   if(ncpdg == kPdgClusterNP || ncpdg == kPdgClusterPP) {
	     allowed = true;
	   }
	 }
       if(allowed) {
	 Interaction * interaction =
	   Interaction::MECCC(tgtpdg,ncpdg,nupdg,0);
	 intlist->push_back(interaction);
       }
     }//CC?
     else
       if(fIsNC)
     {
       LOG("IntLst", pWARN) << "fIsNC = " << fIsNC;
       Interaction * interaction =
	 Interaction::MECNC(tgtpdg,ncpdg,nupdg,0);
       intlist->push_back(interaction);
     }//NC?
     else
       if(fIsEM&&fSetDiNucleonCodeEM) {
	 LOG("IntLst", pWARN) << "fIsEM = " << fIsEM << "  ncpdg = " << ncpdg;
	 Interaction * interaction =
	   Interaction::MECEM(tgtpdg,ncpdg,nupdg,0);
	 intlist->push_back(interaction);
       }//EM?
  }
  return intlist;

}
//___________________________________________________________________________
void MECInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void MECInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void MECInteractionListGenerator::LoadConfigData(void)
{
  GetParamDef( "is-CC", fIsCC, false ) ;
  GetParamDef( "is-NC", fIsNC, false ) ;
  GetParamDef( "is-EM", fIsEM, false ) ;

  GetParam( "SetDiNucleonCode", fSetDiNucleonCode ) ;
  GetParam( "SetDiNucleonCodeEM", fSetDiNucleonCodeEM ) ;

}
//____________________________________________________________________________
