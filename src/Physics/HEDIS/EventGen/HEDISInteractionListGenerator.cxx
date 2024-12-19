//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/EventGen/HEDISInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
HEDISInteractionListGenerator::HEDISInteractionListGenerator() :
InteractionListGeneratorI("genie::HEDISInteractionListGenerator")
{

}
//___________________________________________________________________________
HEDISInteractionListGenerator::HEDISInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::HEDISInteractionListGenerator", config)
{

}
//___________________________________________________________________________
HEDISInteractionListGenerator::~HEDISInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * HEDISInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  vector<InteractionType_t> inttype;
  if      (fIsCC) inttype.push_back(kIntWeakCC);
  else if (fIsNC) inttype.push_back(kIntWeakNC);
  else {
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  int ppdg = init_state.ProbePdg();
  if( !pdg::IsLepton(ppdg) ) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  vector<InitialState> init;
  init.push_back(init_state);
  InteractionList * intlist = this->CreateHEDISlist(init,inttype);

  if(intlist->size() == 0) {
     LOG("IntLst", pERROR)
         << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete intlist;
     return 0;
  }
  return intlist;

}
//____________________________________________________________________________
InteractionList * HEDISInteractionListGenerator::CreateHEDISlist(
    vector<InitialState> vinit, vector<InteractionType_t> vinttype) const
{

  InteractionList * intlist = new InteractionList;

  vector<InitialState>::const_iterator init = vinit.begin();
  for( ; init != vinit.end(); ++init) {

    vector<int> nucl;
    if (init->Tgt().Z()>0)                nucl.push_back(kPdgProton);
    if (init->Tgt().A()-init->Tgt().Z()>0) nucl.push_back(kPdgNeutron);

    vector<int>::const_iterator inucl = nucl.begin();
    for( ; inucl != nucl.end(); ++inucl) {

      vector<InteractionType_t>::const_iterator inttype = vinttype.begin();
      for( ; inttype != vinttype.end(); ++inttype) {
      
        ProcessInfo proc(kScDeepInelastic,*inttype);
        Interaction * interaction = new Interaction(*init, proc);
        interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(*inucl);
        multimap<int,bool> hq = this->GetHitQuarks(interaction);
        multimap<int,bool>::const_iterator hqi = hq.begin();
        for( ; hqi != hq.end(); ++hqi) {
          int  quark_code = hqi->first;
          bool from_sea   = hqi->second;
          interaction->InitStatePtr()->TgtPtr()->SetHitQrkPdg(quark_code);
          interaction->InitStatePtr()->TgtPtr()->SetHitSeaQrk(from_sea);
          vector<int> fq = this->GetFinalQuarks(interaction);
          vector<int>::const_iterator fqi = fq.begin();
          for( ; fqi != fq.end(); ++fqi) {
            XclsTag exclusive_tag;
            exclusive_tag.SetFinalQuark (*fqi);
            interaction->SetExclTag(exclusive_tag);
            Interaction * intq = new Interaction(*interaction);
            intlist->push_back(intq);
          }   
        }
        delete interaction;
      }
    }
  }

  return intlist;

}
//____________________________________________________________________________
multimap<int,bool> HEDISInteractionListGenerator::GetHitQuarks(
                                        const Interaction * interaction) const
{
// Set (PDG code, from-sea flag) for all possible hit and final quarks for the 
//  input interaction

  multimap<int,bool> hq;

  const ProcessInfo & proc = interaction->ProcInfo();

  if(proc.IsWeakNC()) {
    //
    // NC - includes both v+N, vbar+N
    //
    hq.insert(pair<int,bool>(kPdgDQuark,     false));
    hq.insert(pair<int,bool>(kPdgUQuark,     false));
    hq.insert(pair<int,bool>(kPdgDQuark,     true ));
    hq.insert(pair<int,bool>(kPdgUQuark,     true ));
    hq.insert(pair<int,bool>(kPdgSQuark,     true ));
    hq.insert(pair<int,bool>(kPdgCQuark,     true ));
    hq.insert(pair<int,bool>(kPdgBQuark,     true ));
    hq.insert(pair<int,bool>(kPdgAntiDQuark, true ));
    hq.insert(pair<int,bool>(kPdgAntiUQuark, true ));
    hq.insert(pair<int,bool>(kPdgAntiSQuark, true ));
    hq.insert(pair<int,bool>(kPdgAntiCQuark, true ));
    hq.insert(pair<int,bool>(kPdgAntiBQuark, true ));

  } else if (proc.IsWeakCC()) {
    //
    // CC - only I=-1/2 quarks for v+N & I=1/2 quarks for vbar+N
    //
    int ppdg = interaction->InitState().ProbePdg();

    if (pdg::IsNeutrino(ppdg)){
      hq.insert(pair<int,bool>(kPdgDQuark,     false));
      hq.insert(pair<int,bool>(kPdgDQuark,     true ));
      hq.insert(pair<int,bool>(kPdgSQuark,     true ));
      hq.insert(pair<int,bool>(kPdgBQuark,     true ));
      hq.insert(pair<int,bool>(kPdgAntiUQuark, true ));
      hq.insert(pair<int,bool>(kPdgAntiCQuark, true ));
    }
    else if (pdg::IsAntiNeutrino(ppdg)){
      hq.insert(pair<int,bool>(kPdgUQuark,     false));
      hq.insert(pair<int,bool>(kPdgUQuark,     true ));
      hq.insert(pair<int,bool>(kPdgCQuark,     true ));
      hq.insert(pair<int,bool>(kPdgAntiDQuark, true ));
      hq.insert(pair<int,bool>(kPdgAntiSQuark, true ));
      hq.insert(pair<int,bool>(kPdgAntiBQuark, true ));
    }
  }//CC or NC

  return hq;
}
//____________________________________________________________________________
vector<int> HEDISInteractionListGenerator::GetFinalQuarks(
                                 const Interaction * interaction) const
{

  const ProcessInfo & proc = interaction->ProcInfo();
  int hitquark = interaction->InitState().Tgt().HitQrkPdg();

  vector<int> fq;
  if(proc.IsWeakCC()) {
    int abshq = TMath::Abs(hitquark);
    int sign = hitquark/abshq;
    if (abshq==kPdgUQuark || abshq==kPdgCQuark) {
      fq.push_back(sign*kPdgDQuark); 
      fq.push_back(sign*kPdgSQuark); 
      fq.push_back(sign*kPdgBQuark); 
    }
    else if (abshq==kPdgDQuark || abshq==kPdgSQuark || abshq==kPdgBQuark) {
      fq.push_back(sign*kPdgUQuark); 
      fq.push_back(sign*kPdgCQuark); 
      fq.push_back(sign*kPdgTQuark); 
    }
  }
  else if (proc.IsWeakNC()) {
    fq.push_back(hitquark);
  }

  return fq;
}
//___________________________________________________________________________
void HEDISInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void HEDISInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void HEDISInteractionListGenerator::LoadConfigData(void)
{

  GetParamDef("is-CC", fIsCC, false ) ;
  GetParamDef("is-NC", fIsNC, false ) ;

}


