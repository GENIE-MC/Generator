//____________________________________________________________________________
/*

 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

#include "Physics/BoostedDarkMatter/EventGen/DMDISInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
DMDISInteractionListGenerator::DMDISInteractionListGenerator() :
InteractionListGeneratorI("genie::DMDISInteractionListGenerator")
{

}
//___________________________________________________________________________
DMDISInteractionListGenerator::DMDISInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::DMDISInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DMDISInteractionListGenerator::~DMDISInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DMDISInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  InteractionType_t inttype;
  // Only accept DM scattering here
  if (fIsDM) inttype = kIntDarkMatter;
  else {
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  int ppdg = init_state.ProbePdg();
  if( !pdg::IsDarkMatter(ppdg)) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  bool hasP = (init_state.Tgt().Z() > 0);
  bool hasN = (init_state.Tgt().N() > 0);

  InteractionList * intlist = new InteractionList;

  int nuclpdg[2] = { kPdgProton, kPdgNeutron };

  for(int inucl=0; inucl<2; inucl++) {

    int struck_nucleon = nuclpdg[inucl];

    if( (struck_nucleon == kPdgProton  && hasP) ||
        (struck_nucleon == kPdgNeutron && hasN) ) {

      ProcessInfo proc_info(kScDarkMatterDeepInelastic, inttype);

      Interaction * interaction = new Interaction(init_state, proc_info);
      Target * target = interaction->InitStatePtr()->TgtPtr();
      target->SetHitNucPdg(struck_nucleon);

      if(fIsCharm) {
         XclsTag exclusive_tag;
         exclusive_tag.SetCharm();
         interaction->SetExclTag(exclusive_tag);
      }

      if(fSetHitQuark) {
        // Add interactions for all possible hit (valence or sea) quarks

        multimap<int,bool> hq = this->GetHitQuarks(interaction);
        multimap<int,bool>::const_iterator hqi = hq.begin();

        for( ; hqi != hq.end(); ++hqi) {

          int  quark_code = hqi->first;
          bool from_sea   = hqi->second;

          target->SetHitQrkPdg(quark_code);
          target->SetHitSeaQrk(from_sea);

          Interaction * intq = new Interaction(*interaction);
          intlist->push_back(intq);  
        }
        delete interaction; 
      }
      else {
        intlist->push_back(interaction);
      }//set hit quark?

    }//current N exists in nuclear target
  }//N

  if(intlist->size() == 0) {
     LOG("IntLst", pERROR)
         << "Returning NULL InteractionList for init-state: "
                                                  << init_state.AsString();
     delete intlist;
     return 0;
  }
  return intlist;
}
//___________________________________________________________________________
void DMDISInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DMDISInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DMDISInteractionListGenerator::LoadConfigData(void)
{
        GetParamDef( "is-DM", fIsDM, false ) ;
	GetParamDef( "is-Charm", fIsCharm, false ) ;
	GetParamDef( "set-hit-quark", fSetHitQuark, false ) ;
        
}
//____________________________________________________________________________
multimap<int,bool> DMDISInteractionListGenerator::GetHitQuarks(
                                        const Interaction * interaction) const
{
// Set (PDG code, from-sea flag) for all possible hit quarks for the input
// interaction

  multimap<int,bool> hq;

  const ProcessInfo & proc = interaction->ProcInfo();

  if(!fIsCharm) {
    hq.insert(pair<int,bool>(kPdgUQuark,     false));
    hq.insert(pair<int,bool>(kPdgUQuark,     true ));
    hq.insert(pair<int,bool>(kPdgAntiUQuark, true ));
    hq.insert(pair<int,bool>(kPdgDQuark,     false));
    hq.insert(pair<int,bool>(kPdgDQuark,     true ));
    hq.insert(pair<int,bool>(kPdgAntiDQuark, true ));
    hq.insert(pair<int,bool>(kPdgSQuark,     true ));
    hq.insert(pair<int,bool>(kPdgAntiSQuark, true ));
  }
  
  return hq;
}
//____________________________________________________________________________

