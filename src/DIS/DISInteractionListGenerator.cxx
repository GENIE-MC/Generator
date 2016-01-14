//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new DIS package from its previous location (EVGModules).
 @ Sep 21, 2009 - CA
   Generate interaction lists for charge lepton scattering.

*/
//____________________________________________________________________________

#include "DIS/DISInteractionListGenerator.h"
#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
DISInteractionListGenerator::DISInteractionListGenerator() :
InteractionListGeneratorI("genie::DISInteractionListGenerator")
{

}
//___________________________________________________________________________
DISInteractionListGenerator::DISInteractionListGenerator(string config) :
InteractionListGeneratorI("genie::DISInteractionListGenerator", config)
{

}
//___________________________________________________________________________
DISInteractionListGenerator::~DISInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * DISInteractionListGenerator::CreateInteractionList(
                                      const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();

  InteractionType_t inttype;
  if      (fIsCC) inttype = kIntWeakCC;
  else if (fIsNC) inttype = kIntWeakNC;
  else if (fIsEM) inttype = kIntEM;
  else {
     LOG("IntLst", pWARN)
       << "Unknown InteractionType! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  int      ppdg   = init_state.ProbePdg();
  Target * target = init_state.TgtPtr();

  if( !pdg::IsLepton(ppdg) ) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }

  bool hasP = (target->Z() > 0);
  bool hasN = (target->N() > 0);

  InteractionList * intlist = new InteractionList;

  int nuclpdg[2] = { kPdgProton, kPdgNeutron };

  for(int inucl=0; inucl<2; inucl++) {

    int struck_nucleon = nuclpdg[inucl];

    if( (struck_nucleon == kPdgProton  && hasP) ||
        (struck_nucleon == kPdgNeutron && hasN) ) {

      ProcessInfo proc_info(kScDeepInelastic, inttype);

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
void DISInteractionListGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DISInteractionListGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void DISInteractionListGenerator::LoadConfigData(void)
{
  fIsCC        = fConfig->GetBoolDef("is-CC",         false);
  fIsNC        = fConfig->GetBoolDef("is-NC",         false);
  fIsEM        = fConfig->GetBoolDef("is-EM",         false);
  fIsCharm     = fConfig->GetBoolDef("is-Charm",      false);
  fSetHitQuark = fConfig->GetBoolDef("set-hit-quark", false);
}
//____________________________________________________________________________
multimap<int,bool> DISInteractionListGenerator::GetHitQuarks(
                                        const Interaction * interaction) const
{
// Set (PDG code, from-sea flag) for all possible hit quarks for the input
// interaction

  multimap<int,bool> hq;

  const ProcessInfo & proc = interaction->ProcInfo();

  if(proc.IsWeakNC() || proc.IsEM()) {
    //
    // NC - includes both v+N, vbar+N
    //
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

  } else if (proc.IsWeakCC()) {
    //
    // CC - only I=-1/2 quarks for v+N & I=1/2 quarks for vbar+N
    //
    int ppdg = interaction->InitState().ProbePdg();

    if (pdg::IsNeutrino(ppdg)){
       if(!fIsCharm) hq.insert(pair<int,bool>(kPdgAntiUQuark, true ));
                     hq.insert(pair<int,bool>(kPdgDQuark,     false));
                     hq.insert(pair<int,bool>(kPdgDQuark,     true ));
                     hq.insert(pair<int,bool>(kPdgSQuark,     true ));
    } 
    else if (pdg::IsAntiNeutrino(ppdg)){
       if(!fIsCharm) hq.insert(pair<int,bool>(kPdgUQuark,     false));
       if(!fIsCharm) hq.insert(pair<int,bool>(kPdgUQuark,     true ));
                     hq.insert(pair<int,bool>(kPdgAntiDQuark, true ));
                     hq.insert(pair<int,bool>(kPdgAntiSQuark, true ));
    }
  }//CC or NC

  return hq;
}
//____________________________________________________________________________

