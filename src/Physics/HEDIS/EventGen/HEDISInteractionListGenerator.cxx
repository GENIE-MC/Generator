//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/EventGen/HEDISInteractionListGenerator.h"
#include "Framework/EventGen/InteractionList.h"
#include "Framework/Interaction/Interaction.h"
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

  int nupdg = init_state.ProbePdg();
  if( !pdg::IsLepton(nupdg) ) {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
       << "for init-state: " << init_state.AsString();
     return 0;
  }


  const int n_nucc_channels = 36;
  const int n_nunc_channels = 24;
  HEDISQrkChannel_t nucc_channels[n_nucc_channels] = {kHEDISQrkNull};
  HEDISQrkChannel_t nunc_channels[n_nunc_channels] = {kHEDISQrkNull};

  if( pdg::IsNeutrino(nupdg) ) {    
    nucc_channels[0]  = kHEDISQrk_v_cc_p_dval_u,
    nucc_channels[1]  = kHEDISQrk_v_cc_p_dval_c,
    nucc_channels[2]  = kHEDISQrk_v_cc_p_dval_t;
    nucc_channels[3]  = kHEDISQrk_v_cc_p_dsea_u;
    nucc_channels[4]  = kHEDISQrk_v_cc_p_dsea_c;
    nucc_channels[5]  = kHEDISQrk_v_cc_p_dsea_t;
    nucc_channels[6]  = kHEDISQrk_v_cc_p_ssea_u;
    nucc_channels[7]  = kHEDISQrk_v_cc_p_ssea_c;
    nucc_channels[8]  = kHEDISQrk_v_cc_p_ssea_t;
    nucc_channels[9]  = kHEDISQrk_v_cc_p_bsea_u;
    nucc_channels[10] = kHEDISQrk_v_cc_p_bsea_c;
    nucc_channels[11] = kHEDISQrk_v_cc_p_bsea_t;
    nucc_channels[12] = kHEDISQrk_v_cc_p_ubarsea_dbar;
    nucc_channels[13] = kHEDISQrk_v_cc_p_ubarsea_sbar;
    nucc_channels[14] = kHEDISQrk_v_cc_p_ubarsea_bbar;
    nucc_channels[15] = kHEDISQrk_v_cc_p_cbarsea_dbar;
    nucc_channels[16] = kHEDISQrk_v_cc_p_cbarsea_sbar;
    nucc_channels[17] = kHEDISQrk_v_cc_p_cbarsea_bbar;
    nucc_channels[18] = kHEDISQrk_v_cc_n_dval_u;
    nucc_channels[19] = kHEDISQrk_v_cc_n_dval_c;
    nucc_channels[20] = kHEDISQrk_v_cc_n_dval_t;
    nucc_channels[21] = kHEDISQrk_v_cc_n_dsea_u;
    nucc_channels[22] = kHEDISQrk_v_cc_n_dsea_c;
    nucc_channels[23] = kHEDISQrk_v_cc_n_dsea_t;
    nucc_channels[24] = kHEDISQrk_v_cc_n_ssea_u;
    nucc_channels[25] = kHEDISQrk_v_cc_n_ssea_c;
    nucc_channels[26] = kHEDISQrk_v_cc_n_ssea_t;
    nucc_channels[27] = kHEDISQrk_v_cc_n_bsea_u;
    nucc_channels[28] = kHEDISQrk_v_cc_n_bsea_c;
    nucc_channels[29] = kHEDISQrk_v_cc_n_bsea_t;
    nucc_channels[30] = kHEDISQrk_v_cc_n_ubarsea_dbar;
    nucc_channels[31] = kHEDISQrk_v_cc_n_ubarsea_sbar;
    nucc_channels[32] = kHEDISQrk_v_cc_n_ubarsea_bbar;
    nucc_channels[33] = kHEDISQrk_v_cc_n_cbarsea_dbar;
    nucc_channels[34] = kHEDISQrk_v_cc_n_cbarsea_sbar;
    nucc_channels[35] = kHEDISQrk_v_cc_n_cbarsea_bbar;
    nunc_channels[0]  = kHEDISQrk_v_nc_p_dval_d;
    nunc_channels[1]  = kHEDISQrk_v_nc_p_uval_u;
    nunc_channels[2]  = kHEDISQrk_v_nc_p_dsea_d;
    nunc_channels[3]  = kHEDISQrk_v_nc_p_usea_u;
    nunc_channels[4]  = kHEDISQrk_v_nc_p_ssea_s;
    nunc_channels[5]  = kHEDISQrk_v_nc_p_csea_c;
    nunc_channels[6]  = kHEDISQrk_v_nc_p_bsea_b;
    nunc_channels[7]  = kHEDISQrk_v_nc_p_dbarsea_dbar;
    nunc_channels[8]  = kHEDISQrk_v_nc_p_ubarsea_ubar;
    nunc_channels[9]  = kHEDISQrk_v_nc_p_sbarsea_sbar;
    nunc_channels[10] = kHEDISQrk_v_nc_p_cbarsea_cbar;
    nunc_channels[11] = kHEDISQrk_v_nc_p_bbarsea_bbar;
    nunc_channels[12] = kHEDISQrk_v_nc_n_dval_d;
    nunc_channels[13] = kHEDISQrk_v_nc_n_uval_u;
    nunc_channels[14] = kHEDISQrk_v_nc_n_dsea_d;
    nunc_channels[15] = kHEDISQrk_v_nc_n_usea_u;
    nunc_channels[16] = kHEDISQrk_v_nc_n_ssea_s;
    nunc_channels[17] = kHEDISQrk_v_nc_n_csea_c;
    nunc_channels[18] = kHEDISQrk_v_nc_n_bsea_b;
    nunc_channels[19] = kHEDISQrk_v_nc_n_dbarsea_dbar;
    nunc_channels[20] = kHEDISQrk_v_nc_n_ubarsea_ubar;
    nunc_channels[21] = kHEDISQrk_v_nc_n_sbarsea_sbar;
    nunc_channels[22] = kHEDISQrk_v_nc_n_cbarsea_cbar;
    nunc_channels[23] = kHEDISQrk_v_nc_n_bbarsea_bbar;
  } 
  else if ( pdg::IsAntiNeutrino(nupdg) ) {
    nucc_channels[0]  = kHEDISQrk_vbar_cc_p_uval_d;
    nucc_channels[1]  = kHEDISQrk_vbar_cc_p_uval_s;
    nucc_channels[2]  = kHEDISQrk_vbar_cc_p_uval_b;
    nucc_channels[3]  = kHEDISQrk_vbar_cc_p_usea_d;
    nucc_channels[4]  = kHEDISQrk_vbar_cc_p_usea_s;
    nucc_channels[5]  = kHEDISQrk_vbar_cc_p_usea_b;
    nucc_channels[6]  = kHEDISQrk_vbar_cc_p_csea_d;
    nucc_channels[7]  = kHEDISQrk_vbar_cc_p_csea_s;
    nucc_channels[8]  = kHEDISQrk_vbar_cc_p_csea_b;
    nucc_channels[9]  = kHEDISQrk_vbar_cc_p_dbarsea_ubar;
    nucc_channels[10] = kHEDISQrk_vbar_cc_p_dbarsea_cbar;
    nucc_channels[11] = kHEDISQrk_vbar_cc_p_dbarsea_tbar;
    nucc_channels[12] = kHEDISQrk_vbar_cc_p_sbarsea_ubar;
    nucc_channels[13] = kHEDISQrk_vbar_cc_p_sbarsea_cbar;
    nucc_channels[14] = kHEDISQrk_vbar_cc_p_sbarsea_tbar;
    nucc_channels[15] = kHEDISQrk_vbar_cc_p_bbarsea_ubar;
    nucc_channels[16] = kHEDISQrk_vbar_cc_p_bbarsea_cbar;
    nucc_channels[17] = kHEDISQrk_vbar_cc_p_bbarsea_tbar;
    nucc_channels[18] = kHEDISQrk_vbar_cc_n_uval_d;
    nucc_channels[19] = kHEDISQrk_vbar_cc_n_uval_s;
    nucc_channels[20] = kHEDISQrk_vbar_cc_n_uval_b;
    nucc_channels[21] = kHEDISQrk_vbar_cc_n_usea_d;
    nucc_channels[22] = kHEDISQrk_vbar_cc_n_usea_s;
    nucc_channels[23] = kHEDISQrk_vbar_cc_n_usea_b;
    nucc_channels[24] = kHEDISQrk_vbar_cc_n_csea_d;
    nucc_channels[25] = kHEDISQrk_vbar_cc_n_csea_s;
    nucc_channels[26] = kHEDISQrk_vbar_cc_n_csea_b;
    nucc_channels[27] = kHEDISQrk_vbar_cc_n_dbarsea_ubar;
    nucc_channels[28] = kHEDISQrk_vbar_cc_n_dbarsea_cbar;
    nucc_channels[29] = kHEDISQrk_vbar_cc_n_dbarsea_tbar;
    nucc_channels[30] = kHEDISQrk_vbar_cc_n_sbarsea_ubar;
    nucc_channels[31] = kHEDISQrk_vbar_cc_n_sbarsea_cbar;
    nucc_channels[32] = kHEDISQrk_vbar_cc_n_sbarsea_tbar;
    nucc_channels[33] = kHEDISQrk_vbar_cc_n_bbarsea_ubar;
    nucc_channels[34] = kHEDISQrk_vbar_cc_n_bbarsea_cbar;
    nucc_channels[35] = kHEDISQrk_vbar_cc_n_bbarsea_tbar;
    nunc_channels[0]  = kHEDISQrk_vbar_nc_p_dval_d;
    nunc_channels[1]  = kHEDISQrk_vbar_nc_p_uval_u;
    nunc_channels[2]  = kHEDISQrk_vbar_nc_p_dsea_d;
    nunc_channels[3]  = kHEDISQrk_vbar_nc_p_usea_u;
    nunc_channels[4]  = kHEDISQrk_vbar_nc_p_ssea_s;
    nunc_channels[5]  = kHEDISQrk_vbar_nc_p_csea_c;
    nunc_channels[6]  = kHEDISQrk_vbar_nc_p_bsea_b;
    nunc_channels[7]  = kHEDISQrk_vbar_nc_p_dbarsea_dbar;
    nunc_channels[8]  = kHEDISQrk_vbar_nc_p_ubarsea_ubar;
    nunc_channels[9]  = kHEDISQrk_vbar_nc_p_sbarsea_sbar;
    nunc_channels[10] = kHEDISQrk_vbar_nc_p_cbarsea_cbar;
    nunc_channels[11] = kHEDISQrk_vbar_nc_p_bbarsea_bbar;
    nunc_channels[12] = kHEDISQrk_vbar_nc_n_dval_d;
    nunc_channels[13] = kHEDISQrk_vbar_nc_n_uval_u;
    nunc_channels[14] = kHEDISQrk_vbar_nc_n_dsea_d;
    nunc_channels[15] = kHEDISQrk_vbar_nc_n_usea_u;
    nunc_channels[16] = kHEDISQrk_vbar_nc_n_ssea_s;
    nunc_channels[17] = kHEDISQrk_vbar_nc_n_csea_c;
    nunc_channels[18] = kHEDISQrk_vbar_nc_n_bsea_b;
    nunc_channels[19] = kHEDISQrk_vbar_nc_n_dbarsea_dbar;
    nunc_channels[20] = kHEDISQrk_vbar_nc_n_ubarsea_ubar;
    nunc_channels[21] = kHEDISQrk_vbar_nc_n_sbarsea_sbar;
    nunc_channels[22] = kHEDISQrk_vbar_nc_n_cbarsea_cbar;
    nunc_channels[23] = kHEDISQrk_vbar_nc_n_bbarsea_bbar;
  } 
  else {
     LOG("IntLst", pWARN)
       << "Can not handle probe! Returning NULL InteractionList "
                         << "for init-state: " << init_state.AsString();
     return 0;
  }

  bool hasP = (init_state.Tgt().Z() > 0);
  bool hasN = (init_state.Tgt().N() > 0);

  InteractionList * intlist = new InteractionList;

  if (fIsCC) {
    for(int i=0; i<n_nucc_channels; i++) {
      int struck_nucleon = HEDISChannel::HitNuclPdg(nucc_channels[i]);
      if( (struck_nucleon == kPdgProton && hasP) || (struck_nucleon == kPdgNeutron && hasN) ) {
        ProcessInfo proc_info(kScHEDIS, kIntWeakCC);
        Interaction * interaction = new Interaction(init_state, proc_info);
        Target * target = interaction->InitStatePtr()->TgtPtr();
        target->SetHitNucPdg(struck_nucleon);
        this->AddFinalStateInfo(interaction, nucc_channels[i]);
        intlist->push_back(interaction);
      }
    }
  }
  else if (fIsNC) {
    // NC
    for(int i=0; i<n_nunc_channels; i++) {
      int struck_nucleon = HEDISChannel::HitNuclPdg(nunc_channels[i]);
      if( (struck_nucleon == kPdgProton && hasP) || (struck_nucleon == kPdgNeutron && hasN) ) {
        ProcessInfo proc_info(kScHEDIS, kIntWeakNC);
        Interaction * interaction = new Interaction(init_state, proc_info);
        interaction->InitStatePtr()->TgtPtr()->SetHitNucPdg(struck_nucleon);
        this->AddFinalStateInfo(interaction, nunc_channels[i]);
        intlist->push_back(interaction);
      }
    }//nc channels   
  }

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
void HEDISInteractionListGenerator::AddFinalStateInfo(
                       Interaction * interaction, HEDISQrkChannel_t hedischan) const
{

  bool iq_sea = HEDISChannel::HitQuarkSea(hedischan);
  int  iq_pdg = HEDISChannel::HitQuarkPdg(hedischan);
  int  fq_pdg = HEDISChannel::FnlQuarkPdg(hedischan);

  interaction->InitStatePtr()->TgtPtr()->SetHitSeaQrk(iq_sea);
  interaction->InitStatePtr()->TgtPtr()->SetHitQrkPdg(iq_pdg);

  XclsTag exclusive_tag;
  exclusive_tag.SetFinalQuark (fq_pdg);
  exclusive_tag.SetHEDISQrkChannel (hedischan);
  interaction->SetExclTag(exclusive_tag);

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
