#ifndef _HEDIS_CHANNEL_H_
#define _HEDIS_CHANNEL_H_

#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Interaction/InteractionType.h"

namespace genie {

  typedef enum EHEDISNucChannel {

    kHEDISNucNull = 0,

    kHEDISNuc_v_cc_p,
    kHEDISNuc_v_cc_n,
    kHEDISNuc_vbar_cc_p,
    kHEDISNuc_vbar_cc_n,
    kHEDISNuc_v_nc_p,
    kHEDISNuc_v_nc_n,
    kHEDISNuc_vbar_nc_p,
    kHEDISNuc_vbar_nc_n,
    kHEDISNuc_numofchannels

  } HEDISNucChannel_t;

  typedef enum EHEDISQrkChannel {

    kHEDISQrkNull = 0,

    kHEDISQrk_v_cc_p_dval_u,
    kHEDISQrk_v_cc_p_dval_c,
    kHEDISQrk_v_cc_p_dval_t,
    kHEDISQrk_v_cc_p_dsea_u,
    kHEDISQrk_v_cc_p_dsea_c,
    kHEDISQrk_v_cc_p_dsea_t,
    kHEDISQrk_v_cc_p_ssea_u,
    kHEDISQrk_v_cc_p_ssea_c,
    kHEDISQrk_v_cc_p_ssea_t,
    kHEDISQrk_v_cc_p_bsea_u,
    kHEDISQrk_v_cc_p_bsea_c,
    kHEDISQrk_v_cc_p_bsea_t,
    kHEDISQrk_v_cc_p_ubarsea_dbar,
    kHEDISQrk_v_cc_p_ubarsea_sbar,
    kHEDISQrk_v_cc_p_ubarsea_bbar,
    kHEDISQrk_v_cc_p_cbarsea_dbar,
    kHEDISQrk_v_cc_p_cbarsea_sbar,
    kHEDISQrk_v_cc_p_cbarsea_bbar,

    kHEDISQrk_v_cc_n_dval_u,
    kHEDISQrk_v_cc_n_dval_c,
    kHEDISQrk_v_cc_n_dval_t,
    kHEDISQrk_v_cc_n_dsea_u,
    kHEDISQrk_v_cc_n_dsea_c,
    kHEDISQrk_v_cc_n_dsea_t,
    kHEDISQrk_v_cc_n_ssea_u,
    kHEDISQrk_v_cc_n_ssea_c,
    kHEDISQrk_v_cc_n_ssea_t,
    kHEDISQrk_v_cc_n_bsea_u,
    kHEDISQrk_v_cc_n_bsea_c,
    kHEDISQrk_v_cc_n_bsea_t,
    kHEDISQrk_v_cc_n_ubarsea_dbar,
    kHEDISQrk_v_cc_n_ubarsea_sbar,
    kHEDISQrk_v_cc_n_ubarsea_bbar,
    kHEDISQrk_v_cc_n_cbarsea_dbar,
    kHEDISQrk_v_cc_n_cbarsea_sbar,
    kHEDISQrk_v_cc_n_cbarsea_bbar,

    kHEDISQrk_vbar_cc_p_uval_d,
    kHEDISQrk_vbar_cc_p_uval_s,
    kHEDISQrk_vbar_cc_p_uval_b,
    kHEDISQrk_vbar_cc_p_usea_d,
    kHEDISQrk_vbar_cc_p_usea_s,
    kHEDISQrk_vbar_cc_p_usea_b,
    kHEDISQrk_vbar_cc_p_csea_d,
    kHEDISQrk_vbar_cc_p_csea_s,
    kHEDISQrk_vbar_cc_p_csea_b,
    kHEDISQrk_vbar_cc_p_dbarsea_ubar,
    kHEDISQrk_vbar_cc_p_dbarsea_cbar,
    kHEDISQrk_vbar_cc_p_dbarsea_tbar,
    kHEDISQrk_vbar_cc_p_sbarsea_ubar,
    kHEDISQrk_vbar_cc_p_sbarsea_cbar,
    kHEDISQrk_vbar_cc_p_sbarsea_tbar,
    kHEDISQrk_vbar_cc_p_bbarsea_ubar,
    kHEDISQrk_vbar_cc_p_bbarsea_cbar,
    kHEDISQrk_vbar_cc_p_bbarsea_tbar,

    kHEDISQrk_vbar_cc_n_uval_d,
    kHEDISQrk_vbar_cc_n_uval_s,
    kHEDISQrk_vbar_cc_n_uval_b,
    kHEDISQrk_vbar_cc_n_usea_d,
    kHEDISQrk_vbar_cc_n_usea_s,
    kHEDISQrk_vbar_cc_n_usea_b,
    kHEDISQrk_vbar_cc_n_csea_d,
    kHEDISQrk_vbar_cc_n_csea_s,
    kHEDISQrk_vbar_cc_n_csea_b,
    kHEDISQrk_vbar_cc_n_dbarsea_ubar,
    kHEDISQrk_vbar_cc_n_dbarsea_cbar,
    kHEDISQrk_vbar_cc_n_dbarsea_tbar,
    kHEDISQrk_vbar_cc_n_sbarsea_ubar,
    kHEDISQrk_vbar_cc_n_sbarsea_cbar,
    kHEDISQrk_vbar_cc_n_sbarsea_tbar,
    kHEDISQrk_vbar_cc_n_bbarsea_ubar,
    kHEDISQrk_vbar_cc_n_bbarsea_cbar,
    kHEDISQrk_vbar_cc_n_bbarsea_tbar,

    kHEDISQrk_v_nc_p_dval_d,
    kHEDISQrk_v_nc_p_uval_u,
    kHEDISQrk_v_nc_p_dsea_d,
    kHEDISQrk_v_nc_p_usea_u,
    kHEDISQrk_v_nc_p_ssea_s,
    kHEDISQrk_v_nc_p_csea_c,
    kHEDISQrk_v_nc_p_bsea_b,
    kHEDISQrk_v_nc_p_dbarsea_dbar,
    kHEDISQrk_v_nc_p_ubarsea_ubar,
    kHEDISQrk_v_nc_p_sbarsea_sbar,
    kHEDISQrk_v_nc_p_cbarsea_cbar,
    kHEDISQrk_v_nc_p_bbarsea_bbar,

    kHEDISQrk_v_nc_n_dval_d,
    kHEDISQrk_v_nc_n_uval_u,
    kHEDISQrk_v_nc_n_dsea_d,
    kHEDISQrk_v_nc_n_usea_u,
    kHEDISQrk_v_nc_n_ssea_s,
    kHEDISQrk_v_nc_n_csea_c,
    kHEDISQrk_v_nc_n_bsea_b,
    kHEDISQrk_v_nc_n_dbarsea_dbar,
    kHEDISQrk_v_nc_n_ubarsea_ubar,
    kHEDISQrk_v_nc_n_sbarsea_sbar,
    kHEDISQrk_v_nc_n_cbarsea_cbar,
    kHEDISQrk_v_nc_n_bbarsea_bbar,

    kHEDISQrk_vbar_nc_p_dval_d,
    kHEDISQrk_vbar_nc_p_uval_u,
    kHEDISQrk_vbar_nc_p_dsea_d,
    kHEDISQrk_vbar_nc_p_usea_u,
    kHEDISQrk_vbar_nc_p_ssea_s,
    kHEDISQrk_vbar_nc_p_csea_c,
    kHEDISQrk_vbar_nc_p_bsea_b,
    kHEDISQrk_vbar_nc_p_dbarsea_dbar,
    kHEDISQrk_vbar_nc_p_ubarsea_ubar,
    kHEDISQrk_vbar_nc_p_sbarsea_sbar,
    kHEDISQrk_vbar_nc_p_cbarsea_cbar,
    kHEDISQrk_vbar_nc_p_bbarsea_bbar,

    kHEDISQrk_vbar_nc_n_dval_d,
    kHEDISQrk_vbar_nc_n_uval_u,
    kHEDISQrk_vbar_nc_n_dsea_d,
    kHEDISQrk_vbar_nc_n_usea_u,
    kHEDISQrk_vbar_nc_n_ssea_s,
    kHEDISQrk_vbar_nc_n_csea_c,
    kHEDISQrk_vbar_nc_n_bsea_b,
    kHEDISQrk_vbar_nc_n_dbarsea_dbar,
    kHEDISQrk_vbar_nc_n_ubarsea_ubar,
    kHEDISQrk_vbar_nc_n_sbarsea_sbar,
    kHEDISQrk_vbar_nc_n_cbarsea_cbar,
    kHEDISQrk_vbar_nc_n_bbarsea_bbar,
    kHEDISQrk_numofchannels

  } HEDISQrkChannel_t;


  class HEDISChannel
  {
    public:


      //__________________________________________________________________________
      static string AsString(HEDISNucChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNuc_v_cc_p)       : return "nu_CC_p";      break;
          case (kHEDISNuc_v_cc_n)       : return "nu_CC_n";      break;
          case (kHEDISNuc_vbar_cc_p)    : return "nubar_CC_p";   break;
          case (kHEDISNuc_vbar_cc_n)    : return "nubar_CC_n";   break;
          case (kHEDISNuc_v_nc_p)       : return "nu_NC_p";      break;
          case (kHEDISNuc_v_nc_n)       : return "nu_NC_n";      break;
          case (kHEDISNuc_vbar_nc_p)    : return "nubar_NC_p";   break;
          case (kHEDISNuc_vbar_nc_n)    : return "nubar_NC_n";   break;

          default : return "Unknown";  break;
        }
        return "Unknown";
      }
      //__________________________________________________________________________
      static HEDISQrkChannel_t GetFirstHEDISQrkChannel(HEDISNucChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNuc_v_cc_p)       : return kHEDISQrk_v_cc_p_dval_u;      break;
          case (kHEDISNuc_v_cc_n)       : return kHEDISQrk_v_cc_n_dval_u;      break;
          case (kHEDISNuc_vbar_cc_p)    : return kHEDISQrk_vbar_cc_p_uval_d;   break;
          case (kHEDISNuc_vbar_cc_n)    : return kHEDISQrk_vbar_cc_n_uval_d;   break;
          case (kHEDISNuc_v_nc_p)       : return kHEDISQrk_v_nc_p_dval_d;      break;
          case (kHEDISNuc_v_nc_n)       : return kHEDISQrk_v_nc_n_dval_d;      break;
          case (kHEDISNuc_vbar_nc_p)    : return kHEDISQrk_vbar_nc_p_dval_d;   break;
          case (kHEDISNuc_vbar_nc_n)    : return kHEDISQrk_vbar_nc_n_dval_d;   break;

          default : return kHEDISQrkNull;  break;
        }
        return kHEDISQrkNull;
      }
      //__________________________________________________________________________
      static HEDISQrkChannel_t GetLastHEDISQrkChannel(HEDISNucChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNuc_v_cc_p)       : return kHEDISQrk_v_cc_p_cbarsea_bbar;      break;
          case (kHEDISNuc_v_cc_n)       : return kHEDISQrk_v_cc_n_cbarsea_bbar;      break;
          case (kHEDISNuc_vbar_cc_p)    : return kHEDISQrk_vbar_cc_p_bbarsea_tbar;   break;
          case (kHEDISNuc_vbar_cc_n)    : return kHEDISQrk_vbar_cc_n_bbarsea_tbar;   break;
          case (kHEDISNuc_v_nc_p)       : return kHEDISQrk_v_nc_p_bbarsea_bbar;      break;
          case (kHEDISNuc_v_nc_n)       : return kHEDISQrk_v_nc_n_bbarsea_bbar;      break;
          case (kHEDISNuc_vbar_nc_p)    : return kHEDISQrk_vbar_nc_p_bbarsea_bbar;   break;
          case (kHEDISNuc_vbar_nc_n)    : return kHEDISQrk_vbar_nc_n_bbarsea_bbar;   break;

          default : return kHEDISQrkNull;  break;
        }
        return kHEDISQrkNull;
      }
      static InteractionType_t InteractionType(HEDISNucChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNuc_v_nc_p)       : return kIntWeakNC;   break;
          case (kHEDISNuc_v_nc_n)       : return kIntWeakNC;   break;
          case (kHEDISNuc_vbar_nc_p)    : return kIntWeakNC;   break;
          case (kHEDISNuc_vbar_nc_n)    : return kIntWeakNC;   break;

          default : return kIntWeakCC;  break;
        }
        return kIntNull;
      }
      //__________________________________________________________________________
      static bool IsNu(HEDISNucChannel_t channel)
      {
        switch (channel) {
          case (kHEDISNuc_vbar_cc_p)    : return false;   break;
          case (kHEDISNuc_vbar_cc_n)    : return false;   break;
          case (kHEDISNuc_vbar_nc_p)    : return false;   break;
          case (kHEDISNuc_vbar_nc_n)    : return false;   break;

          default : return true;  break;
        }
        return true;
      }
      //__________________________________________________________________________
      static int HitNuclPdg(HEDISNucChannel_t channel)
      {

        switch (channel) {
          case (kHEDISNuc_v_cc_n)       : return kPdgNeutron;      break;
          case (kHEDISNuc_vbar_cc_n)    : return kPdgNeutron;   break;
          case (kHEDISNuc_v_nc_n)       : return kPdgNeutron;      break;
          case (kHEDISNuc_vbar_nc_n)    : return kPdgNeutron;   break;

          default : return kPdgProton;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      static string AsString(HEDISQrkChannel_t channel)
      {

        switch (channel) {
          case (kHEDISQrk_v_cc_p_dval_u)       : return "nu_CC_p_dval_u";   break;
          case (kHEDISQrk_v_cc_p_dval_c)       : return "nu_CC_p_dval_c";   break;
          case (kHEDISQrk_v_cc_p_dval_t)       : return "nu_CC_p_dval_t";   break;
          case (kHEDISQrk_v_cc_p_dsea_u)       : return "nu_CC_p_dsea_u";   break;
          case (kHEDISQrk_v_cc_p_dsea_c)       : return "nu_CC_p_dsea_c";   break;
          case (kHEDISQrk_v_cc_p_dsea_t)       : return "nu_CC_p_dsea_t";   break;
          case (kHEDISQrk_v_cc_p_ssea_u)       : return "nu_CC_p_ssea_u";   break;
          case (kHEDISQrk_v_cc_p_ssea_c)       : return "nu_CC_p_ssea_c";   break;
          case (kHEDISQrk_v_cc_p_ssea_t)       : return "nu_CC_p_ssea_t";   break;
          case (kHEDISQrk_v_cc_p_bsea_u)       : return "nu_CC_p_bsea_u";   break;
          case (kHEDISQrk_v_cc_p_bsea_c)       : return "nu_CC_p_bsea_c";   break;
          case (kHEDISQrk_v_cc_p_bsea_t)       : return "nu_CC_p_bsea_t";   break;
          case (kHEDISQrk_v_cc_p_ubarsea_dbar) : return "nu_CC_p_ubarsea_dbar";   break;
          case (kHEDISQrk_v_cc_p_ubarsea_sbar) : return "nu_CC_p_ubarsea_sbar";   break;
          case (kHEDISQrk_v_cc_p_ubarsea_bbar) : return "nu_CC_p_ubarsea_bbar";   break;
          case (kHEDISQrk_v_cc_p_cbarsea_dbar) : return "nu_CC_p_cbarsea_dbar";   break;
          case (kHEDISQrk_v_cc_p_cbarsea_sbar) : return "nu_CC_p_cbarsea_sbar";   break;
          case (kHEDISQrk_v_cc_p_cbarsea_bbar) : return "nu_CC_p_cbarsea_bbar";   break;

          case (kHEDISQrk_v_cc_n_dval_u)       : return "nu_CC_n_dval_u";   break;
          case (kHEDISQrk_v_cc_n_dval_c)       : return "nu_CC_n_dval_c";   break;
          case (kHEDISQrk_v_cc_n_dval_t)       : return "nu_CC_n_dval_t";   break;
          case (kHEDISQrk_v_cc_n_dsea_u)       : return "nu_CC_n_dsea_u";   break;
          case (kHEDISQrk_v_cc_n_dsea_c)       : return "nu_CC_n_dsea_c";   break;
          case (kHEDISQrk_v_cc_n_dsea_t)       : return "nu_CC_n_dsea_t";   break;
          case (kHEDISQrk_v_cc_n_ssea_u)       : return "nu_CC_n_ssea_u";   break;
          case (kHEDISQrk_v_cc_n_ssea_c)       : return "nu_CC_n_ssea_c";   break;
          case (kHEDISQrk_v_cc_n_ssea_t)       : return "nu_CC_n_ssea_t";   break;
          case (kHEDISQrk_v_cc_n_bsea_u)       : return "nu_CC_n_bsea_u";   break;
          case (kHEDISQrk_v_cc_n_bsea_c)       : return "nu_CC_n_bsea_c";   break;
          case (kHEDISQrk_v_cc_n_bsea_t)       : return "nu_CC_n_bsea_t";   break;
          case (kHEDISQrk_v_cc_n_ubarsea_dbar) : return "nu_CC_n_ubarsea_dbar";   break;
          case (kHEDISQrk_v_cc_n_ubarsea_sbar) : return "nu_CC_n_ubarsea_sbar";   break;
          case (kHEDISQrk_v_cc_n_ubarsea_bbar) : return "nu_CC_n_ubarsea_bbar";   break;
          case (kHEDISQrk_v_cc_n_cbarsea_dbar) : return "nu_CC_n_cbarsea_dbar";   break;
          case (kHEDISQrk_v_cc_n_cbarsea_sbar) : return "nu_CC_n_cbarsea_sbar";   break;
          case (kHEDISQrk_v_cc_n_cbarsea_bbar) : return "nu_CC_n_cbarsea_bbar";   break;

          case (kHEDISQrk_vbar_cc_p_uval_d)       : return "nubar_CC_p_uval_d";   break;
          case (kHEDISQrk_vbar_cc_p_uval_s)       : return "nubar_CC_p_uval_s";   break;
          case (kHEDISQrk_vbar_cc_p_uval_b)       : return "nubar_CC_p_uval_b";   break;
          case (kHEDISQrk_vbar_cc_p_usea_d)       : return "nubar_CC_p_usea_d";   break;
          case (kHEDISQrk_vbar_cc_p_usea_s)       : return "nubar_CC_p_usea_s";   break;
          case (kHEDISQrk_vbar_cc_p_usea_b)       : return "nubar_CC_p_usea_b";   break;
          case (kHEDISQrk_vbar_cc_p_csea_d)       : return "nubar_CC_p_csea_d";   break;
          case (kHEDISQrk_vbar_cc_p_csea_s)       : return "nubar_CC_p_csea_s";   break;
          case (kHEDISQrk_vbar_cc_p_csea_b)       : return "nubar_CC_p_csea_b";   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_ubar) : return "nubar_CC_p_dbarsea_ubar";   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_cbar) : return "nubar_CC_p_dbarsea_cbar";   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_tbar) : return "nubar_CC_p_dbarsea_tbar";   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_ubar) : return "nubar_CC_p_sbarsea_ubar";   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_cbar) : return "nubar_CC_p_sbarsea_cbar";   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_tbar) : return "nubar_CC_p_sbarsea_tbar";   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_ubar) : return "nubar_CC_p_bbarsea_ubar";   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_cbar) : return "nubar_CC_p_bbarsea_cbar";   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_tbar) : return "nubar_CC_p_bbarsea_tbar";   break;

          case (kHEDISQrk_vbar_cc_n_uval_d)       : return "nubar_CC_n_uval_d";   break;
          case (kHEDISQrk_vbar_cc_n_uval_s)       : return "nubar_CC_n_uval_s";   break;
          case (kHEDISQrk_vbar_cc_n_uval_b)       : return "nubar_CC_n_uval_b";   break;
          case (kHEDISQrk_vbar_cc_n_usea_d)       : return "nubar_CC_n_usea_d";   break;
          case (kHEDISQrk_vbar_cc_n_usea_s)       : return "nubar_CC_n_usea_s";   break;
          case (kHEDISQrk_vbar_cc_n_usea_b)       : return "nubar_CC_n_usea_b";   break;
          case (kHEDISQrk_vbar_cc_n_csea_d)       : return "nubar_CC_n_csea_d";   break;
          case (kHEDISQrk_vbar_cc_n_csea_s)       : return "nubar_CC_n_csea_s";   break;
          case (kHEDISQrk_vbar_cc_n_csea_b)       : return "nubar_CC_n_csea_b";   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_ubar) : return "nubar_CC_n_dbarsea_ubar";   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_cbar) : return "nubar_CC_n_dbarsea_cbar";   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_tbar) : return "nubar_CC_n_dbarsea_tbar";   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_ubar) : return "nubar_CC_n_sbarsea_ubar";   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_cbar) : return "nubar_CC_n_sbarsea_cbar";   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_tbar) : return "nubar_CC_n_sbarsea_tbar";   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_ubar) : return "nubar_CC_n_bbarsea_ubar";   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_cbar) : return "nubar_CC_n_bbarsea_cbar";   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_tbar) : return "nubar_CC_n_bbarsea_tbar";   break;

          case (kHEDISQrk_v_nc_p_dval_d)       : return "nu_NC_p_dval_d";   break;
          case (kHEDISQrk_v_nc_p_uval_u)       : return "nu_NC_p_uval_u";   break;
          case (kHEDISQrk_v_nc_p_dsea_d)       : return "nu_NC_p_dsea_d";   break;
          case (kHEDISQrk_v_nc_p_usea_u)       : return "nu_NC_p_usea_u";   break;
          case (kHEDISQrk_v_nc_p_ssea_s)       : return "nu_NC_p_ssea_s";   break;
          case (kHEDISQrk_v_nc_p_csea_c)       : return "nu_NC_p_csea_c";   break;
          case (kHEDISQrk_v_nc_p_bsea_b)       : return "nu_NC_p_bsea_b";   break;
          case (kHEDISQrk_v_nc_p_dbarsea_dbar) : return "nu_NC_p_dbarsea_dbar";   break;
          case (kHEDISQrk_v_nc_p_ubarsea_ubar) : return "nu_NC_p_ubarsea_ubar";   break;
          case (kHEDISQrk_v_nc_p_sbarsea_sbar) : return "nu_NC_p_sbarsea_sbar";   break;
          case (kHEDISQrk_v_nc_p_cbarsea_cbar) : return "nu_NC_p_cbarsea_cbar";   break;
          case (kHEDISQrk_v_nc_p_bbarsea_bbar) : return "nu_NC_p_bbarsea_bbar";   break;

          case (kHEDISQrk_v_nc_n_dval_d)       : return "nu_NC_n_dval_d";   break;
          case (kHEDISQrk_v_nc_n_uval_u)       : return "nu_NC_n_uval_u";   break;
          case (kHEDISQrk_v_nc_n_dsea_d)       : return "nu_NC_n_dsea_d";   break;
          case (kHEDISQrk_v_nc_n_usea_u)       : return "nu_NC_n_usea_u";   break;
          case (kHEDISQrk_v_nc_n_ssea_s)       : return "nu_NC_n_ssea_s";   break;
          case (kHEDISQrk_v_nc_n_csea_c)       : return "nu_NC_n_csea_c";   break;
          case (kHEDISQrk_v_nc_n_bsea_b)       : return "nu_NC_n_bsea_b";   break;
          case (kHEDISQrk_v_nc_n_dbarsea_dbar) : return "nu_NC_n_dbarsea_dbar";   break;
          case (kHEDISQrk_v_nc_n_ubarsea_ubar) : return "nu_NC_n_ubarsea_ubar";   break;
          case (kHEDISQrk_v_nc_n_sbarsea_sbar) : return "nu_NC_n_sbarsea_sbar";   break;
          case (kHEDISQrk_v_nc_n_cbarsea_cbar) : return "nu_NC_n_cbarsea_cbar";   break;
          case (kHEDISQrk_v_nc_n_bbarsea_bbar) : return "nu_NC_n_bbarsea_bbar";   break;

          case (kHEDISQrk_vbar_nc_p_dval_d)       : return "nubar_NC_p_dval_d";   break;
          case (kHEDISQrk_vbar_nc_p_uval_u)       : return "nubar_NC_p_uval_u";   break;
          case (kHEDISQrk_vbar_nc_p_dsea_d)       : return "nubar_NC_p_dsea_d";   break;
          case (kHEDISQrk_vbar_nc_p_usea_u)       : return "nubar_NC_p_usea_u";   break;
          case (kHEDISQrk_vbar_nc_p_ssea_s)       : return "nubar_NC_p_ssea_s";   break;
          case (kHEDISQrk_vbar_nc_p_csea_c)       : return "nubar_NC_p_csea_c";   break;
          case (kHEDISQrk_vbar_nc_p_bsea_b)       : return "nubar_NC_p_bsea_b";   break;
          case (kHEDISQrk_vbar_nc_p_dbarsea_dbar) : return "nubar_NC_p_dbarsea_dbar";   break;
          case (kHEDISQrk_vbar_nc_p_ubarsea_ubar) : return "nubar_NC_p_ubarsea_ubar";   break;
          case (kHEDISQrk_vbar_nc_p_sbarsea_sbar) : return "nubar_NC_p_sbarsea_sbar";   break;
          case (kHEDISQrk_vbar_nc_p_cbarsea_cbar) : return "nubar_NC_p_cbarsea_cbar";   break;
          case (kHEDISQrk_vbar_nc_p_bbarsea_bbar) : return "nubar_NC_p_bbarsea_bbar";   break;

          case (kHEDISQrk_vbar_nc_n_dval_d)       : return "nubar_NC_n_dval_d";   break;
          case (kHEDISQrk_vbar_nc_n_uval_u)       : return "nubar_NC_n_uval_u";   break;
          case (kHEDISQrk_vbar_nc_n_dsea_d)       : return "nubar_NC_n_dsea_d";   break;
          case (kHEDISQrk_vbar_nc_n_usea_u)       : return "nubar_NC_n_usea_u";   break;
          case (kHEDISQrk_vbar_nc_n_ssea_s)       : return "nubar_NC_n_ssea_s";   break;
          case (kHEDISQrk_vbar_nc_n_csea_c)       : return "nubar_NC_n_csea_c";   break;
          case (kHEDISQrk_vbar_nc_n_bsea_b)       : return "nubar_NC_n_bsea_b";   break;
          case (kHEDISQrk_vbar_nc_n_dbarsea_dbar) : return "nubar_NC_n_dbarsea_dbar";   break;
          case (kHEDISQrk_vbar_nc_n_ubarsea_ubar) : return "nubar_NC_n_ubarsea_ubar";   break;
          case (kHEDISQrk_vbar_nc_n_sbarsea_sbar) : return "nubar_NC_n_sbarsea_sbar";   break;
          case (kHEDISQrk_vbar_nc_n_cbarsea_cbar) : return "nubar_NC_n_cbarsea_cbar";   break;
          case (kHEDISQrk_vbar_nc_n_bbarsea_bbar) : return "nubar_NC_n_bbarsea_bbar";   break;


          default : return "Unknown";  break;
        }
        return "Unknown";
      }
      //__________________________________________________________________________
      static InteractionType_t InteractionType(HEDISQrkChannel_t channel)
      {

        switch (channel) {
          case (kHEDISQrk_v_nc_p_dval_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_uval_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_usea_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_csea_c)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_p_bbarsea_bbar) : return kIntWeakNC;   break;

          case (kHEDISQrk_v_nc_n_dval_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_uval_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_usea_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_csea_c)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_v_nc_n_bbarsea_bbar) : return kIntWeakNC;   break;

          case (kHEDISQrk_vbar_nc_p_dval_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_uval_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_usea_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_csea_c)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_p_bbarsea_bbar) : return kIntWeakNC;   break;

          case (kHEDISQrk_vbar_nc_n_dval_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_uval_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_dsea_d)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_usea_u)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_ssea_s)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_csea_c)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_bsea_b)       : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_dbarsea_dbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_ubarsea_ubar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_sbarsea_sbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_cbarsea_cbar) : return kIntWeakNC;   break;
          case (kHEDISQrk_vbar_nc_n_bbarsea_bbar) : return kIntWeakNC;   break;

          default : return kIntWeakCC;  break;
        }
        return kIntNull;
      }
      //__________________________________________________________________________
      static bool IsNu(HEDISQrkChannel_t channel)
      {
        switch (channel) {
          case (kHEDISQrk_vbar_cc_p_uval_d)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_uval_s)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_uval_b)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_usea_d)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_usea_s)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_usea_b)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_csea_d)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_csea_s)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_csea_b)       : return false;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_tbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_tbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_tbar) : return false;   break;

          case (kHEDISQrk_vbar_cc_n_uval_d)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_uval_s)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_uval_b)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_usea_d)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_usea_s)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_usea_b)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_csea_d)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_csea_s)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_csea_b)       : return false;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_tbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_tbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_tbar) : return false;   break;

          case (kHEDISQrk_vbar_nc_p_dval_d)       : return false;   break;
          case (kHEDISQrk_vbar_nc_p_uval_u)       : return false;   break;
          case (kHEDISQrk_vbar_nc_p_dsea_d)       : return false;   break;
          case (kHEDISQrk_vbar_nc_p_usea_u)       : return false;   break;
          case (kHEDISQrk_vbar_nc_p_ssea_s)       : return false;   break;
          case (kHEDISQrk_vbar_nc_p_csea_c)       : return false;   break;
          case (kHEDISQrk_vbar_nc_p_bsea_b)       : return false;   break;
          case (kHEDISQrk_vbar_nc_p_dbarsea_dbar) : return false;   break;
          case (kHEDISQrk_vbar_nc_p_ubarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_nc_p_sbarsea_sbar) : return false;   break;
          case (kHEDISQrk_vbar_nc_p_cbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_nc_p_bbarsea_bbar) : return false;   break;

          case (kHEDISQrk_vbar_nc_n_dval_d)       : return false;   break;
          case (kHEDISQrk_vbar_nc_n_uval_u)       : return false;   break;
          case (kHEDISQrk_vbar_nc_n_dsea_d)       : return false;   break;
          case (kHEDISQrk_vbar_nc_n_usea_u)       : return false;   break;
          case (kHEDISQrk_vbar_nc_n_ssea_s)       : return false;   break;
          case (kHEDISQrk_vbar_nc_n_csea_c)       : return false;   break;
          case (kHEDISQrk_vbar_nc_n_bsea_b)       : return false;   break;
          case (kHEDISQrk_vbar_nc_n_dbarsea_dbar) : return false;   break;
          case (kHEDISQrk_vbar_nc_n_ubarsea_ubar) : return false;   break;
          case (kHEDISQrk_vbar_nc_n_sbarsea_sbar) : return false;   break;
          case (kHEDISQrk_vbar_nc_n_cbarsea_cbar) : return false;   break;
          case (kHEDISQrk_vbar_nc_n_bbarsea_bbar) : return false;   break;

          default : return true;  break;
        }
        return true;
      }

      //__________________________________________________________________________
      static int HitNuclPdg(HEDISQrkChannel_t channel)
      {

        switch (channel) {

          case (kHEDISQrk_v_cc_n_dval_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_dval_c)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_dval_t)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_dsea_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_dsea_c)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_dsea_t)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_ssea_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_ssea_c)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_ssea_t)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_bsea_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_bsea_c)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_bsea_t)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_bbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_bbar) : return kPdgNeutron;   break;

          case (kHEDISQrk_vbar_cc_n_uval_d)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_uval_s)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_uval_b)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_usea_d)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_usea_s)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_usea_b)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_csea_d)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_csea_s)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_csea_b)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_tbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_tbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_tbar) : return kPdgNeutron;   break;

          case (kHEDISQrk_v_nc_n_dval_d)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_uval_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_dsea_d)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_usea_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_ssea_s)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_csea_c)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_bsea_b)       : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_dbarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_ubarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_sbarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_cbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_v_nc_n_bbarsea_bbar) : return kPdgNeutron;   break;

          case (kHEDISQrk_vbar_nc_n_dval_d)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_uval_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_dsea_d)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_usea_u)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_ssea_s)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_csea_c)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_bsea_b)       : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_dbarsea_dbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_ubarsea_ubar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_sbarsea_sbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_cbarsea_cbar) : return kPdgNeutron;   break;
          case (kHEDISQrk_vbar_nc_n_bbarsea_bbar) : return kPdgNeutron;   break;

          default : return kPdgProton;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      static bool HitQuarkSea(HEDISQrkChannel_t channel)
      {
        switch (channel) {
          case (kHEDISQrk_v_cc_p_dval_u)       : return false;  break;
          case (kHEDISQrk_v_cc_p_dval_c)       : return false;  break;
          case (kHEDISQrk_v_cc_p_dval_t)       : return false;  break;
          case (kHEDISQrk_v_cc_n_dval_u)       : return false;  break;
          case (kHEDISQrk_v_cc_n_dval_c)       : return false;  break;
          case (kHEDISQrk_v_cc_n_dval_t)       : return false;  break;
          case (kHEDISQrk_vbar_cc_p_uval_d)    : return false;  break;
          case (kHEDISQrk_vbar_cc_p_uval_s)    : return false;  break;
          case (kHEDISQrk_vbar_cc_p_uval_b)    : return false;  break;
          case (kHEDISQrk_vbar_cc_n_uval_d)    : return false;  break;
          case (kHEDISQrk_vbar_cc_n_uval_s)    : return false;  break;
          case (kHEDISQrk_vbar_cc_n_uval_b)    : return false;  break;
          case (kHEDISQrk_v_nc_p_dval_d)       : return false;  break;
          case (kHEDISQrk_v_nc_p_uval_u)       : return false;  break;
          case (kHEDISQrk_v_nc_n_dval_d)       : return false;  break;
          case (kHEDISQrk_v_nc_n_uval_u)       : return false;  break;
          case (kHEDISQrk_vbar_nc_p_dval_d)    : return false;  break;
          case (kHEDISQrk_vbar_nc_p_uval_u)    : return false;  break;
          case (kHEDISQrk_vbar_nc_n_dval_d)    : return false;  break;
          case (kHEDISQrk_vbar_nc_n_uval_u)    : return false;  break;

          default : return true;  break;
        }
        return true;
      }
      //__________________________________________________________________________
      static int HitQuarkPdg(HEDISQrkChannel_t channel)
      {
        switch (channel) {
          case (kHEDISQrk_v_cc_p_dval_u)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_p_dval_c)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_p_dval_t)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_p_dsea_u)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_p_dsea_c)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_p_dsea_t)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_p_ssea_u)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_cc_p_ssea_c)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_cc_p_ssea_t)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_cc_p_bsea_u)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_cc_p_bsea_c)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_cc_p_bsea_t)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_dbar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_sbar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_bbar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_dbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_sbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_bbar) : return kPdgAntiCQuark;   break;

          case (kHEDISQrk_v_cc_n_dval_u)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_n_dval_c)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_n_dval_t)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_n_dsea_u)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_n_dsea_c)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_n_dsea_t)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_cc_n_ssea_u)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_cc_n_ssea_c)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_cc_n_ssea_t)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_cc_n_bsea_u)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_cc_n_bsea_c)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_cc_n_bsea_t)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_dbar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_sbar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_bbar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_dbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_sbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_bbar) : return kPdgAntiCQuark;   break;

          case (kHEDISQrk_vbar_cc_p_uval_d)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_uval_s)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_uval_b)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_usea_d)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_usea_s)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_usea_b)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_csea_d)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_cc_p_csea_s)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_cc_p_csea_b)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_ubar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_cbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_tbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_ubar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_cbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_tbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_ubar) : return kPdgAntiBQuark;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_cbar) : return kPdgAntiBQuark;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_tbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_vbar_cc_n_uval_d)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_uval_s)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_uval_b)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_usea_d)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_usea_s)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_usea_b)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_csea_d)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_cc_n_csea_s)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_cc_n_csea_b)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_ubar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_cbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_tbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_ubar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_cbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_tbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_ubar) : return kPdgAntiBQuark;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_cbar) : return kPdgAntiBQuark;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_tbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_v_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_v_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_vbar_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_vbar_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          default : return 0;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      static int FnlQuarkPdg(HEDISQrkChannel_t channel)
      {
        switch (channel) {
          case (kHEDISQrk_v_cc_p_dval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_p_dval_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_p_dval_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_p_dsea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_p_dsea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_p_dsea_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_p_ssea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_p_ssea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_p_ssea_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_p_bsea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_p_bsea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_p_bsea_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_bbar) : return kPdgAntiBQuark;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_v_cc_n_dval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_n_dval_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_n_dval_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_n_dsea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_n_dsea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_n_dsea_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_n_ssea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_n_ssea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_n_ssea_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_n_bsea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_cc_n_bsea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_cc_n_bsea_t)       : return kPdgTQuark;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_bbar) : return kPdgAntiBQuark;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_vbar_cc_p_uval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_cc_p_uval_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_cc_p_uval_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_cc_p_usea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_cc_p_usea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_cc_p_usea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_cc_p_csea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_cc_p_csea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_cc_p_csea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_tbar) : return kPdgAntiTQuark;   break;

          case (kHEDISQrk_vbar_cc_n_uval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_cc_n_uval_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_cc_n_uval_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_cc_n_usea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_cc_n_usea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_cc_n_usea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_cc_n_csea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_cc_n_csea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_cc_n_csea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_tbar) : return kPdgAntiTQuark;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_tbar) : return kPdgAntiTQuark;   break;

          case (kHEDISQrk_v_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_v_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_v_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_v_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_v_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_v_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_v_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_v_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_v_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_v_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_v_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_vbar_nc_p_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_p_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_p_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_p_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_p_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_nc_p_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_nc_p_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_nc_p_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_nc_p_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_nc_p_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_nc_p_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_nc_p_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          case (kHEDISQrk_vbar_nc_n_dval_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_n_uval_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_n_dsea_d)       : return kPdgDQuark;   break;
          case (kHEDISQrk_vbar_nc_n_usea_u)       : return kPdgUQuark;   break;
          case (kHEDISQrk_vbar_nc_n_ssea_s)       : return kPdgSQuark;   break;
          case (kHEDISQrk_vbar_nc_n_csea_c)       : return kPdgCQuark;   break;
          case (kHEDISQrk_vbar_nc_n_bsea_b)       : return kPdgBQuark;   break;
          case (kHEDISQrk_vbar_nc_n_dbarsea_dbar) : return kPdgAntiDQuark;   break;
          case (kHEDISQrk_vbar_nc_n_ubarsea_ubar) : return kPdgAntiUQuark;   break;
          case (kHEDISQrk_vbar_nc_n_sbarsea_sbar) : return kPdgAntiSQuark;   break;
          case (kHEDISQrk_vbar_nc_n_cbarsea_cbar) : return kPdgAntiCQuark;   break;
          case (kHEDISQrk_vbar_nc_n_bbarsea_bbar) : return kPdgAntiBQuark;   break;

          default : return 0;  break;
        }
        return 0;
      }
      //__________________________________________________________________________
      //__________________________________________________________________________
      static HEDISNucChannel_t HEDISNucChannel(HEDISQrkChannel_t channel)
      {
        switch (channel) {
          case (kHEDISQrk_v_cc_p_dval_u)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_dval_c)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_dval_t)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_dsea_u)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_dsea_c)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_dsea_t)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_ssea_u)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_ssea_c)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_ssea_t)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_bsea_u)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_bsea_c)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_bsea_t)       : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_dbar) : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_sbar) : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_ubarsea_bbar) : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_dbar) : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_sbar) : return kHEDISNuc_v_cc_p;   break;
          case (kHEDISQrk_v_cc_p_cbarsea_bbar) : return kHEDISNuc_v_cc_p;   break;

          case (kHEDISQrk_v_cc_n_dval_u)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_dval_c)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_dval_t)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_dsea_u)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_dsea_c)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_dsea_t)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_ssea_u)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_ssea_c)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_ssea_t)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_bsea_u)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_bsea_c)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_bsea_t)       : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_dbar) : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_sbar) : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_ubarsea_bbar) : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_dbar) : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_sbar) : return kHEDISNuc_v_cc_n;   break;
          case (kHEDISQrk_v_cc_n_cbarsea_bbar) : return kHEDISNuc_v_cc_n;   break;

          case (kHEDISQrk_vbar_cc_p_uval_d)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_uval_s)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_uval_b)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_usea_d)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_usea_s)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_usea_b)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_csea_d)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_csea_s)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_csea_b)       : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_ubar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_cbar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_dbarsea_tbar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_ubar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_cbar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_sbarsea_tbar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_ubar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_cbar) : return kHEDISNuc_vbar_cc_p;   break;
          case (kHEDISQrk_vbar_cc_p_bbarsea_tbar) : return kHEDISNuc_vbar_cc_p;   break;

          case (kHEDISQrk_vbar_cc_n_uval_d)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_uval_s)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_uval_b)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_usea_d)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_usea_s)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_usea_b)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_csea_d)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_csea_s)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_csea_b)       : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_ubar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_cbar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_dbarsea_tbar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_ubar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_cbar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_sbarsea_tbar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_ubar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_cbar) : return kHEDISNuc_vbar_cc_n;   break;
          case (kHEDISQrk_vbar_cc_n_bbarsea_tbar) : return kHEDISNuc_vbar_cc_n;   break;

          case (kHEDISQrk_v_nc_p_dval_d)       : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_uval_u)       : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_dsea_d)       : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_usea_u)       : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_ssea_s)       : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_csea_c)       : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_bsea_b)       : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_dbarsea_dbar) : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_ubarsea_ubar) : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_sbarsea_sbar) : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_cbarsea_cbar) : return kHEDISNuc_v_nc_p;   break;
          case (kHEDISQrk_v_nc_p_bbarsea_bbar) : return kHEDISNuc_v_nc_p;   break;

          case (kHEDISQrk_v_nc_n_dval_d)       : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_uval_u)       : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_dsea_d)       : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_usea_u)       : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_ssea_s)       : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_csea_c)       : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_bsea_b)       : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_dbarsea_dbar) : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_ubarsea_ubar) : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_sbarsea_sbar) : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_cbarsea_cbar) : return kHEDISNuc_v_nc_n;   break;
          case (kHEDISQrk_v_nc_n_bbarsea_bbar) : return kHEDISNuc_v_nc_n;   break;

          case (kHEDISQrk_vbar_nc_p_dval_d)       : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_uval_u)       : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_dsea_d)       : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_usea_u)       : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_ssea_s)       : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_csea_c)       : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_bsea_b)       : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_dbarsea_dbar) : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_ubarsea_ubar) : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_sbarsea_sbar) : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_cbarsea_cbar) : return kHEDISNuc_vbar_nc_p;   break;
          case (kHEDISQrk_vbar_nc_p_bbarsea_bbar) : return kHEDISNuc_vbar_nc_p;   break;

          case (kHEDISQrk_vbar_nc_n_dval_d)       : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_uval_u)       : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_dsea_d)       : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_usea_u)       : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_ssea_s)       : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_csea_c)       : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_bsea_b)       : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_dbarsea_dbar) : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_ubarsea_ubar) : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_sbarsea_sbar) : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_cbarsea_cbar) : return kHEDISNuc_vbar_nc_n;   break;
          case (kHEDISQrk_vbar_nc_n_bbarsea_bbar) : return kHEDISNuc_vbar_nc_n;   break;

          default : return kHEDISNucNull;  break;
        }
        return kHEDISNucNull;
      }
      //__________________________________________________________________________





  };

}      // genie namespace

#endif // _HEDIS_CHANNEL_H_
