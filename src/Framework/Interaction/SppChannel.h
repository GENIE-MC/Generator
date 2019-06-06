//____________________________________________________________________________
/*!

\class    genie::SppChannel

\brief    Enumeration of single pion production channels

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  December 16, 2004

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _SPP_CHANNEL_H_
#define _SPP_CHANNEL_H_

#include <string>

#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using std::string;

namespace genie {

typedef enum ESppChannel {

  kSppNull = 0,

         /* [p,n,pi+,pi0,pi-] */

  kSpp_vp_cc_10100,  /* neutrino  CC */
  kSpp_vn_cc_10010,
  kSpp_vn_cc_01100,

  kSpp_vp_nc_10010,  /* neutrino NC */
  kSpp_vp_nc_01100,
  kSpp_vn_nc_01010,
  kSpp_vn_nc_10001,

  kSpp_vbn_cc_01001,  /* anti-neutrino  CC */
  kSpp_vbp_cc_01010,
  kSpp_vbp_cc_10001,

  kSpp_vbp_nc_10010,  /* anti-neutrino NC */
  kSpp_vbp_nc_01100,
  kSpp_vbn_nc_01010,
  kSpp_vbn_nc_10001

} SppChannel_t;


class SppChannel
{
public:

  //__________________________________________________________________________
  static string AsString(SppChannel_t channel)
  {
    switch (channel) {

      case (kSpp_vp_cc_10100) : return "v p -> l- p pi+";   break;
      case (kSpp_vn_cc_10010) : return "v n -> l- p pi0";   break;
      case (kSpp_vn_cc_01100) : return "v n -> l- n pi+";   break; 

      case (kSpp_vp_nc_10010) : return "v p -> v p pi0";    break;
      case (kSpp_vp_nc_01100) : return "v p -> v n pi+";    break;
      case (kSpp_vn_nc_01010) : return "v n -> v n pi0";    break;
      case (kSpp_vn_nc_10001) : return "v n -> v p pi-";    break;

      case (kSpp_vbn_cc_01001): return "vb n -> l+ n pi-";  break;
      case (kSpp_vbp_cc_01010): return "vb p -> l+ n pi0";  break;
      case (kSpp_vbp_cc_10001): return "vb p -> l+ p pi-";  break;

      case (kSpp_vbp_nc_10010): return "vb p -> vb p pi0";  break;
      case (kSpp_vbp_nc_01100): return "vb p -> vb n pi+";  break;
      case (kSpp_vbn_nc_01010): return "vb n -> vb n pi0";  break;
      case (kSpp_vbn_nc_10001): return "vb n -> vb p pi-";  break;

      default : return "Unknown";  break;
    }
    return "Unknown";
  }
  //__________________________________________________________________________
  static int InitStateNucleon(SppChannel_t channel)
  {
    switch (channel) {

      case (kSpp_vp_cc_10100) : return kPdgProton;   break;
      case (kSpp_vn_cc_10010) : return kPdgNeutron;  break;
      case (kSpp_vn_cc_01100) : return kPdgNeutron;  break;

      case (kSpp_vp_nc_10010) : return kPdgProton;   break;
      case (kSpp_vp_nc_01100) : return kPdgProton;   break;
      case (kSpp_vn_nc_01010) : return kPdgNeutron;  break;
      case (kSpp_vn_nc_10001) : return kPdgNeutron;  break;

      case (kSpp_vbn_cc_01001): return kPdgNeutron;  break;
      case (kSpp_vbp_cc_01010): return kPdgProton;   break;
      case (kSpp_vbp_cc_10001): return kPdgProton;   break;

      case (kSpp_vbp_nc_10010): return kPdgProton;   break;
      case (kSpp_vbp_nc_01100): return kPdgProton;   break;
      case (kSpp_vbn_nc_01010): return kPdgNeutron;  break;
      case (kSpp_vbn_nc_10001): return kPdgNeutron;  break;

      default : return 0;  break;
    }
    return 0;
  }
  //__________________________________________________________________________
  static int FinStateNucleon(SppChannel_t channel)
  {
    switch (channel) {

      case (kSpp_vp_cc_10100) : return kPdgProton;   break;
      case (kSpp_vn_cc_10010) : return kPdgProton;   break;
      case (kSpp_vn_cc_01100) : return kPdgNeutron;  break;

      case (kSpp_vp_nc_10010) : return kPdgProton;   break;
      case (kSpp_vp_nc_01100) : return kPdgNeutron;  break;
      case (kSpp_vn_nc_01010) : return kPdgNeutron;  break;
      case (kSpp_vn_nc_10001) : return kPdgProton;   break;

      case (kSpp_vbn_cc_01001): return kPdgNeutron;  break;
      case (kSpp_vbp_cc_01010): return kPdgNeutron;  break;
      case (kSpp_vbp_cc_10001): return kPdgProton;   break;

      case (kSpp_vbp_nc_10010): return kPdgProton;   break;
      case (kSpp_vbp_nc_01100): return kPdgNeutron;  break;
      case (kSpp_vbn_nc_01010): return kPdgNeutron;  break;
      case (kSpp_vbn_nc_10001): return kPdgProton;   break;

      default : return 0;  break;
    }
    return 0;
  }
  //__________________________________________________________________________
  static int FinStatePion(SppChannel_t channel)
  {
    switch (channel) {

      case (kSpp_vp_cc_10100) : return kPdgPiP;  break;
      case (kSpp_vn_cc_10010) : return kPdgPi0;  break;
      case (kSpp_vn_cc_01100) : return kPdgPiP;  break;

      case (kSpp_vp_nc_10010) : return kPdgPi0;  break;
      case (kSpp_vp_nc_01100) : return kPdgPiP;  break;
      case (kSpp_vn_nc_01010) : return kPdgPi0;  break;
      case (kSpp_vn_nc_10001) : return kPdgPiM;  break;

      case (kSpp_vbn_cc_01001): return kPdgPiM;  break;
      case (kSpp_vbp_cc_01010): return kPdgPi0;  break;
      case (kSpp_vbp_cc_10001): return kPdgPiM;  break;

      case (kSpp_vbp_nc_10010): return kPdgPi0;  break;
      case (kSpp_vbp_nc_01100): return kPdgPiP;  break;
      case (kSpp_vbn_nc_01010): return kPdgPi0;  break;
      case (kSpp_vbn_nc_10001): return kPdgPiM;  break;

      default : return 0;  break;
    }
    return 0;
  }
  //__________________________________________________________________________
  static int ResonanceCharge(SppChannel_t channel)
  {
    switch (channel) {

      case (kSpp_vp_cc_10100) : return  2;   break;
      case (kSpp_vn_cc_10010) : return  1;   break;
      case (kSpp_vn_cc_01100) : return  1;   break;

      case (kSpp_vp_nc_10010) : return  1;   break;
      case (kSpp_vp_nc_01100) : return  1;   break;
      case (kSpp_vn_nc_01010) : return  0;   break;
      case (kSpp_vn_nc_10001) : return  0;   break;

      case (kSpp_vbn_cc_01001): return -1;   break;
      case (kSpp_vbp_cc_01010): return  0;   break;
      case (kSpp_vbp_cc_10001): return  0;   break;

      case (kSpp_vbp_nc_10010): return  1;   break;
      case (kSpp_vbp_nc_01100): return  1;   break;
      case (kSpp_vbn_nc_01010): return  0;   break;
      case (kSpp_vbn_nc_10001): return  0;   break;

      default : return 0;  break;
    }
    return 0;
  }
  //__________________________________________________________________________
  static double IsospinWeight(SppChannel_t channel, Resonance_t res)
  {
    // return the isospin Glebsch Gordon coefficient for the input resonance
    // contribution to the input exclusive channel

    bool is_delta = utils::res::IsDelta(res);

    double iw_1_3 = 0.33333333;
    double iw_2_3 = 0.66666666;

    switch (channel) {

      //-- v CC
      case (kSpp_vp_cc_10100) : return (is_delta) ? (1.0)    : (0.0);    break;
      case (kSpp_vn_cc_10010) : return (is_delta) ? (iw_2_3) : (iw_1_3); break;
      case (kSpp_vn_cc_01100) : return (is_delta) ? (iw_1_3) : (iw_2_3); break;

      //-- v NC
      case (kSpp_vp_nc_10010) : return (is_delta) ? (iw_2_3) : (iw_1_3); break;
      case (kSpp_vp_nc_01100) : return (is_delta) ? (iw_1_3) : (iw_2_3); break;
      case (kSpp_vn_nc_01010) : return (is_delta) ? (iw_2_3) : (iw_1_3); break;
      case (kSpp_vn_nc_10001) : return (is_delta) ? (iw_1_3) : (iw_2_3); break;

      //-- same as for neutrinos (? - check)

      //-- vbar CC
      case (kSpp_vbn_cc_01001): return (is_delta) ? (1.0)    : (0.0);    break;
      case (kSpp_vbp_cc_01010): return (is_delta) ? (iw_2_3) : (iw_1_3); break;
      case (kSpp_vbp_cc_10001): return (is_delta) ? (iw_1_3) : (iw_2_3); break;

      //-- vbar NC
      case (kSpp_vbp_nc_10010): return (is_delta) ? (iw_2_3) : (iw_1_3); break;
      case (kSpp_vbp_nc_01100): return (is_delta) ? (iw_1_3) : (iw_2_3); break;
      case (kSpp_vbn_nc_01010): return (is_delta) ? (iw_2_3) : (iw_1_3); break;
      case (kSpp_vbn_nc_10001): return (is_delta) ? (iw_1_3) : (iw_2_3); break;

      default : return 0;  break;
    }

    return 0;
  }
  //__________________________________________________________________________
  // The values of resonance branching ratio is taken from 
  //  K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014)
  static double BranchingRatio(SppChannel_t /*channel*/, Resonance_t res)
  {
    // return the BR for the decay of the input resonance to the final state
    // hadronic system of the input exclusive channel.

    // get list of TDecayChannels, match one with the input channel and get
    // the branching ratio.

    switch(res) {
      case kP33_1232  : return 0.994; break;
      case kS11_1535  : return 0.450; break;
      case kD13_1520  : return 0.600; break; // REMOVE HARDCODED DATA FROM
      case kS11_1650  : return 0.700; break; // HERE AND GET BR's from PDG
      case kD13_1700  : return 0.120; break; // TABLES via TDatabasePDG
      case kD15_1675  : return 0.400; break;
      case kS31_1620  : return 0.250; break;
      case kD33_1700  : return 0.150; break;
      case kP11_1440  : return 0.650; break;
      case kP33_1600  : return 0.175; break;
      case kP13_1720  : return 0.110; break;
      case kF15_1680  : return 0.675; break;
      case kP31_1910  : return 0.225; break;
      case kP33_1920  : return 0.125; break;
      case kF35_1905  : return 0.120; break;
      case kF37_1950  : return 0.400; break;
      case kP11_1710  : return 0.125; break;
      case kF17_1970  : return 0.150; break;
      default: break;
    }
    return 0;
  }
  //__________________________________________________________________________
  static SppChannel_t FromInteraction(const Interaction * interaction)
  {
    const InitialState & init_state = interaction->InitState();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    const XclsTag &      xcls_tag   = interaction->ExclTag();

    if( xcls_tag.NPions()    != 1 ) return kSppNull;
    if( xcls_tag.NNucleons() != 1 ) return kSppNull;

    // get struck nucleon
    int hit_nucl_pdgc = init_state.Tgt().HitNucPdg();
    if( ! pdg::IsNeutronOrProton(hit_nucl_pdgc) ) return kSppNull;
    bool hit_p = pdg::IsProton(hit_nucl_pdgc);
    bool hit_n = !hit_p;

    // the final state hadronic sytem has 1 pi and 1 nucleon
    bool fs_pi_plus  = ( xcls_tag.NPiPlus()   == 1 );
    bool fs_pi_minus = ( xcls_tag.NPiMinus()  == 1 );
    bool fs_pi_0     = ( xcls_tag.NPi0()      == 1 );
    bool fs_p        = ( xcls_tag.NProtons()  == 1 );
    bool fs_n        = ( xcls_tag.NNeutrons() == 1 );

    // get probe
    int probe = init_state.ProbePdg();

    // figure out spp channel
    if( pdg::IsNeutrino(probe) ) { 

       if ( proc_info.IsWeakCC() ) {
          if      (hit_p && fs_p && fs_pi_plus ) return kSpp_vp_cc_10100;
          else if (hit_n && fs_p && fs_pi_0    ) return kSpp_vn_cc_10010;
          else if (hit_n && fs_n && fs_pi_plus ) return kSpp_vn_cc_01100;
          else                                   return kSppNull;
       } else if ( proc_info.IsWeakNC() ) {
          if      (hit_p && fs_p && fs_pi_0    ) return kSpp_vp_nc_10010;
          else if (hit_p && fs_n && fs_pi_plus ) return kSpp_vp_nc_01100;
          else if (hit_n && fs_n && fs_pi_0    ) return kSpp_vn_nc_01010;
          else if (hit_n && fs_p && fs_pi_minus) return kSpp_vn_nc_10001;
          else                                   return kSppNull;
       } else return kSppNull;

    } else if( pdg::IsAntiNeutrino(probe) ) { 

       if ( proc_info.IsWeakCC() ) {
          if      (hit_n && fs_n && fs_pi_minus) return kSpp_vbn_cc_01001;
          else if (hit_p && fs_n && fs_pi_0    ) return kSpp_vbp_cc_01010;
          else if (hit_p && fs_p && fs_pi_minus) return kSpp_vbp_cc_10001;
          else                                   return kSppNull;
       } else if ( proc_info.IsWeakNC() ) {
          if      (hit_p && fs_p && fs_pi_0    ) return kSpp_vbp_nc_10010;
          else if (hit_p && fs_n && fs_pi_plus ) return kSpp_vbp_nc_01100;
          else if (hit_n && fs_n && fs_pi_0    ) return kSpp_vbn_nc_01010;
          else if (hit_n && fs_p && fs_pi_minus) return kSpp_vbn_nc_10001;
          else                                   return kSppNull;
       } else return kSppNull;
    }

    return kSppNull;
  }
  //__________________________________________________________________________
};

}      // genie namespace

#endif // _SPP_CHANNEL_H_
