//____________________________________________________________________________
/*!

\class    genie::SppChannel

\brief    Enumeration of single pion production channels

\authors   Costas Andreopoulos <c.andreopoulos \at cern.ch>
           University of Liverpool \n
           Igor Kakorin <kakorin@jinr.ru> Joint Institute for Nuclear Research \n

\created  December 16, 2004

\update   November 12, 2019
          Added extra functions for MK model. \n
          Branching ratios are looked in particle database now.

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _SPP_CHANNEL_H_
#define _SPP_CHANNEL_H_

#include <string>

#include <TDecayChannel.h>

#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Conventions/Constants.h"

using std::string;
using namespace genie::constants;

namespace genie {

typedef enum ESppChannel {

  kSppNull = 0,

  // [p,n,pi+,pi0,pi-]

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
  static int FinStateIsospin(SppChannel_t channel)
  {
    switch (channel) {

      case (kSpp_vp_cc_10100) : return  3;   break;
      case (kSpp_vn_cc_10010) : return  1;   break;
      case (kSpp_vn_cc_01100) : return  1;   break;

      case (kSpp_vp_nc_10010) : return  1;   break;
      case (kSpp_vp_nc_01100) : return  1;   break;
      case (kSpp_vn_nc_01010) : return  1;   break;
      case (kSpp_vn_nc_10001) : return  1;   break;

      case (kSpp_vbn_cc_01001): return  3;   break;
      case (kSpp_vbp_cc_01010): return  1;   break;
      case (kSpp_vbp_cc_10001): return  1;   break;

      case (kSpp_vbp_nc_10010): return  1;   break;
      case (kSpp_vbp_nc_01100): return  1;   break;
      case (kSpp_vbn_nc_01010): return  1;   break;
      case (kSpp_vbn_nc_10001): return  1;   break;

      default : return 0;  break;
    }
    return 0;
  }
  //__________________________________________________________________________
  static double IsospinWeight(SppChannel_t channel, Resonance_t res)
  {
    // return the square of isospin Glebsch Gordon coefficient for the input resonance
    // contribution to the input exclusive channel

    bool is_delta = utils::res::IsDelta(res);

    double iw_1_3 = 1./3;
    double iw_2_3 = 2./3;

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
  static double Isospin3Coefficients(SppChannel_t channel)
  {
    // return the isospin coefficients for the channel
    // with final state isospin = 3/2

    // [p,n,pi+,pi0,pi-]
    switch (channel) {

      //-- v CC
      case (kSpp_vp_cc_10100) : return  kSqrt3;
      case (kSpp_vn_cc_10010) : return -kSqrt2_3;
      case (kSpp_vn_cc_01100) : return  k1_Sqrt3;

      //-- v NC
      case (kSpp_vp_nc_10010) : return  kSqrt2_3;
      case (kSpp_vp_nc_01100) : return -k1_Sqrt3;
      case (kSpp_vn_nc_01010) : return  kSqrt2_3;
      case (kSpp_vn_nc_10001) : return  k1_Sqrt3;



      //-- vbar CC
      case (kSpp_vbn_cc_01001): return  kSqrt3;
      case (kSpp_vbp_cc_01010): return -kSqrt2_3;
      case (kSpp_vbp_cc_10001): return  k1_Sqrt3;
                                
      //-- vbar NC              
      case (kSpp_vbp_nc_10010): return  kSqrt2_3;
      case (kSpp_vbp_nc_01100): return -k1_Sqrt3;
      case (kSpp_vbn_nc_01010): return  kSqrt2_3;
      case (kSpp_vbn_nc_10001): return  k1_Sqrt3;

      default : return 0;
    }

  }
  //__________________________________________________________________________
  static double Isospin1Coefficients(SppChannel_t channel)
  {
    // return the isospin coefficients for the channel
    // with final state isospin = 1/2

    // [p,n,pi+,pi0,pi-]
    switch (channel) {

      //-- v CC
      case (kSpp_vp_cc_10100) : return  0.;
      case (kSpp_vn_cc_10010) : return  k1_Sqrt3;
      case (kSpp_vn_cc_01100) : return  kSqrt2_3;

      //-- v NC
      case (kSpp_vp_nc_10010) : return -k1_Sqrt3;
      case (kSpp_vp_nc_01100) : return -kSqrt2_3;
      case (kSpp_vn_nc_01010) : return  k1_Sqrt3;
      case (kSpp_vn_nc_10001) : return -kSqrt2_3;



      //-- vbar CC
      case (kSpp_vbn_cc_01001): return  0.;
      case (kSpp_vbp_cc_01010): return  k1_Sqrt3;
      case (kSpp_vbp_cc_10001): return  kSqrt2_3;
                                
      //-- vbar NC              
      case (kSpp_vbp_nc_10010): return -k1_Sqrt3;
      case (kSpp_vbp_nc_01100): return -kSqrt2_3;
      case (kSpp_vbn_nc_01010): return  k1_Sqrt3;
      case (kSpp_vbn_nc_10001): return -kSqrt2_3;

      default : return 0;
    }

  }
  //__________________________________________________________________________
  // The values of resonance mass and width is taken from
  // M. Tanabashi et al. (Particle Data Group) Phys. Rev. D 98, 030001
  // Hardcoded data are removed, now they are taken from PDG table via TDatabasePDG and cached
  static double BranchingRatio(Resonance_t res)
  {
    // return the BR for the decay of the input resonance to the final state: nucleon + pion.
    // get list of TDecayChannels, match one with the input channel and get
    // the branching ratio.
    static std::map<Resonance_t, double> cache ;

    auto it = cache.find( res ) ;
    if ( it != cache.end() ) {
      return it -> second ;
    }

    double BR = 0. ;
    
    PDGLibrary * pdglib = PDGLibrary::Instance();
    // the charge of resonance does not matter
    int pdg = genie::utils::res::PdgCode(res, 0);
    TParticlePDG * res_pdg = pdglib->Find( pdg );
    if (res_pdg != 0)
    {
      for (int nch = 0; nch < res_pdg->NDecayChannels(); nch++)
      {
        TDecayChannel * ch = res_pdg->DecayChannel(nch);
        if (ch->NDaughters() == 2)
        {
          int first_daughter_pdg  = ch->DaughterPdgCode (0);
          int second_daughter_pdg = ch->DaughterPdgCode (1);
          if ((genie::pdg::IsNucleon(first_daughter_pdg ) && genie::pdg::IsPion(second_daughter_pdg)) ||
              (genie::pdg::IsNucleon(second_daughter_pdg) && genie::pdg::IsPion(first_daughter_pdg )))
	    {
	      BR += ch->BranchingRatio();
          }
        }
      }
      cache[res] = BR;
      return BR;
    }

    // should not be here - meaningless to return anything
    gAbortingInErr = true;
    LOG("SppChannel", pFATAL) << "Unknown resonance " << res;
    exit(1);

  }
  //__________________________________________________________________________
  static SppChannel_t FromInteraction(const Interaction * interaction)
  {
    const InitialState & init_state = interaction->InitState();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    if ( !proc_info.IsSinglePion() ) return kSppNull;
    
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
