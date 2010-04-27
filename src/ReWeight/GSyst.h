//____________________________________________________________________________
/*!

\class    genie::rew::GSyst_t

\brief    An enumeration of systematic parameters

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_SYSTEMATIC_PARAM_H_
#define _G_SYSTEMATIC_PARAM_H_

#include <string>

#include "PDG/PDGUtils.h"
#include "Interaction/InteractionType.h"
#include "HadronTransport/INukeHadroFates.h"

using std::string;

namespace genie {
namespace rew   {

typedef enum EGSyst {

  kSystNull = 0,

  //
  // Neutrino cross section systematics
  // 
  //
  //
  //
  //

  kSystNuXSec_NormCCQE,          ///< tweak CCQE normalization
  kSystNuXSec_MaCCQEshape,       ///< Ma CCQE, affects dsigma(CCQE)/dQ2 - in shape only (normalized to constant integral)
  kSystNuXSec_MaCCQE,            ///< Ma CCQE, affects dsigma(CCQE)/dQ2 - both in shape and normalization
  kSystNuXSec_MvCCQE,            ///< Mv CCQE, affects dsigma(CCQE)/dQ2 - both in shape and normalization
  kSystNuXSec_NormCCRES,         ///< tweak CCRES normalization
  kSystNuXSec_MaCCRESshape,      ///< Ma CCRES, affects d2sigma(CCRES)/dWdQ2 - in shape only (normalized to constant integral)
  kSystNuXSec_MvCCRESshape,      ///< Mv CCRES, affects d2sigma(CCRES)/dWdQ2 - in shape only (normalized to constant integral)
  kSystNuXSec_MaCCRES,           ///< Ma CCRES, affects d2sigma(CCRES)/dWdQ2 - both in shape and normalization
  kSystNuXSec_MvCCRES,           ///< Mv CCRES, affects d2sigma(CCRES)/dWdQ2 - both in shape and normalization
  kSystNuXSec_MaCOHPi,           ///< Ma COHPi
  kSystNuXSec_RvpCC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+p CC
  kSystNuXSec_RvpCC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+p CC
  kSystNuXSec_RvpNC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+p NC
  kSystNuXSec_RvpNC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+p NC
  kSystNuXSec_RvnCC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+n CC
  kSystNuXSec_RvnCC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+n CC
  kSystNuXSec_RvnNC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+n NC
  kSystNuXSec_RvnNC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+n NC
  kSystNuXSec_RvbarpCC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+p CC
  kSystNuXSec_RvbarpCC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+p CC
  kSystNuXSec_RvbarpNC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+p NC
  kSystNuXSec_RvbarpNC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+p NC
  kSystNuXSec_RvbarnCC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+n CC
  kSystNuXSec_RvbarnCC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+n CC
  kSystNuXSec_RvbarnNC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+n NC
  kSystNuXSec_RvbarnNC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+n NC
  kSystNuXSec_NormCCSafeDIS,     ///< tweak CC Safe (Q2>Q2o, W>Wo) DIS normalization, typically Q2o=1GeV^2, Wo=1.7-2.0GeV
  //...

  //
  // Intranuclear rescattering systematics
  // 

  // Parameters controlling the total rescattering probability
  kSystINuke_MFPTwk_pi,          ///< mean free path tweaking factor for pions
  kSystINuke_MFPTwk_N,           ///< mean free path tweaking factor for nucleons

  // Parameters controlling the rescattered hadron fate
  kSystINuke_CExTwk_pi,          ///< charge exchange probability tweaking factor for pions
  kSystINuke_ElTwk_pi,           ///< elastic         probability tweaking factor for pions
  kSystINuke_InelTwk_pi,         ///< inelastic       probability tweaking factor for pions
  kSystINuke_AbsTwk_pi,          ///< absorption      probability tweaking factor for pions
  kSystINuke_PiProdTwk_pi,       ///< pion production probability tweaking factor for pions
  kSystINuke_CExTwk_N,           ///< charge exchange probability tweaking factor for nucleons
  kSystINuke_ElTwk_N,            ///< elastic         probability tweaking factor for nucleons
  kSystINuke_InelTwk_N,          ///< inelastic       probability tweaking factor for nucleons
  kSystINuke_AbsTwk_N,           ///< absorption      probability tweaking factor for nucleons
  kSystINuke_PiProdTwk_N         ///< pion production probability tweaking factor for nucleons
  //...

} GSyst_t;


class GSyst {
public:
 //......................................................................................
 static string AsString(GSyst_t syst) 
 {
   switch(syst) {
     case ( kSystNuXSec_NormCCQE      ) : return "NormCCQE";             break;
     case ( kSystNuXSec_MaCCQE        ) : return "MaCCQE";               break;
     case ( kSystNuXSec_MaCCQEshape   ) : return "MaCCQEshape";          break;
     case ( kSystNuXSec_MvCCQE        ) : return "MvCCQE";               break;
     case ( kSystNuXSec_NormCCRES     ) : return "NormCCRES";            break;
     case ( kSystNuXSec_MaCCRESshape  ) : return "MaCCRESshape";         break;
     case ( kSystNuXSec_MvCCRESshape  ) : return "MvCCRESshape";         break;
     case ( kSystNuXSec_MaCCRES       ) : return "MaCCRES";              break;
     case ( kSystNuXSec_MvCCRES       ) : return "MvCCRES";              break;
     case ( kSystNuXSec_MaCOHPi       ) : return "MaCOHPi";              break;
     case ( kSystNuXSec_RvpCC1pi      ) : return "NonRESBGvpCC1pi";      break;
     case ( kSystNuXSec_RvpCC2pi      ) : return "NonRESBGvpCC2pi";      break;
     case ( kSystNuXSec_RvpNC1pi      ) : return "NonRESBGvpNC1pi";      break;
     case ( kSystNuXSec_RvpNC2pi      ) : return "NonRESBGvpNC2pi";      break;
     case ( kSystNuXSec_RvnCC1pi      ) : return "NonRESBGvpCC1pi";      break;
     case ( kSystNuXSec_RvnCC2pi      ) : return "NonRESBGvpCC2pi";      break;
     case ( kSystNuXSec_RvnNC1pi      ) : return "NonRESBGvpNC1pi";      break;
     case ( kSystNuXSec_RvnNC2pi      ) : return "NonRESBGvpNC2pi";      break;
     case ( kSystNuXSec_RvbarpCC1pi   ) : return "NonRESBGvbarpCC1pi";   break;
     case ( kSystNuXSec_RvbarpCC2pi   ) : return "NonRESBGvbarpCC2pi";   break;
     case ( kSystNuXSec_RvbarpNC1pi   ) : return "NonRESBGvbarpNC1pi";   break;
     case ( kSystNuXSec_RvbarpNC2pi   ) : return "NonRESBGvbarpNC2pi";   break;
     case ( kSystNuXSec_RvbarnCC1pi   ) : return "NonRESBGvbarpCC1pi";   break;
     case ( kSystNuXSec_RvbarnCC2pi   ) : return "NonRESBGvbarpCC2pi";   break;
     case ( kSystNuXSec_RvbarnNC1pi   ) : return "NonRESBGvbarpNC1pi";   break;
     case ( kSystNuXSec_RvbarnNC2pi   ) : return "NonRESBGvbarpNC2pi";   break;
     case ( kSystNuXSec_NormCCSafeDIS ) : return "NormCCSafeDIS";        break;
     case ( kSystINuke_MFPTwk_pi      ) : return "PiMeanFreePathTwk";    break;
     case ( kSystINuke_MFPTwk_N       ) : return "NucMeanFreePathTwk";   break;
     case ( kSystINuke_CExTwk_pi      ) : return "PiCExTwk";             break;
     case ( kSystINuke_ElTwk_pi       ) : return "PiElTwk";              break;
     case ( kSystINuke_InelTwk_pi     ) : return "PiInelTwk";            break;
     case ( kSystINuke_AbsTwk_pi      ) : return "PiAbsTwk";             break;
     case ( kSystINuke_PiProdTwk_pi   ) : return "PiPiProdTwk";          break;
     case ( kSystINuke_CExTwk_N       ) : return "NucCExTwk";            break;
     case ( kSystINuke_ElTwk_N        ) : return "NucElTwk";             break;
     case ( kSystINuke_InelTwk_N      ) : return "NucInelTwk";           break;
     case ( kSystINuke_AbsTwk_N       ) : return "NucAbsTwk";            break;
     case ( kSystINuke_PiProdTwk_N    ) : return "NucPiProdTwk";         break;
     default: 
       return "-";
   }
   return "";
 }
 //......................................................................................
 static bool IsINukePionFateSystematic(GSyst_t syst) 
 {
   switch(syst) {
     case ( kSystINuke_CExTwk_pi   ) : 
     case ( kSystINuke_ElTwk_pi    ) : 
     case ( kSystINuke_InelTwk_pi  ) : 
     case ( kSystINuke_AbsTwk_pi   ) : 
     case ( kSystINuke_PiProdTwk_pi) : 
        return true;
        break;
     default: 
        return false;
        break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeNuclFateSystematic(GSyst_t syst) 
 {
   switch(syst) {
     case ( kSystINuke_CExTwk_N   ) : 
     case ( kSystINuke_ElTwk_N    ) : 
     case ( kSystINuke_InelTwk_N  ) : 
     case ( kSystINuke_AbsTwk_N   ) : 
     case ( kSystINuke_PiProdTwk_N) : 
        return true;
        break;
     default: 
        return false;
        break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeFateSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kSystINuke_CExTwk_pi    ) : 
     case ( kSystINuke_ElTwk_pi     ) :
     case ( kSystINuke_InelTwk_pi   ) :
     case ( kSystINuke_AbsTwk_pi    ) :
     case ( kSystINuke_PiProdTwk_pi ) :
     case ( kSystINuke_CExTwk_N     ) : 
     case ( kSystINuke_ElTwk_N      ) :
     case ( kSystINuke_InelTwk_N    ) :
     case ( kSystINuke_AbsTwk_N     ) :
     case ( kSystINuke_PiProdTwk_N  ) :
       return true;
       break;
     
     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukePionMeanFreePathSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kSystINuke_MFPTwk_pi ) : 
       return true;
       break;
     
     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeNuclMeanFreePathSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kSystINuke_MFPTwk_N  ) : 
       return true;
       break;
     
     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeMeanFreePathSystematic(GSyst_t syst)
 {
   switch(syst) {
     case ( kSystINuke_MFPTwk_pi ) : 
     case ( kSystINuke_MFPTwk_N  ) : 
       return true;
       break;
     
     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static GSyst_t NextPionFateSystematic(int i)
 {
    if(i==0) return kSystINuke_CExTwk_pi;
    if(i==1) return kSystINuke_ElTwk_pi;
    if(i==2) return kSystINuke_InelTwk_pi;
    if(i==3) return kSystINuke_AbsTwk_pi;
    if(i==4) return kSystINuke_PiProdTwk_pi;

    return kSystNull;
 }
 //......................................................................................
 static GSyst_t NextNuclFateSystematic(int i)
 {
    if(i==0) return kSystINuke_CExTwk_N;
    if(i==1) return kSystINuke_ElTwk_N;
    if(i==2) return kSystINuke_InelTwk_N;
    if(i==3) return kSystINuke_AbsTwk_N;
    if(i==4) return kSystINuke_PiProdTwk_N;

    return kSystNull;
 }
 //......................................................................................
 static GSyst_t INukeFate2GSyst(INukeFateHA_t fate, int pdgc)
 {
  // get the corresponding GSyst_t systematic parameter enumeration from the
  // input intranuke fate enumeration and PDG code
  //
  if(pdg::IsPion(pdgc)) {
     switch (fate) {
      case kIHAFtUndefined : return kSystNull;                break;
      case kIHAFtCEx       : return kSystINuke_CExTwk_pi;     break;
      case kIHAFtElas      : return kSystINuke_ElTwk_pi;      break;
      case kIHAFtInelas    : return kSystINuke_InelTwk_pi;    break;
      case kIHAFtAbsNP     : return kSystINuke_AbsTwk_pi;     break;
      case kIHAFtAbsPP     : return kSystINuke_AbsTwk_pi;     break;
      case kIHAFtAbsNPP    : return kSystINuke_AbsTwk_pi;     break;
      case kIHAFtAbsNNP    : return kSystINuke_AbsTwk_pi;     break;
      case kIHAFtAbs2N2P   : return kSystINuke_AbsTwk_pi;     break;
      case kIHAFtAbs2N3P   : return kSystINuke_AbsTwk_pi;     break;
      case kIHAFtNPip      : return kSystINuke_PiProdTwk_pi;  break;
      case kIHAFtNPipPi0   : return kSystINuke_PiProdTwk_pi;  break;
      default              : return kSystNull;                break;
     }
  } else
  if(pdg::IsNucleon(pdgc)) {
     switch (fate) {
      case kIHAFtUndefined : return kSystNull;               break;
      case kIHAFtCEx       : return kSystINuke_CExTwk_N;     break;
      case kIHAFtElas      : return kSystINuke_ElTwk_N;      break;
      case kIHAFtInelas    : return kSystINuke_InelTwk_N;    break;
      case kIHAFtAbsNP     : return kSystINuke_AbsTwk_N;     break;
      case kIHAFtAbsPP     : return kSystINuke_AbsTwk_N;     break;
      case kIHAFtAbsNPP    : return kSystINuke_AbsTwk_N;     break;
      case kIHAFtAbsNNP    : return kSystINuke_AbsTwk_N;     break;
      case kIHAFtAbs2N2P   : return kSystINuke_AbsTwk_N;     break;
      case kIHAFtAbs2N3P   : return kSystINuke_AbsTwk_N;     break;
      case kIHAFtNPip      : return kSystINuke_PiProdTwk_N;  break;
      case kIHAFtNPipPi0   : return kSystINuke_PiProdTwk_N;  break;
      default              : return kSystNull;               break;
     }
  }
  return kSystNull;
 }
 //......................................................................................
 static GSyst_t RBkg(InteractionType_t itype, int probe, int hitnuc, int npi)
 {
   bool is_v    = pdg::IsNeutrino     (probe);
   bool is_vbar = pdg::IsAntiNeutrino (probe);
   bool is_p    = pdg::IsProton       (hitnuc);
   bool is_n    = pdg::IsNeutron      (hitnuc);
  
   // CC
   bool is_cc = (itype == kIntWeakCC);
   if(is_cc) {
     if(is_v && is_p) {
       if(npi==1) return kSystNuXSec_RvpCC1pi;
       if(npi==2) return kSystNuXSec_RvpCC2pi;   
     }
     if(is_v && is_n) {
       if(npi==1) return kSystNuXSec_RvnCC1pi;
       if(npi==2) return kSystNuXSec_RvnCC2pi;   
     }
     if(is_vbar && is_p) {
       if(npi==1) return kSystNuXSec_RvbarpCC1pi;
       if(npi==2) return kSystNuXSec_RvbarpCC2pi;   
     }
     if(is_vbar && is_n) {
       if(npi==1) return kSystNuXSec_RvbarnCC1pi;
       if(npi==2) return kSystNuXSec_RvbarnCC2pi;   
     }
   }//cc

   // NC
   bool is_nc = (itype == kIntWeakNC);
   if(is_nc) {
     if(is_v && is_p) {
       if(npi==1) return kSystNuXSec_RvpNC1pi;
       if(npi==2) return kSystNuXSec_RvpNC2pi;   
     }
     if(is_v && is_n) {
       if(npi==1) return kSystNuXSec_RvnNC1pi;
       if(npi==2) return kSystNuXSec_RvnNC2pi;   
     }
     if(is_vbar && is_p) {
       if(npi==1) return kSystNuXSec_RvbarpNC1pi;
       if(npi==2) return kSystNuXSec_RvbarpNC2pi;   
     }
     if(is_vbar && is_n) {
       if(npi==1) return kSystNuXSec_RvbarnNC1pi;
       if(npi==2) return kSystNuXSec_RvbarnNC2pi;   
     }
   }//nc

   return kSystNull;
 }
 //......................................................................................

};

} // rew   namespace
} // genie namespace

#endif 

