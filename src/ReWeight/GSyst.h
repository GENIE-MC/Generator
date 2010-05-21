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

  kNullSystematic = 0,

  //
  // Neutrino cross section systematics
  // 
  // Note: 
  // 
  //
  //

  kXSecTwkDial_NormCCQE,          ///< tweak CCQE normalization
  kXSecTwkDial_MaCCQEshape,       ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 - in shape only (normalized to constant integral)
  kXSecTwkDial_MaCCQE,            ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 - both in shape and normalization
  kXSecTwkDial_NormCCRES,         ///< tweak CCRES normalization
  kXSecTwkDial_MaCCRESshape,      ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 - in shape only (normalized to constant integral)
  kXSecTwkDial_MvCCRESshape,      ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 - in shape only (normalized to constant integral)
  kXSecTwkDial_MaCCRES,           ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 - both in shape and normalization
  kXSecTwkDial_MvCCRES,           ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 - both in shape and normalization
  kXSecTwkDial_MaCOHpi,           ///< tweak Ma for COH pion production
  kXSecTwkDial_R0COHpi,           ///< tweak R0 for COH pion production
  kXSecTwkDial_RvpCC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+p CC
  kXSecTwkDial_RvpCC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+p CC
  kXSecTwkDial_RvpNC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+p NC
  kXSecTwkDial_RvpNC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+p NC
  kXSecTwkDial_RvnCC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+n CC
  kXSecTwkDial_RvnCC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+n CC
  kXSecTwkDial_RvnNC1pi,          ///< controls the 1pi non-RES bkg in the RES region, for v+n NC
  kXSecTwkDial_RvnNC2pi,          ///< controls the 2pi non-RES bkg in the RES region, for v+n NC
  kXSecTwkDial_RvbarpCC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+p CC
  kXSecTwkDial_RvbarpCC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+p CC
  kXSecTwkDial_RvbarpNC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+p NC
  kXSecTwkDial_RvbarpNC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+p NC
  kXSecTwkDial_RvbarnCC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+n CC
  kXSecTwkDial_RvbarnCC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+n CC
  kXSecTwkDial_RvbarnNC1pi,       ///< controls the 1pi non-RES bkg in the RES region, for vbar+n NC
  kXSecTwkDial_RvbarnNC2pi,       ///< controls the 2pi non-RES bkg in the RES region, for vbar+n NC
  kXSecTwkDial_NormCCSafeDIS,     ///< tweak CC Safe (Q2>Q2o, W>Wo) DIS normalization, typically Q2o=1GeV^2, Wo=1.7-2.0GeV
  kXSecTwkDial_DISNuclMod,        ///< tweak DIS nuclear modification (shadowing, anti-shadowing, EMC)


  //
  // Hadronization (free nucleon target)
  // 

  kHadrAGKYTwkDial_xF1pi,         ///< tweak xF distribution for low multiplicity (N + pi) DIS f/s produced by AGKY
  kHadrAGKYTwkDial_pT1pi,         ///< tweak pT distribution for low multiplicity (N + pi) DIS f/s produced by AGKY


  //
  // Medium-effects to hadronization
  // 

  kHadrNuclTwkDial_FormZone,         ///< tweak formation zone


  //
  // Intranuclear rescattering systematics.
  // There are 2 sets of parameters:
  // - parameters that control the total rescattering probability, P(total)
  // - parameters that control the fraction of each process (`fate'), given a total rescat. prob., P(fate|total)
  // These parameters are considered separately for pions and nucleons.
  //

  kINukeTwkDial_MFP_pi,      ///< tweak mean free path for pions
  kINukeTwkDial_MFP_N,       ///< tweak mean free path for nucleons
  kINukeTwkDial_FrCEx_pi,    ///< tweak charge exchange probability for pions, for given total rescattering probability
  kINukeTwkDial_FrElas_pi,   ///< tweak elastic         probability for pions, for given total rescattering probability
  kINukeTwkDial_FrInel_pi,   ///< tweak inelastic       probability for pions, for given total rescattering probability
  kINukeTwkDial_FrAbs_pi,    ///< tweak absorption      probability for pions, for given total rescattering probability
  kINukeTwkDial_FrPiProd_pi, ///< tweak pion production probability for pions, for given total rescattering probability
  kINukeTwkDial_FrCEx_N,     ///< tweak charge exchange probability for nucleons, for given total rescattering probability
  kINukeTwkDial_FrElas_N,    ///< tweak elastic         probability for nucleons, for given total rescattering probability
  kINukeTwkDial_FrInel_N,    ///< tweak inelastic       probability for nucleons, for given total rescattering probability
  kINukeTwkDial_FrAbs_N,     ///< tweak absorption      probability for nucleons, for given total rescattering probability
  kINukeTwkDial_FrPiProd_N,  ///< tweak pion production probability for nucleons, for given total rescattering probability

  //
  // Nuclear model
  // 

  kSystNucl_CCQEPauliSupViaKF,   ///<
  kSystNucl_CCQEMomDistroFGtoSF, ///<

  //
  // Resonance decays
  // 

  kRDcyTwkDial_BR1gamma,        ///< tweak Resonance -> X + gamma branching ratio, eg Delta+(1232) -> p gamma
  kRDcyTwkDial_BR1eta,          ///< tweak Resonance -> X + eta   branching ratio, eg N+(1440) -> p eta
  kRDcyTwkDial_Theta_Delta2Npi  ///< distort pi angular distribution in Delta -> N + pi


  //
  // Misc
  // 


} GSyst_t;


class GSyst {
public:
 //......................................................................................
 static string AsString(GSyst_t syst) 
 {
   switch(syst) {
     case ( kXSecTwkDial_NormCCQE         ) : return "NormCCQE";             break;
     case ( kXSecTwkDial_MaCCQE           ) : return "MaCCQE";               break;
     case ( kXSecTwkDial_MaCCQEshape      ) : return "MaCCQEshape";          break;
     case ( kXSecTwkDial_NormCCRES        ) : return "NormCCRES";            break;
     case ( kXSecTwkDial_MaCCRESshape     ) : return "MaCCRESshape";         break;
     case ( kXSecTwkDial_MvCCRESshape     ) : return "MvCCRESshape";         break;
     case ( kXSecTwkDial_MaCCRES          ) : return "MaCCRES";              break;
     case ( kXSecTwkDial_MvCCRES          ) : return "MvCCRES";              break;
     case ( kXSecTwkDial_MaCOHpi          ) : return "MaCOHpi";              break;
     case ( kXSecTwkDial_R0COHpi          ) : return "R0COHpi";              break;
     case ( kXSecTwkDial_RvpCC1pi         ) : return "NonRESBGvpCC1pi";      break;
     case ( kXSecTwkDial_RvpCC2pi         ) : return "NonRESBGvpCC2pi";      break;
     case ( kXSecTwkDial_RvpNC1pi         ) : return "NonRESBGvpNC1pi";      break;
     case ( kXSecTwkDial_RvpNC2pi         ) : return "NonRESBGvpNC2pi";      break;
     case ( kXSecTwkDial_RvnCC1pi         ) : return "NonRESBGvpCC1pi";      break;
     case ( kXSecTwkDial_RvnCC2pi         ) : return "NonRESBGvpCC2pi";      break;
     case ( kXSecTwkDial_RvnNC1pi         ) : return "NonRESBGvpNC1pi";      break;
     case ( kXSecTwkDial_RvnNC2pi         ) : return "NonRESBGvpNC2pi";      break;
     case ( kXSecTwkDial_RvbarpCC1pi      ) : return "NonRESBGvbarpCC1pi";   break;
     case ( kXSecTwkDial_RvbarpCC2pi      ) : return "NonRESBGvbarpCC2pi";   break;
     case ( kXSecTwkDial_RvbarpNC1pi      ) : return "NonRESBGvbarpNC1pi";   break;
     case ( kXSecTwkDial_RvbarpNC2pi      ) : return "NonRESBGvbarpNC2pi";   break;
     case ( kXSecTwkDial_RvbarnCC1pi      ) : return "NonRESBGvbarpCC1pi";   break;
     case ( kXSecTwkDial_RvbarnCC2pi      ) : return "NonRESBGvbarpCC2pi";   break;
     case ( kXSecTwkDial_RvbarnNC1pi      ) : return "NonRESBGvbarpNC1pi";   break;
     case ( kXSecTwkDial_RvbarnNC2pi      ) : return "NonRESBGvbarpNC2pi";   break;
     case ( kXSecTwkDial_NormCCSafeDIS    ) : return "NormCCSafeDIS";        break;
     case ( kXSecTwkDial_DISNuclMod       ) : return "DISNuclMod";           break;
     case ( kHadrAGKYTwkDial_xF1pi        ) : return "AGKYxF1pi";            break;
     case ( kHadrAGKYTwkDial_pT1pi        ) : return "AGKYpT1pi";            break;
     case ( kHadrNuclTwkDial_FormZone     ) : return "FormZone";             break;
     case ( kINukeTwkDial_MFP_pi          ) : return "MFP_pi";               break;
     case ( kINukeTwkDial_MFP_N           ) : return "MFP_N";                break;
     case ( kINukeTwkDial_FrCEx_pi        ) : return "FrCEx_pi";             break;
     case ( kINukeTwkDial_FrElas_pi       ) : return "FrElas_pi";            break;
     case ( kINukeTwkDial_FrInel_pi       ) : return "FrInel_pi";            break;
     case ( kINukeTwkDial_FrAbs_pi        ) : return "FrAbs_pi";             break;
     case ( kINukeTwkDial_FrPiProd_pi     ) : return "FrPiProd_pi";          break;
     case ( kINukeTwkDial_FrCEx_N         ) : return "FrCEx_N";              break;
     case ( kINukeTwkDial_FrElas_N        ) : return "FrElas_N";             break;
     case ( kINukeTwkDial_FrInel_N        ) : return "FrInel_N";             break;
     case ( kINukeTwkDial_FrAbs_N         ) : return "FrAbs_N";              break;
     case ( kINukeTwkDial_FrPiProd_N      ) : return "FrPiProd_N";           break;
     case ( kSystNucl_CCQEPauliSupViaKF   ) : return "CCQEPauliSupViaKF";    break;
     case ( kSystNucl_CCQEMomDistroFGtoSF ) : return "CCQEMomDistroFGtoSF";  break;
     case ( kRDcyTwkDial_BR1gamma         ) : return "RDecBR1gamma";         break;
     case ( kRDcyTwkDial_BR1eta           ) : return "RDecBR1eta";           break;
     case ( kRDcyTwkDial_Theta_Delta2Npi  ) : return "Theta_Delta2Npi";      break;

     default: 
       return "-";
   }
   return "";
 }
 //......................................................................................
 static bool IsINukePionFateSystematic(GSyst_t syst) 
 {
   switch(syst) {
     case ( kINukeTwkDial_FrCEx_pi   ) : 
     case ( kINukeTwkDial_FrElas_pi  ) : 
     case ( kINukeTwkDial_FrInel_pi  ) : 
     case ( kINukeTwkDial_FrAbs_pi   ) : 
     case ( kINukeTwkDial_FrPiProd_pi) : 
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
     case ( kINukeTwkDial_FrCEx_N   ) : 
     case ( kINukeTwkDial_FrElas_N  ) : 
     case ( kINukeTwkDial_FrInel_N  ) : 
     case ( kINukeTwkDial_FrAbs_N   ) : 
     case ( kINukeTwkDial_FrPiProd_N) : 
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
     case ( kINukeTwkDial_FrCEx_pi    ) : 
     case ( kINukeTwkDial_FrElas_pi   ) :
     case ( kINukeTwkDial_FrInel_pi   ) :
     case ( kINukeTwkDial_FrAbs_pi    ) :
     case ( kINukeTwkDial_FrPiProd_pi ) :
     case ( kINukeTwkDial_FrCEx_N     ) : 
     case ( kINukeTwkDial_FrElas_N    ) :
     case ( kINukeTwkDial_FrInel_N    ) :
     case ( kINukeTwkDial_FrAbs_N     ) :
     case ( kINukeTwkDial_FrPiProd_N  ) :
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
     case ( kINukeTwkDial_MFP_pi ) : 
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
     case ( kINukeTwkDial_MFP_N  ) : 
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
     case ( kINukeTwkDial_MFP_pi ) : 
     case ( kINukeTwkDial_MFP_N  ) : 
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
    if(i==0) return kINukeTwkDial_FrCEx_pi;
    if(i==1) return kINukeTwkDial_FrElas_pi;
    if(i==2) return kINukeTwkDial_FrInel_pi;
    if(i==3) return kINukeTwkDial_FrAbs_pi;
    if(i==4) return kINukeTwkDial_FrPiProd_pi;

    return kNullSystematic;
 }
 //......................................................................................
 static GSyst_t NextNuclFateSystematic(int i)
 {
    if(i==0) return kINukeTwkDial_FrCEx_N;
    if(i==1) return kINukeTwkDial_FrElas_N;
    if(i==2) return kINukeTwkDial_FrInel_N;
    if(i==3) return kINukeTwkDial_FrAbs_N;
    if(i==4) return kINukeTwkDial_FrPiProd_N;

    return kNullSystematic;
 }
 //......................................................................................
 static GSyst_t INukeFate2GSyst(INukeFateHA_t fate, int pdgc)
 {
  // get the corresponding GSyst_t systematic parameter enumeration from the
  // input intranuke fate enumeration and PDG code
  //
  if(pdg::IsPion(pdgc)) {
     switch (fate) {
      case kIHAFtUndefined : return kNullSystematic;                  break;
      case kIHAFtCEx       : return kINukeTwkDial_FrCEx_pi;     break;
      case kIHAFtElas      : return kINukeTwkDial_FrElas_pi;    break;
      case kIHAFtInelas    : return kINukeTwkDial_FrInel_pi;    break;
      case kIHAFtAbsNP     : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtAbsPP     : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtAbsNPP    : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtAbsNNP    : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtAbs2N2P   : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtAbs2N3P   : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtNPip      : return kINukeTwkDial_FrPiProd_pi;  break;
      case kIHAFtNPipPi0   : return kINukeTwkDial_FrPiProd_pi;  break;
      default              : return kNullSystematic;                  break;
     }
  } else
  if(pdg::IsNucleon(pdgc)) {
     switch (fate) {
      case kIHAFtUndefined : return kNullSystematic;                 break;
      case kIHAFtCEx       : return kINukeTwkDial_FrCEx_N;     break;
      case kIHAFtElas      : return kINukeTwkDial_FrElas_N;    break;
      case kIHAFtInelas    : return kINukeTwkDial_FrInel_N;    break;
      case kIHAFtAbsNP     : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtAbsPP     : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtAbsNPP    : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtAbsNNP    : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtAbs2N2P   : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtAbs2N3P   : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtNPip      : return kINukeTwkDial_FrPiProd_N;  break;
      case kIHAFtNPipPi0   : return kINukeTwkDial_FrPiProd_N;  break;
      default              : return kNullSystematic;                 break;
     }
  }
  return kNullSystematic;
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
       if(npi==1) return kXSecTwkDial_RvpCC1pi;
       if(npi==2) return kXSecTwkDial_RvpCC2pi;   
     }
     if(is_v && is_n) {
       if(npi==1) return kXSecTwkDial_RvnCC1pi;
       if(npi==2) return kXSecTwkDial_RvnCC2pi;   
     }
     if(is_vbar && is_p) {
       if(npi==1) return kXSecTwkDial_RvbarpCC1pi;
       if(npi==2) return kXSecTwkDial_RvbarpCC2pi;   
     }
     if(is_vbar && is_n) {
       if(npi==1) return kXSecTwkDial_RvbarnCC1pi;
       if(npi==2) return kXSecTwkDial_RvbarnCC2pi;   
     }
   }//cc

   // NC
   bool is_nc = (itype == kIntWeakNC);
   if(is_nc) {
     if(is_v && is_p) {
       if(npi==1) return kXSecTwkDial_RvpNC1pi;
       if(npi==2) return kXSecTwkDial_RvpNC2pi;   
     }
     if(is_v && is_n) {
       if(npi==1) return kXSecTwkDial_RvnNC1pi;
       if(npi==2) return kXSecTwkDial_RvnNC2pi;   
     }
     if(is_vbar && is_p) {
       if(npi==1) return kXSecTwkDial_RvbarpNC1pi;
       if(npi==2) return kXSecTwkDial_RvbarpNC2pi;   
     }
     if(is_vbar && is_n) {
       if(npi==1) return kXSecTwkDial_RvbarnNC1pi;
       if(npi==2) return kXSecTwkDial_RvbarnNC2pi;   
     }
   }//nc

   return kNullSystematic;
 }
 //......................................................................................

};

} // rew   namespace
} // genie namespace

#endif 

