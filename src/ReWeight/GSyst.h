//____________________________________________________________________________
/*!

\class    genie::rew::GSyst_t

\brief    An enumeration of systematic parameters

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _G_SYSTEMATIC_PARAM_H_
#define _G_SYSTEMATIC_PARAM_H_

#include <string>

#include "ReWeight/GSystType.h"

using std::string;

namespace genie {
namespace rew   {

typedef enum EGSyst {

  kSystNull = 0,

  //
  // *** Neutrino cross section systematics
  // 

  kSystNuXSec_MaQEL,             ///< Ma QEL
  kSystNuXSec_MvQEL,             ///< Mv QEL
  kSystNuXSec_MaRES,             ///< Ma RES
  kSystNuXSec_MvRES,             ///< Mv RES
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
  //...

  //
  // *** Intranuclear rescattering systematics
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
     case ( kSystNuXSec_MaQEL       ) : return "MaQEL";                break;
     case ( kSystNuXSec_MvQEL       ) : return "MvQEL";                break;
     case ( kSystNuXSec_MaRES       ) : return "MaRES";                break;
     case ( kSystNuXSec_MvRES       ) : return "MvRES";                break;
     case ( kSystNuXSec_MaCOHPi     ) : return "MaCOHPi";              break;
     case ( kSystNuXSec_RvpCC1pi    ) : return "NonRESBGvpCC1pi";      break;
     case ( kSystNuXSec_RvpCC2pi    ) : return "NonRESBGvpCC2pi";      break;
     case ( kSystNuXSec_RvpNC1pi    ) : return "NonRESBGvpNC1pi";      break;
     case ( kSystNuXSec_RvpNC2pi    ) : return "NonRESBGvpNC2pi";      break;
     case ( kSystNuXSec_RvnCC1pi    ) : return "NonRESBGvpCC1pi";      break;
     case ( kSystNuXSec_RvnCC2pi    ) : return "NonRESBGvpCC2pi";      break;
     case ( kSystNuXSec_RvnNC1pi    ) : return "NonRESBGvpNC1pi";      break;
     case ( kSystNuXSec_RvnNC2pi    ) : return "NonRESBGvpNC2pi";      break;
     case ( kSystNuXSec_RvbarpCC1pi ) : return "NonRESBGvbarpCC1pi";   break;
     case ( kSystNuXSec_RvbarpCC2pi ) : return "NonRESBGvbarpCC2pi";   break;
     case ( kSystNuXSec_RvbarpNC1pi ) : return "NonRESBGvbarpNC1pi";   break;
     case ( kSystNuXSec_RvbarpNC2pi ) : return "NonRESBGvbarpNC2pi";   break;
     case ( kSystNuXSec_RvbarnCC1pi ) : return "NonRESBGvbarpCC1pi";   break;
     case ( kSystNuXSec_RvbarnCC2pi ) : return "NonRESBGvbarpCC2pi";   break;
     case ( kSystNuXSec_RvbarnNC1pi ) : return "NonRESBGvbarpNC1pi";   break;
     case ( kSystNuXSec_RvbarnNC2pi ) : return "NonRESBGvbarpNC2pi";   break;
     case ( kSystINuke_MFPTwk_pi    ) : return "PiMeanFreePathTwk";    break;
     case ( kSystINuke_MFPTwk_N     ) : return "NucMeanFreePathTwk";   break;
     case ( kSystINuke_CExTwk_pi    ) : return "PiCExTwk";             break;
     case ( kSystINuke_ElTwk_pi     ) : return "PiElTwk";              break;
     case ( kSystINuke_InelTwk_pi   ) : return "PiInelTwk";            break;
     case ( kSystINuke_AbsTwk_pi    ) : return "PiAbsTwk";             break;
     case ( kSystINuke_PiProdTwk_pi ) : return "PiPiProdTwk";          break;
     case ( kSystINuke_CExTwk_N     ) : return "NucCExTwk";            break;
     case ( kSystINuke_ElTwk_N      ) : return "NucElTwk";             break;
     case ( kSystINuke_InelTwk_N    ) : return "NucInelTwk";           break;
     case ( kSystINuke_AbsTwk_N     ) : return "NucAbsTwk";            break;
     case ( kSystINuke_PiProdTwk_N  ) : return "NucPiProdTwk";         break;
     default: 
       return "-";
   }
   return "";
 }
 //......................................................................................
 static GSystType_t Type(GSyst_t syst)
 {
   // report the 'type of systematic'
   //
   switch(syst) {
     case ( kSystNuXSec_MaQEL       ) : 
     case ( kSystNuXSec_MvQEL       ) : 
     case ( kSystNuXSec_MaRES       ) : 
     case ( kSystNuXSec_MvRES       ) : 
     case ( kSystNuXSec_MaCOHPi     ) : 
     case ( kSystNuXSec_RvpCC1pi    ) : 
     case ( kSystNuXSec_RvpCC2pi    ) : 
     case ( kSystNuXSec_RvpNC1pi    ) : 
     case ( kSystNuXSec_RvpNC2pi    ) : 
     case ( kSystNuXSec_RvnCC1pi    ) : 
     case ( kSystNuXSec_RvnCC2pi    ) : 
     case ( kSystNuXSec_RvnNC1pi    ) : 
     case ( kSystNuXSec_RvnNC2pi    ) : 
     case ( kSystNuXSec_RvbarpCC1pi ) : 
     case ( kSystNuXSec_RvbarpCC2pi ) : 
     case ( kSystNuXSec_RvbarpNC1pi ) : 
     case ( kSystNuXSec_RvbarpNC2pi ) : 
     case ( kSystNuXSec_RvbarnCC1pi ) :  
     case ( kSystNuXSec_RvbarnCC2pi ) :  
     case ( kSystNuXSec_RvbarnNC1pi ) :  
     case ( kSystNuXSec_RvbarnNC2pi ) :  

       return kSystType_NuXSec;
       break;

     case ( kSystINuke_MFPTwk_pi    ) : 
     case ( kSystINuke_MFPTwk_N     ) : 
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

       return kSystType_INuke;
       break;
     
     default:
       return kSystType_Null;
       break;
   }
   return kSystType_Null;
 }
 //......................................................................................
 static bool IntranukePionFateSystematic(GSyst_t syst) 
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
 static bool IntranukeNucleonFateSystematic(GSyst_t syst) 
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

};

} // rew   namespace
} // genie namespace

#endif 

