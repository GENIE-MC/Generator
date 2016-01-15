//____________________________________________________________________________
/*!

\class    genie::rew::GSyst_t

\brief    An enumeration of systematic parameters

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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

  // NCEL tweaking parameters:
  kXSecTwkDial_MaNCEL,            ///< tweak Ma NCEL, affects dsigma(NCEL)/dQ2 both in shape and normalization
  kXSecTwkDial_EtaNCEL,           ///< tweak NCEL strange axial form factor eta, affects dsigma(NCEL)/dQ2 both in shape and normalization
  // CCQE tweaking parameters:
  kXSecTwkDial_NormCCQE,          ///< tweak CCQE normalization (energy independent)
  kXSecTwkDial_NormCCQEenu,       ///< tweak CCQE normalization (maintains dependence on neutrino energy)
  kXSecTwkDial_MaCCQEshape,       ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaCCQE,            ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 both in shape and normalization
  kXSecTwkDial_VecFFCCQEshape,    ///< tweak elastic nucleon form factors (BBA/default -> dipole) - shape only effect of dsigma(CCQE)/dQ2
  // Resonance neutrino-production tweaking parameters:
  kXSecTwkDial_NormCCRES,         ///< tweak CCRES normalization
  kXSecTwkDial_MaCCRESshape,      ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MvCCRESshape,      ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaCCRES,           ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_MvCCRES,           ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_NormNCRES,         ///< tweak NCRES normalization
  kXSecTwkDial_MaNCRESshape,      ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MvNCRESshape,      ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaNCRES,           ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_MvNCRES,           ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
  // Coherent pion production tweaking parameters:
  kXSecTwkDial_MaCOHpi,           ///< tweak Ma for COH pion production
  kXSecTwkDial_R0COHpi,           ///< tweak R0 for COH pion production
  // Non-resonance background tweaking parameters:
  kXSecTwkDial_RvpCC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+p CC
  kXSecTwkDial_RvpCC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+p CC
  kXSecTwkDial_RvpNC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+p NC
  kXSecTwkDial_RvpNC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+p NC
  kXSecTwkDial_RvnCC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+n CC
  kXSecTwkDial_RvnCC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+n CC
  kXSecTwkDial_RvnNC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+n NC
  kXSecTwkDial_RvnNC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+n NC
  kXSecTwkDial_RvbarpCC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p CC
  kXSecTwkDial_RvbarpCC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p CC
  kXSecTwkDial_RvbarpNC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p NC
  kXSecTwkDial_RvbarpNC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p NC
  kXSecTwkDial_RvbarnCC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n CC
  kXSecTwkDial_RvbarnCC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n CC
  kXSecTwkDial_RvbarnNC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n NC
  kXSecTwkDial_RvbarnNC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n NC
  // DIS tweaking parameters - applied for DIS events with (Q2>Q2o, W>Wo), typically Q2o=1GeV^2, Wo=1.7-2.0GeV
  kXSecTwkDial_AhtBY,             ///< tweak the Bodek-Yang model parameter A_{ht} - incl. both shape and normalization effect
  kXSecTwkDial_BhtBY,             ///< tweak the Bodek-Yang model parameter B_{ht} - incl. both shape and normalization effect 
  kXSecTwkDial_CV1uBY,            ///< tweak the Bodek-Yang model parameter CV1u - incl. both shape and normalization effect 
  kXSecTwkDial_CV2uBY,            ///< tweak the Bodek-Yang model parameter CV2u - incl. both shape and normalization effect 
  kXSecTwkDial_AhtBYshape,        ///< tweak the Bodek-Yang model parameter A_{ht} - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_BhtBYshape,        ///< tweak the Bodek-Yang model parameter B_{ht} - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_CV1uBYshape,       ///< tweak the Bodek-Yang model parameter CV1u - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_CV2uBYshape,       ///< tweak the Bodek-Yang model parameter CV2u - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_NormDISCC,         ///< tweak the inclusive DIS CC normalization
  kXSecTwkDial_RnubarnuCC,        ///< tweak the ratio of \sigma(\bar\nu CC) / \sigma(\nu CC)
  kXSecTwkDial_DISNuclMod,        ///< tweak DIS nuclear modification (shadowing, anti-shadowing, EMC)
  //
  kXSecTwkDial_NC,                ///< 


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
     case ( kXSecTwkDial_MaNCEL           ) : return "MaNCEL";               break;
     case ( kXSecTwkDial_EtaNCEL          ) : return "EtaNCEL";              break;
     case ( kXSecTwkDial_NormCCQE         ) : return "NormCCQE";             break;
     case ( kXSecTwkDial_NormCCQEenu      ) : return "NormCCQEenu";          break;
     case ( kXSecTwkDial_MaCCQE           ) : return "MaCCQE";               break;
     case ( kXSecTwkDial_MaCCQEshape      ) : return "MaCCQEshape";          break;
     case ( kXSecTwkDial_VecFFCCQEshape   ) : return "VecFFCCQEshape";       break;
     case ( kXSecTwkDial_NormCCRES        ) : return "NormCCRES";            break;
     case ( kXSecTwkDial_MaCCRESshape     ) : return "MaCCRESshape";         break;
     case ( kXSecTwkDial_MvCCRESshape     ) : return "MvCCRESshape";         break;
     case ( kXSecTwkDial_MaCCRES          ) : return "MaCCRES";              break;
     case ( kXSecTwkDial_MvCCRES          ) : return "MvCCRES";              break;
     case ( kXSecTwkDial_NormNCRES        ) : return "NormNCRES";            break;
     case ( kXSecTwkDial_MaNCRESshape     ) : return "MaNCRESshape";         break;
     case ( kXSecTwkDial_MvNCRESshape     ) : return "MvNCRESshape";         break;
     case ( kXSecTwkDial_MaNCRES          ) : return "MaNCRES";              break;
     case ( kXSecTwkDial_MvNCRES          ) : return "MvNCRES";              break;
     case ( kXSecTwkDial_MaCOHpi          ) : return "MaCOHpi";              break;
     case ( kXSecTwkDial_R0COHpi          ) : return "R0COHpi";              break;
     case ( kXSecTwkDial_RvpCC1pi         ) : return "NonRESBGvpCC1pi";      break;
     case ( kXSecTwkDial_RvpCC2pi         ) : return "NonRESBGvpCC2pi";      break;
     case ( kXSecTwkDial_RvpNC1pi         ) : return "NonRESBGvpNC1pi";      break;
     case ( kXSecTwkDial_RvpNC2pi         ) : return "NonRESBGvpNC2pi";      break;
     case ( kXSecTwkDial_RvnCC1pi         ) : return "NonRESBGvnCC1pi";      break;
     case ( kXSecTwkDial_RvnCC2pi         ) : return "NonRESBGvnCC2pi";      break;
     case ( kXSecTwkDial_RvnNC1pi         ) : return "NonRESBGvnNC1pi";      break;
     case ( kXSecTwkDial_RvnNC2pi         ) : return "NonRESBGvnNC2pi";      break;
     case ( kXSecTwkDial_RvbarpCC1pi      ) : return "NonRESBGvbarpCC1pi";   break;
     case ( kXSecTwkDial_RvbarpCC2pi      ) : return "NonRESBGvbarpCC2pi";   break;
     case ( kXSecTwkDial_RvbarpNC1pi      ) : return "NonRESBGvbarpNC1pi";   break;
     case ( kXSecTwkDial_RvbarpNC2pi      ) : return "NonRESBGvbarpNC2pi";   break;
     case ( kXSecTwkDial_RvbarnCC1pi      ) : return "NonRESBGvbarnCC1pi";   break;
     case ( kXSecTwkDial_RvbarnCC2pi      ) : return "NonRESBGvbarnCC2pi";   break;
     case ( kXSecTwkDial_RvbarnNC1pi      ) : return "NonRESBGvbarnNC1pi";   break;
     case ( kXSecTwkDial_RvbarnNC2pi      ) : return "NonRESBGvbarnNC2pi";   break;
     case ( kXSecTwkDial_AhtBY            ) : return "AhtBY";                break;
     case ( kXSecTwkDial_BhtBY            ) : return "BhtBY";                break;
     case ( kXSecTwkDial_CV1uBY           ) : return "CV1uBY";               break;
     case ( kXSecTwkDial_CV2uBY           ) : return "CV2uBY";               break;
     case ( kXSecTwkDial_AhtBYshape       ) : return "AhtBYshape";           break;
     case ( kXSecTwkDial_BhtBYshape       ) : return "BhtBYshape";           break;
     case ( kXSecTwkDial_CV1uBYshape      ) : return "CV1uBYshape";          break;
     case ( kXSecTwkDial_CV2uBYshape      ) : return "CV2uBYshape";          break;
     case ( kXSecTwkDial_NormDISCC        ) : return "NormDISCC";            break;
     case ( kXSecTwkDial_RnubarnuCC       ) : return "RnubarnuCC";           break;
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
 static GSyst_t FromString(string syst)
 {
   GSyst_t systematics[] = 
   {
       kXSecTwkDial_MaNCEL,
       kXSecTwkDial_EtaNCEL,
       kXSecTwkDial_NormCCQE,   
       kXSecTwkDial_NormCCQEenu,   
       kXSecTwkDial_MaCCQE,        
       kXSecTwkDial_MaCCQEshape,   
       kXSecTwkDial_VecFFCCQEshape,
       kXSecTwkDial_NormCCRES,    
       kXSecTwkDial_MaCCRESshape, 
       kXSecTwkDial_MvCCRESshape, 
       kXSecTwkDial_MaCCRES,      
       kXSecTwkDial_MvCCRES,      
       kXSecTwkDial_NormNCRES,    
       kXSecTwkDial_MaNCRESshape, 
       kXSecTwkDial_MvNCRESshape, 
       kXSecTwkDial_MaNCRES,      
       kXSecTwkDial_MvNCRES,      
       kXSecTwkDial_MaCOHpi,      
       kXSecTwkDial_R0COHpi,    
       kXSecTwkDial_RvpCC1pi,   
       kXSecTwkDial_RvpCC2pi,   
       kXSecTwkDial_RvpNC1pi,    
       kXSecTwkDial_RvpNC2pi,     
       kXSecTwkDial_RvnCC1pi,     
       kXSecTwkDial_RvnCC2pi,     
       kXSecTwkDial_RvnNC1pi,     
       kXSecTwkDial_RvnNC2pi,     
       kXSecTwkDial_RvbarpCC1pi,  
       kXSecTwkDial_RvbarpCC2pi,  
       kXSecTwkDial_RvbarpNC1pi,  
       kXSecTwkDial_RvbarpNC2pi,  
       kXSecTwkDial_RvbarnCC1pi,  
       kXSecTwkDial_RvbarnCC2pi,  
       kXSecTwkDial_RvbarnNC1pi,  
       kXSecTwkDial_RvbarnNC2pi,  
       kXSecTwkDial_AhtBY,        
       kXSecTwkDial_BhtBY,        
       kXSecTwkDial_CV1uBY,       
       kXSecTwkDial_CV2uBY,       
       kXSecTwkDial_AhtBYshape,   
       kXSecTwkDial_BhtBYshape,   
       kXSecTwkDial_CV1uBYshape,  
       kXSecTwkDial_CV2uBYshape,  
       kXSecTwkDial_NormDISCC,    
       kXSecTwkDial_RnubarnuCC,   
       kXSecTwkDial_DISNuclMod,   
       kHadrAGKYTwkDial_xF1pi,    
       kHadrAGKYTwkDial_pT1pi,    
       kHadrNuclTwkDial_FormZone, 
       kINukeTwkDial_MFP_pi,      
       kINukeTwkDial_MFP_N,       
       kINukeTwkDial_FrCEx_pi,    
       kINukeTwkDial_FrElas_pi,   
       kINukeTwkDial_FrInel_pi,   
       kINukeTwkDial_FrAbs_pi,    
       kINukeTwkDial_FrPiProd_pi, 
       kINukeTwkDial_FrCEx_N,    
       kINukeTwkDial_FrElas_N,   
       kINukeTwkDial_FrInel_N,   
       kINukeTwkDial_FrAbs_N,    
       kINukeTwkDial_FrPiProd_N, 
       kSystNucl_CCQEPauliSupViaKF,   
       kSystNucl_CCQEMomDistroFGtoSF,
       kRDcyTwkDial_BR1gamma,       
       kRDcyTwkDial_BR1eta,         
       kRDcyTwkDial_Theta_Delta2Npi,
       kNullSystematic
   };

   int isyst=0;
   while(systematics[isyst]!=kNullSystematic) {
     if( AsString(systematics[isyst]) == syst ) {
        return systematics[isyst];
     }
     isyst++;
   }

   return kNullSystematic;
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
      case kIHAFtUndefined : return kNullSystematic;            break;
      case kIHAFtCEx       : return kINukeTwkDial_FrCEx_pi;     break;
      case kIHAFtElas      : return kINukeTwkDial_FrElas_pi;    break;
      case kIHAFtInelas    : return kINukeTwkDial_FrInel_pi;    break;
      case kIHAFtAbs       : return kINukeTwkDial_FrAbs_pi;     break;
      case kIHAFtPiProd    : return kINukeTwkDial_FrPiProd_pi;  break;
      default              : return kNullSystematic;            break;
     }
  } else
  if(pdg::IsNucleon(pdgc)) {
     switch (fate) {
      case kIHAFtUndefined : return kNullSystematic;           break;
      case kIHAFtCEx       : return kINukeTwkDial_FrCEx_N;     break;
      case kIHAFtElas      : return kINukeTwkDial_FrElas_N;    break;
      case kIHAFtInelas    : return kINukeTwkDial_FrInel_N;    break;
      case kIHAFtAbs       : return kINukeTwkDial_FrAbs_N;     break;
      case kIHAFtPiProd    : return kINukeTwkDial_FrPiProd_N;  break;
      default              : return kNullSystematic;           break;
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

