//____________________________________________________________________________ 
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 11, 2009 - CA
   Was first added in v2.5.1. Adapted from the T2K-specific version of the
   GENIE reweighting tool.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Units.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "ReWeight/GReWeightUtils.h"

using namespace genie;

//____________________________________________________________________________
double genie::utils::rew::MeanFreePathWeight(
  int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
  double mfp_scale_factor, bool interacted,
  double nRpi, double nRnuc, double NR, double R0)
{
   // Get the nominal survival probability
   double pdef = utils::intranuke::ProbSurvival(pdgc,x4,p4,A,1.);
   if(pdef<=0) return 1.;

   // Get the survival probability for the tweaked mean free path
   double ptwk = utils::intranuke::ProbSurvival(pdgc,x4,p4,A,mfp_scale_factor);
   if(ptwk<=0) return 1.;

   // Calculate weight
   double w_mfp = utils::rew::MeanFreePathWeight(pdef, ptwk, interacted);
   return w_mfp;
}
//____________________________________________________________________________
double genie::utils::rew::MeanFreePathWeight(
   double pdef, double ptwk, bool interacted)
{
// Returns a weight to account for a change in hadron mean free path inside
// insidea nuclear medium.
//
// Inputs:
//   pdef : nominal survival probability
//   ptwk : survival probability for the tweaked mean free path
//   interacted : flag indicating whether the hadron interacted or escaped
//
// See utils::intranuke::ProbSurvival() for the calculation of probabilities.
//
  double w_mfp = 1.;

  if(interacted) {
     w_mfp = (1-pdef>0) ?  (1-ptwk)  /   (1-pdef)  : 1;
  } else {
     w_mfp = (pdef>0) ? ptwk / pdef : 1;
  }
  w_mfp = TMath::Max(0.,w_mfp);

  return w_mfp;
}
//____________________________________________________________________________
double genie::utils::rew::FateXSec(genie::rew::GSyst_t syst, double kinE)
{
  double fate_xsec = 0.0;
      
  INukeHadroData * hd = INukeHadroData::Instance();

  // convert to MeV and      
  double ke = kinE / units::MeV;         
  ke = TMath::Max(INukeHadroData::fMinKinEnergy,   ke);
  ke = TMath::Min(INukeHadroData::fMaxKinEnergyHA, ke);
         
  switch (syst) {

    //
    // pions
    //

    case (genie::rew::kSystINuke_CExTwk_pi) :
    {
      fate_xsec = hd->XSecPipA_CEx()->Evaluate(ke);
    }
    break;

    case (genie::rew::kSystINuke_ElTwk_pi) :
    {
      fate_xsec = hd->XSecPipA_Elas()->Evaluate(ke);
    }
    break;

    case (genie::rew::kSystINuke_InelTwk_pi) :
    {
      fate_xsec = hd->XSecPipA_Inel()-> Evaluate(ke);
    }
    break;

    case (genie::rew::kSystINuke_AbsTwk_pi) :
    {   
      double abs_pion_np   = hd->XSecPipA_NP  ()->Evaluate(ke);
      double abs_pion_pp   = hd->XSecPipA_PP  ()->Evaluate(ke);
      double abs_pion_npp  = hd->XSecPipA_NPP ()->Evaluate(ke);
      double abs_pion_nnp  = hd->XSecPipA_NNP ()->Evaluate(ke);
      double abs_pion_nnpp = hd->XSecPipA_NNPP()->Evaluate(ke);

      double abs_pion  = abs_pion_np  + 
                         abs_pion_pp  + 
                         abs_pion_npp + 
                         abs_pion_nnp + 
                         abs_pion_nnpp;

      fate_xsec = abs_pion;
    }
    break;

    case (genie::rew::kSystINuke_PiProdTwk_pi) :
    {
      fate_xsec = hd->XSecPipA_NPipPi0()->Evaluate(ke);
    }
    break;
 
    //
    // nucleons
    //

    case (genie::rew::kSystINuke_CExTwk_N) :
    {
      fate_xsec = hd->XSecPA_CEx()->Evaluate(ke);
    }
    break;

    case (genie::rew::kSystINuke_ElTwk_N) :
    {
      fate_xsec = hd->XSecPA_Elas()->Evaluate(ke);
    }
    break;

    case (genie::rew::kSystINuke_InelTwk_N) :
    {
      fate_xsec = hd->XSecPA_Inel()->Evaluate(ke);
    }
    break;

    case (genie::rew::kSystINuke_AbsTwk_N) :
    {
      double abs_nucl_np   = hd->XSecPA_NP   ()->Evaluate(ke);
      double abs_nucl_pp   = hd->XSecPA_PP   ()->Evaluate(ke);
      double abs_nucl_npp  = hd->XSecPA_NPP  ()->Evaluate(ke);
      double abs_nucl_nnp  = hd->XSecPA_NNP  ()->Evaluate(ke);
      double abs_nucl_nnppp= hd->XSecPA_NNPPP()->Evaluate(ke);

      double abs_nucl = abs_nucl_np  + 
                        abs_nucl_pp  + 
                        abs_nucl_npp + 
                        abs_nucl_nnp + 
                        abs_nucl_nnppp;

      fate_xsec = abs_nucl;
    }
    break;

    case (genie::rew::kSystINuke_PiProdTwk_N) :
    {
      double piprod_NPip    = hd->XSecPA_NPip   ()->Evaluate(ke);
      double piprod_NPipPi0 = hd->XSecPA_NPipPi0()->Evaluate(ke);

      fate_xsec = piprod_NPip + piprod_NPipPi0;
    }   
    break;

    default:
    {
      LOG("ReW", pDEBUG)
        << "Have reached default case and assigning cross_section{fate} = 0";
      fate_xsec = 0;
    }
    break;

  } // hadron_fate?
                 
  // the xsection splines in INukeHadroData return the hadron x-section in mb 
  // -> convert to fm^2
  fate_xsec *= (units::mb / units::fm2);
                 
  return fate_xsec;
}
//____________________________________________________________________________

