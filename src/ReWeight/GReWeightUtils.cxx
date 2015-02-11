//____________________________________________________________________________ 
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 11, 2009 - CA
   Was first added in v2.5.1. Adapted from the T2K-specific version of the
   GENIE reweighting tool.
 @ Dec 17, 2010 - JD
   Added method to calculate weight for a modified formation zone. 
 @ Jul 29, 2011 - SD,AM
   Update INUKE fates. Mean free path is now function of Z too.   
 @ Feb 08, 2013 - CA
   Adjust formation zone reweighting. Mean free path is function of Z too.
*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Conventions/Units.h"
#include "Conventions/Controls.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightUtils.h"

using namespace genie;
using namespace genie::rew;
using namespace genie::controls;

//____________________________________________________________________________
double genie::utils::rew::MeanFreePathWeight(
  int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, 
  double A, double Z,
  double mfp_scale_factor, bool interacted,
  double nRpi, double nRnuc, double NR, double R0)
{
   LOG("ReW", pINFO)  
     << "Calculating mean free path weight: "
     << "A = " << A << ", Z = " << Z << ", mfp_scale = " << mfp_scale_factor
     << ", interacted = " << interacted;
   LOG("ReW", pDEBUG)  
     << "nR_pion = " << nRpi << ", nR_nucleon = " << nRnuc 
     << ", NR = " << NR << ", R0 = " << R0;

   // Get the nominal survival probability
   double pdef = utils::intranuke::ProbSurvival(
      pdgc,x4,p4,A,Z,1.,nRpi,nRnuc,NR,R0);
   LOG("ReW", pINFO)  << "Probability(default mfp) = " << pdef;      
   if(pdef<=0) return 1.;

   // Get the survival probability for the tweaked mean free path
   double ptwk = utils::intranuke::ProbSurvival(
      pdgc,x4,p4,A,Z,mfp_scale_factor,nRpi,nRnuc,NR,R0);
   LOG("ReW", pINFO)  << "Probability(tweaked mfp) = " << ptwk;      
   if(ptwk<=0) return 1.;

   // Calculate weight
   double w_mfp = utils::rew::MeanFreePathWeight(pdef, ptwk, interacted);
   LOG("ReW", pINFO)  << "Mean free path weight = " << w_mfp;     
   return w_mfp;
}
//____________________________________________________________________________
double genie::utils::rew::FZoneWeight(
  int pdgc, const TLorentzVector & vtx, const TLorentzVector & x4, 
  const TLorentzVector & p4, double A, double Z, 
  double fz_scale_factor, bool interacted,
  double nRpi, double nRnuc, double NR, double R0)
{
   // Calculate hadron start assuming tweaked formation zone
   TLorentzVector fz    = x4 - vtx;
   TLorentzVector fztwk = fz_scale_factor*fz;
   TLorentzVector x4twk = x4 + fztwk - fz;

   LOG("ReW", pDEBUG)  << "Formation zone = "<< fz.Vect().Mag() << " fm";

   // Get nominal survival probability.
   double pdef = utils::intranuke::ProbSurvival(
      pdgc,x4,p4,A,Z,1.,nRpi,nRnuc,NR,R0);
   LOG("ReW", pDEBUG)  << "Survival probability (nominal) = "<< pdef;
   if(pdef<=0) return 1.; 
   if(pdef>=1.){
     LOG("ReW", pERROR) 
         <<  "Default formation zone takes hadron outside "
         <<  "nucleus so cannot reweight!" ;
     return 1.;
   }

   // Get tweaked survival probability.
   double ptwk = utils::intranuke::ProbSurvival(
      pdgc,x4twk,p4,A,Z,1.,nRpi,nRnuc,NR,R0);
   if(ptwk<=0) return 1.;
   LOG("ReW", pDEBUG)  << "Survival probability (tweaked) = "<< ptwk;

   // Calculate weight
   double w_fz = utils::rew::MeanFreePathWeight(pdef, ptwk, interacted);
   LOG("ReW", pDEBUG)  
       << "Particle weight for formation zone tweak = "<< ptwk;
   return w_fz;  
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
double genie::utils::rew::FateFraction(
   genie::rew::GSyst_t syst, double kinE, double frac_scale_factor)
{
  double fate_frac = 0.0;
      
  INukeHadroData * hd = INukeHadroData::Instance();

  // convert to MeV and      
  double ke = kinE / units::MeV;         
  ke = TMath::Max(INukeHadroData::fMinKinEnergy,   ke);
  ke = TMath::Min(INukeHadroData::fMaxKinEnergyHA, ke);
         
  switch (syst) {

    //
    // pions
    //

    case (genie::rew::kINukeTwkDial_FrCEx_pi) :
    {
      fate_frac = hd->Frac(kPdgPiP, kIHAFtCEx, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrElas_pi) :
    {
      fate_frac = hd->Frac(kPdgPiP, kIHAFtElas, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrInel_pi) :
    {
      fate_frac = hd->Frac(kPdgPiP, kIHAFtInelas, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrAbs_pi) :
    {  
      fate_frac = hd->Frac(kPdgPiP, kIHAFtAbs, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrPiProd_pi) :
    {
      fate_frac = hd->Frac(kPdgPiP, kIHAFtPiProd,  ke);
    }
    break;
 
    //
    // nucleons
    //

    case (genie::rew::kINukeTwkDial_FrCEx_N) :
    {
      fate_frac = hd->Frac(kPdgProton, kIHAFtCEx, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrElas_N) :
    {
      fate_frac = hd->Frac(kPdgProton, kIHAFtElas, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrInel_N) :
    {
      fate_frac = hd->Frac(kPdgProton, kIHAFtInelas, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrAbs_N) :
    {
      fate_frac = hd->Frac(kPdgProton, kIHAFtAbs,    ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrPiProd_N) :
    {
      fate_frac = hd->Frac(kPdgProton, kIHAFtPiProd,  ke);
    }   
    break;

    default:
    {
      LOG("ReW", pDEBUG)
        << "Have reached default case and assigning fraction{fate} = 0";
      fate_frac = 0;
    }
    break;

  } // hadron_fate?

  fate_frac *= frac_scale_factor;
                 
  return fate_frac;
}
//____________________________________________________________________________
double genie::utils::rew::WhichFateFractionScaleFactor(
    genie::rew::GSyst_t syst, double kinE, double fate_frac)
{  
  double fate_frac_nominal = FateFraction(syst,kinE,1.);

  if(TMath::Abs(fate_frac-fate_frac_nominal) < kASmallNum) return 0;

  if(fate_frac_nominal <= 0) { return -99999; }

  double scale = TMath::Max(0.,fate_frac)/fate_frac_nominal;
  return scale;
}
//____________________________________________________________________________
bool genie::utils::rew::HadronizedByAGKY(const EventRecord & event)
{ 
  Interaction * interaction = event.Summary();
  assert(interaction);                                                     

  bool is_dis  = interaction->ProcInfo().IsDeepInelastic();
  bool charm   = interaction->ExclTag().IsCharmEvent();
  bool by_agky = is_dis && !charm;

  return by_agky;
}
//____________________________________________________________________________
bool genie::utils::rew::HadronizedByAGKYPythia(const EventRecord & event) 
{ 
// Check whether the event was hadronized by AGKY/KNO or AGKY/PYTHIA
   
  GHepStatus_t prefragm = kIStDISPreFragmHadronicState;
  bool found_string  = (event.FindParticle(kPdgString,  prefragm, 0) != 0);
  bool found_cluster = (event.FindParticle(kPdgCluster, prefragm, 0) != 0);
  bool handled_by_pythia = found_string || found_cluster;

  return handled_by_pythia;
}
//____________________________________________________________________________
TLorentzVector genie::utils::rew::Hadronic4pLAB(const EventRecord & event)
{
  GHepParticle * nu = event.Probe();                    // incoming v
  GHepParticle * N  = event.HitNucleon();               // struck nucleon
  GHepParticle * l  = event.FinalStatePrimaryLepton();  // f/s primary lepton
 
  assert(nu);                                                              
  assert(N);
  assert(l);

  // Compute the Final State Hadronic System 4p (PX = Pv + PN - Pl)

  const TLorentzVector & p4nu = *(nu->P4());
  const TLorentzVector & p4N  = *(N ->P4());
  const TLorentzVector & p4l  = *(l ->P4());

  TLorentzVector pX4 = p4nu + p4N - p4l;

  return pX4;
}
//____________________________________________________________________________
double genie::utils::rew::AGKYWeight(int /*pdgc*/, double /*xF*/, double /*pT2*/)
{
  return 1.0;
}
//____________________________________________________________________________
int genie::utils::rew::Sign(double twkdial)
{
  if(twkdial < 0.) return -1;
  if(twkdial > 0.) return +1;
  return 0;
}
//____________________________________________________________________________

