//____________________________________________________________________________ 
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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

#include <cassert>

#include <TMath.h>

#include "Conventions/Units.h"
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

//____________________________________________________________________________
double genie::utils::rew::MeanFreePathWeight(
  int pdgc, const TLorentzVector & x4, const TLorentzVector & p4, double A,
  double mfp_scale_factor, bool interacted,
  double nRpi, double nRnuc, double NR, double R0)
{
   // Get the nominal survival probability
   double pdef = utils::intranuke::ProbSurvival(
      pdgc,x4,p4,A,1.,nRpi,nRnuc,NR,R0);
   if(pdef<=0) return 1.;

   // Get the survival probability for the tweaked mean free path
   double ptwk = utils::intranuke::ProbSurvival(
      pdgc,x4,p4,A,mfp_scale_factor,nRpi,nRnuc,NR,R0);
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
      fate_frac = 0.; 
      fate_frac += hd->Frac(kPdgPiP, kIHAFtAbsNP,   ke);
      fate_frac += hd->Frac(kPdgPiP, kIHAFtAbsPP,   ke);
      fate_frac += hd->Frac(kPdgPiP, kIHAFtAbsNPP,  ke);
      fate_frac += hd->Frac(kPdgPiP, kIHAFtAbsNNP,  ke);
      fate_frac += hd->Frac(kPdgPiP, kIHAFtAbs2N2P, ke);
    //fate_frac += hd->Frac(kPdgPiP, kIHAFtAbs2N3P, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrPiProd_pi) :
    {
      fate_frac = 0.; 
    //fate_frac += hd->Frac(kPdgPiP, kIHAFtNPip,    ke);
      fate_frac += hd->Frac(kPdgPiP, kIHAFtNPipPi0, ke);
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
      fate_frac = 0.; 
      fate_frac += hd->Frac(kPdgProton, kIHAFtAbsNP,   ke);
      fate_frac += hd->Frac(kPdgProton, kIHAFtAbsPP,   ke);
      fate_frac += hd->Frac(kPdgProton, kIHAFtAbsNPP,  ke);
      fate_frac += hd->Frac(kPdgProton, kIHAFtAbsNNP,  ke);
    //fate_frac += hd->Frac(kPdgProton, kIHAFtAbs2N2P, ke);
      fate_frac += hd->Frac(kPdgProton, kIHAFtAbs2N3P, ke);
    }
    break;

    case (genie::rew::kINukeTwkDial_FrPiProd_N) :
    {
      fate_frac = 0.; 
      fate_frac += hd->Frac(kPdgProton, kIHAFtNPip,    ke);
      fate_frac += hd->Frac(kPdgProton, kIHAFtNPipPi0, ke);
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

