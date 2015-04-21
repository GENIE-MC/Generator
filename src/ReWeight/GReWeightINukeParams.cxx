//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
         Imperial College London

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.
 @ Mar 09, 2010 - JD
   Protect against negative hadron fate cross section
 @ Feb 09, 2011 - JD
   Remove condition that require a non-zero tweak dial when calling SetCurTwkDial
   as this can lead to old tweak dial value being used inadvertently.
 @ Feb 08, 2013 - CA
   Almost complete re-write of GReWeightINukeParams::Fates to improve readability
   and make sure it handles all special cases correctly.

*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>

#include <TMath.h>
#include <TLorentzVector.h>

#include "Conventions/Controls.h"
#include "ReWeight/GReWeightINukeParams.h"
#include "ReWeight/GReWeightUtils.h"
#include "Messenger/Messenger.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Utils/NuclearUtils.h"
#include "Numerical/Spline.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//___________________________________________________________________________ 
GReWeightINukeParams::GReWeightINukeParams(void)
{
  this->Init();
}
//___________________________________________________________________________ 
GReWeightINukeParams::~GReWeightINukeParams(void)
{

}
//___________________________________________________________________________ 
void GReWeightINukeParams::Init(void)
{
  fParmPionFates = new Fates (kRwINukePion);
  fParmNuclFates = new Fates (kRwINukeNucl);
  fParmPionMFP   = new MFP   (kRwINukePion);
  fParmNuclMFP   = new MFP   (kRwINukeNucl);
}
//___________________________________________________________________________ 
GReWeightINukeParams::Fates * 
  GReWeightINukeParams::FateParams(int pdgc) const
{
  if(pdg::IsPion   (pdgc)) return fParmPionFates;
  if(pdg::IsNucleon(pdgc)) return fParmNuclFates;
  return 0;
}
//___________________________________________________________________________ 
GReWeightINukeParams::MFP * 
  GReWeightINukeParams::MeanFreePathParams(int pdgc) const
{
  if(pdg::IsPion   (pdgc)) return fParmPionMFP;
  if(pdg::IsNucleon(pdgc)) return fParmNuclMFP;
  return 0;
}
//___________________________________________________________________________ 
void GReWeightINukeParams::Reset(void)
{
  fParmPionFates -> Reset();
  fParmNuclFates -> Reset();
  fParmPionMFP   -> Reset();
  fParmNuclMFP   -> Reset();
}
//___________________________________________________________________________ 
void GReWeightINukeParams::Reconfigure(void)
{
  LOG("ReW", pINFO) 
   << "Reconfiguring weight calculator for pion intranuke fate systematics";
  fParmPionFates -> Reconfigure();

  LOG("ReW", pINFO) 
   << "Reconfiguring weight calculator for nucleon intranuke fate systematics";
  fParmNuclFates -> Reconfigure();
}
//___________________________________________________________________________ 
double GReWeightINukeParams::ChisqPenalty(void) const
{
  double chisq = 0.;

  chisq += (fParmPionFates -> ChisqPenalty());
  chisq += (fParmNuclFates -> ChisqPenalty());
  chisq += (fParmPionMFP   -> ChisqPenalty());
  chisq += (fParmNuclMFP   -> ChisqPenalty());

  return chisq;
}
//___________________________________________________________________________ 
void GReWeightINukeParams::SetTwkDial(GSyst_t syst, double val)
{
  if(GSyst::IsINukePionFateSystematic(syst)) 
  {
    fParmPionFates->SetTwkDial(syst,val);
  }
  else
  if(GSyst::IsINukeNuclFateSystematic(syst))
  {
    fParmNuclFates->SetTwkDial(syst,val);
  }
  else
  if(GSyst::IsINukePionMeanFreePathSystematic(syst))
  {
    fParmPionMFP->SetTwkDial(val);
  }
  else
  if(GSyst::IsINukeNuclMeanFreePathSystematic(syst))
  {
    fParmNuclMFP->SetTwkDial(val);
  } 
}
//___________________________________________________________________________ 
//
// Fates nested class 
//
//___________________________________________________________________________ 
GReWeightINukeParams::Fates::Fates(GReWeightINukeParams::HadronType_t ht) 
{
  if (ht != kRwINukePion && ht != kRwINukeNucl) {
    LOG("ReW", pFATAL) 
      << "Mean-free path reweighting code handles only pions or nucleons";
    exit(1);
  }

  fHadType = ht;

  this->Reset();
}
//___________________________________________________________________________
GReWeightINukeParams::Fates::~Fates(void)
{
  this->Reset();
}
//___________________________________________________________________________
void GReWeightINukeParams::Fates::Reset()
{
  fSystValuesUser.clear();
  fSystValuesActual.clear();
  fIsCushion.clear();
}
//___________________________________________________________________________
void GReWeightINukeParams::Fates::Reconfigure()
{
  this->AddCushionTerms();
}
//___________________________________________________________________________
void GReWeightINukeParams::Fates::SetTwkDial(GSyst_t syst, double val)
{
  // check type of systematic 
  if(!this->IsHandled(syst)) return;

  // can not explicitly set params designated as cushion terms
  if(this->IsIncluded(syst)) {
    if(this->IsCushionTerm(syst)) {
       LOG("ReW", pWARN) 
         << "You may not set the value of cushion term " << GSyst::AsString(syst)
         << ". It is set automatically to maintain unitarity";
       return;
    }
  }

  // update tweaking dial
  fSystValuesUser[syst] = val;
  fIsCushion[syst]   = false;
}
//___________________________________________________________________________
double GReWeightINukeParams::Fates::ScaleFactor(
      GSyst_t syst, const TLorentzVector & p4) const
{
  double KE = p4.Energy() - p4.M(); // kinetic energy
  return this->ScaleFactor(syst, KE);
}
//___________________________________________________________________________
double GReWeightINukeParams::Fates::ScaleFactor(
      GSyst_t syst, double KE) const
{
  GSystUncertainty * uncert = GSystUncertainty::Instance();

  double fractional_error = uncert->OneSigmaErr(syst);
  double twk_dial         = this->ActualTwkDial(syst, KE);

  double fate_fraction_scale = 1. + twk_dial * fractional_error;

  LOG("ReW", pNOTICE) 
      << "Systematic " << GSyst::AsString(syst) 
      << " (1 sigma frac. err = " << 100*fractional_error
      << "%, tweak dial = " << twk_dial
      << "): Weight = " << fate_fraction_scale;

  return fate_fraction_scale;
}
//___________________________________________________________________________
double GReWeightINukeParams::Fates::ActualTwkDial(GSyst_t syst, double KE) const
{
  if(! this->IsIncluded(syst)) {
    LOG("ReW", pWARN) 
      << "Systematic " << GSyst::AsString(syst) << " not included";
    return 0.;
  }
  if(KE < 0.) {
    LOG("ReW", pERROR) << "Negative kinetic energy! (KE = " << KE << ")";
    return 0.;
  }

  GSystUncertainty * uncert = GSystUncertainty::Instance();
  map<GSyst_t, double>::const_iterator iter;  

  //
  // initialize actual systematic parameter values to the ones set by the user
  //
  fSystValuesActual.clear();
  fSystValuesActual.insert(fSystValuesUser.begin(), fSystValuesUser.end());

  //
  // Check the tweaking dials set by the user.
  // If a systematic variation is too extreme so that the scale factor for 
  // the corresponding fate fraction becomes negative, automatically re-adjust
  // the systematic variation so that the scale factor is exactly 0
  //

  for(iter = fSystValuesActual.begin(); iter != fSystValuesActual.end(); ++iter)
  {
     GSyst_t curr_syst     = iter->first;
     double  curr_twk_dial = iter->second;

     bool curr_is_cushion = this->IsCushionTerm(curr_syst);
     if(curr_is_cushion) continue;

     double fractional_frac_err = uncert->OneSigmaErr(curr_syst); // fractional 1 sigma error
     double frac_scale = 1. + curr_twk_dial * fractional_frac_err;

     if(frac_scale < 0) {
       double curr_twk_dial_min = -1/fractional_frac_err;
       fSystValuesActual[curr_syst] = curr_twk_dial_min;
       LOG("ReW", pNOTICE) 
         << "Too large systematic variation for non-cushion systematic " << GSyst::AsString(curr_syst) 
         << " makes the corresponding fate scaling factor negative! Taking corrective action...";
       LOG("ReW", pINFO) 
         << "Tweak dial automatically re-adjusted from " << curr_twk_dial
         << " to " << curr_twk_dial_min;
     }
  }

  //
  // Loop over all non-cushion terms and calculate the change in their total fraction
  //

  double sum_nocushion_fate_fraction_nom = 0.;
  double sum_nocushion_fate_fraction_twk = 0.;

  for(iter = fSystValuesActual.begin(); iter != fSystValuesActual.end(); ++iter)
  {
     GSyst_t curr_syst     = iter->first;
     double  curr_twk_dial = iter->second;

     bool curr_is_cushion = this->IsCushionTerm(curr_syst);
     if(curr_is_cushion) continue;

     double fractional_frac_err = uncert->OneSigmaErr(curr_syst); // fractional 1 sigma error

     double frac_scale = 1. + curr_twk_dial * fractional_frac_err;
     double curr_frac  = genie::utils::rew::FateFraction(curr_syst, KE, frac_scale);
     double nom_frac   = genie::utils::rew::FateFraction(curr_syst, KE, 1.);

     sum_nocushion_fate_fraction_nom += nom_frac; 
     sum_nocushion_fate_fraction_twk += curr_frac;

     LOG("ReW", pDEBUG) 
         << "Non-cushion systematic " << GSyst::AsString(curr_syst) 
         << " (1 sigma frac. err = " << 100*fractional_frac_err
         << "%, tweak dial = " << curr_twk_dial 
         << "): fate fraction (KE = " << KE << " GeV) = "
         << nom_frac  << " (nominal), "
         << curr_frac << " (tweaked)";
  }

  LOG("ReW", pDEBUG) 
    << "Sum of non-cushion fate fractions = "
    << sum_nocushion_fate_fraction_nom << " (nominal), "
    << sum_nocushion_fate_fraction_twk << " (tweaked)";

  //
  // If the systematic variation specified by the user was so extreme that the
  // tweaked total fate fraction of non-cushion terms is above 1, automatically re-scale the
  // tweak dials so that the sum of their fate fraction is 1.
  // In that case, all cushion terms will be set to 0.
  //
  
  if(sum_nocushion_fate_fraction_twk > 1.) {
 
    LOG("ReW", pNOTICE) 
       << "Current sum of non-cushion fate fractions = "
       << sum_nocushion_fate_fraction_twk << " ( > 1 ). "
       << "Can not maintain unitarity with current set of systematic parameter values. "
       << "Taking corrective action...";

    for(iter = fSystValuesActual.begin(); iter != fSystValuesActual.end(); ++iter)
    {
       GSyst_t curr_syst     = iter->first;
       double  curr_twk_dial = iter->second;

       bool curr_is_cushion = this->IsCushionTerm(curr_syst);

       double fractional_frac_err = uncert->OneSigmaErr(curr_syst); // fractional 1 sigma error
  
       double frac_scale = 1. + curr_twk_dial * fractional_frac_err;
       double curr_frac  = genie::utils::rew::FateFraction(curr_syst, KE, frac_scale);

       double curr_frac_new     = (curr_is_cushion) ? 0 : curr_frac * (1./sum_nocushion_fate_fraction_twk);
       double frac_scale_new    = genie::utils::rew::WhichFateFractionScaleFactor(curr_syst, KE, curr_frac_new);
       double curr_twk_dial_new = (frac_scale_new != 0) ? (frac_scale_new - 1.)/fractional_frac_err : 0;

       fSystValuesActual[curr_syst] = curr_twk_dial_new;

       LOG("ReW", pINFO) 
         << ((curr_is_cushion) ? " Cushion " : " Non-cushion ")
         << "systematic " << GSyst::AsString(curr_syst) 
         << ": fate fraction re-adjusted from " << curr_frac << " to " << curr_frac_new
         << " (tweak dial re-adjusted from " << curr_twk_dial
         << " to " << curr_twk_dial_new << ")";
    }

    // update sum of fate fraction (set at limit)
    sum_nocushion_fate_fraction_twk = 1.;
  }

  // 
  // In the normal case where the sum of all non-cushion terms was less than 1,
  // leave them as they are. Adjust all all non-cushion terms accordingly so as to respect unitarity.
  // There are many ways to adjust the cushion terms, if more than one such term exists.
  // Chose to tweak them so that they all change the same amount, in units of the corresponding 1 sigma error.
  //

  else {

      double delta_fate_fraction = 
          sum_nocushion_fate_fraction_twk - sum_nocushion_fate_fraction_nom;

      LOG("ReW", pDEBUG) 
        << "Sum of non-cushion fate fractions changed by " << delta_fate_fraction
        << " - Change to be absorbed by the selected cushion terms...";

      double sum = 0;
      for(iter = fSystValuesActual.begin(); iter != fSystValuesActual.end(); ++iter)
      {
         GSyst_t curr_syst = iter->first;

         bool curr_is_cushion = this->IsCushionTerm(curr_syst);
         if(!curr_is_cushion) continue;

         double fractional_frac_err = uncert->OneSigmaErr(curr_syst); // fractional 1 sigma error
         double nom_frac = genie::utils::rew::FateFraction(curr_syst, KE, 1.);
         sum += (nom_frac * fractional_frac_err);
      }

      double twk_dial = -1.*delta_fate_fraction / sum;

      for(iter = fSystValuesActual.begin(); iter != fSystValuesActual.end(); ++iter)
      {
         GSyst_t curr_syst  = iter->first;
         bool curr_is_cushion = this->IsCushionTerm(curr_syst);
         if(!curr_is_cushion) continue;

         fSystValuesActual[curr_syst] = twk_dial;
      }

      LOG("ReW", pDEBUG) 
        << "To absorb change in the sum of non-cushion fate fractions, "
        << "the tweaking dial for all cushion terms needs to be set to "
        << twk_dial;

  } //check sum_nocushion_fate_fraction_twk < 1.

  //
  // Confirm unitarity
  //

  LOG("ReW", pINFO) << "Confirming unitarity for current fate fractions";

  double sum_fate_fraction_all = 0;
  for(iter = fSystValuesActual.begin(); iter != fSystValuesActual.end(); ++iter)
  {
     GSyst_t curr_syst     = iter->first;
     double  curr_twk_dial = iter->second;

     bool curr_is_cushion = this->IsCushionTerm(curr_syst);

     double fractional_frac_err = uncert->OneSigmaErr(curr_syst); // fractional 1 sigma error

     double frac_scale = 1. + curr_twk_dial * fractional_frac_err;
     double curr_frac  = genie::utils::rew::FateFraction(curr_syst, KE, frac_scale);

     LOG("ReW", pINFO) 
         << "Systematic " << GSyst::AsString(curr_syst) 
         << " (1 sigma frac. err = " << 100*fractional_frac_err
         << "%, tweak dial = " << curr_twk_dial 
         << "): Current fate fraction (KE = " << KE << " GeV) = " << curr_frac
         << ((curr_is_cushion) ? " ** cushion term ** " : "");

     sum_fate_fraction_all += curr_frac;
  }

  LOG("ReW", pINFO) << "Current sum of all fate fractions = " << sum_fate_fraction_all;

  if(TMath::Abs(sum_fate_fraction_all-1) > 0.01) {   
      LOG("ReW", pWARN) << "Unitarity violation level exceeded 1 part in 100.";
      LOG("ReW", pWARN) << "Current sum of all fate fractions = " << sum_fate_fraction_all;
  }

  map<GSyst_t, double>::const_iterator dial_iter = fSystValuesActual.find(syst);
  return dial_iter->second;
}
//___________________________________________________________________________
double GReWeightINukeParams::Fates::ChisqPenalty(void) const
{
// calculate an average chisq over the energy range

  double chisq = 0.0;

  const int    kN      = 200;
  const double kKEmin  =  0;
  const double kKEmax  = 10;
  const double kKEstep = (kKEmax-kKEmin)/(kN-1);

  GSyst_t syst = kNullSystematic;
  int i=0;
  while( (syst = (fHadType == kRwINukePion) ?
            GSyst::NextPionFateSystematic(i++) :
            GSyst::NextNuclFateSystematic(i++)
         ) != kNullSystematic
       )
  {
    for(int j=0; j < kN; j++) {
       double KE = kKEmin + j * kKEstep;
       double twkdial = this->ActualTwkDial(syst,KE);
       chisq += (TMath::Power(twkdial, 2.));
    }//j
  }//i 

  chisq /= kN;

  return chisq;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsIncluded(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fSystValuesUser.find(syst);
  if(iter != fSystValuesUser.end()) return true;
  return false;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsCushionTerm(GSyst_t syst) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsCushion.find(syst);
  if(iter != fIsCushion.end()) return iter->second;

  LOG("ReW", pERROR) 
    << "Cannot query 'is cushion' flag for systematic " 
    << GSyst::AsString(syst);

  return false;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsTweaked(GSyst_t syst) const
{
  if(!this->IsHandled (syst)) return false;
  if(!this->IsIncluded(syst)) return false;

  map<GSyst_t, double>::const_iterator dial_iter = fSystValuesUser.find(syst);
  double dial = dial_iter->second;

  bool tweaked = (TMath::Abs(dial) >= controls::kASmallNum);

  LOG("ReW", pDEBUG) 
    << "Systematic " << GSyst::AsString(syst)
    << ((tweaked) ? " was " : " was not ")
    << "tweaked (tweaking dial = " << dial << ")";

  return tweaked;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsTweaked(void) const
{
  GSyst_t syst = kNullSystematic;
  int i=0;
  while( (syst = (fHadType == kRwINukePion) ?
            GSyst::NextPionFateSystematic(i++) :
            GSyst::NextNuclFateSystematic(i++)
         ) != kNullSystematic
       )
  {
    bool tweaked = this->IsTweaked(syst);
    if(tweaked) return true;
  }
  return false;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsHandled(GSyst_t syst) const
{
//
//
  if(fHadType == kRwINukePion) {
    if(GSyst::IsINukePionFateSystematic(syst)) return true;
  }
  else 
  if(fHadType == kRwINukeNucl) {
    if(GSyst::IsINukeNuclFateSystematic(syst)) return true;
  }

  LOG("ReW", pWARN) 
    << "Can not handle systematic " << GSyst::AsString(syst);

  return false;
}
//___________________________________________________________________________
void GReWeightINukeParams::Fates::AddCushionTerms(void)
{
// When this method is called for the first time, all terms not already 
// included are included as cushion terms.
// On every call, the values of the cushion terms are reset.

  LOG("ReW", pINFO) << "Adding cushion terms...";

  int ncushions = 0;

  GSyst_t syst = kNullSystematic;
  int i=0;
  while( (syst = (fHadType == kRwINukePion) ?
            GSyst::NextPionFateSystematic(i++) :
            GSyst::NextNuclFateSystematic(i++)
         ) != kNullSystematic
       )
  {
      if(this->IsIncluded(syst)) {
        bool is_cushion = this->IsCushionTerm(syst);
	LOG("ReW", pINFO) 
          << "Systematic " << GSyst::AsString(syst) 
          << " was already specified as a"
          << ((is_cushion) ? " cuhsion " : " non-cushion ")
          << "term";
        if(is_cushion) {
           fSystValuesUser[syst]  = 0.; 
           ncushions++;
        }
      }
      else {
	LOG("ReW", pINFO) 
          << "Adding systematic " << GSyst::AsString(syst) 
          << " as a cushion term.";

        fSystValuesUser[syst]  = 0.; 
        fIsCushion[syst] = true;
        ncushions++;
      }
  } // end loop over relevant fate systematics

  LOG("ReW", pINFO) << "Number of cushion terms = " << ncushions;

  if(ncushions <= 0) {
     LOG("ReW", pFATAL) 
       << "There must be at least one cushion term (" 
       << ncushions << " were set)";
     exit(1);
  }
}
//___________________________________________________________________________
//
// MFP nested class 
//
//___________________________________________________________________________ 
GReWeightINukeParams::MFP::MFP(GReWeightINukeParams::HadronType_t ht)
{
  assert(ht != kRwINukeUndefined);

  fHadType = ht;

  // get corresponding GSyst_t param
  if      (ht == kRwINukePion) fSyst = kINukeTwkDial_MFP_pi;
  else if (ht == kRwINukeNucl) fSyst = kINukeTwkDial_MFP_N;
  else {
    LOG("ReW", pFATAL) 
      << "Mean-free path reweighting code handles only pions or nucleons";
    exit(1);
  }

  this->Reset();
}
//___________________________________________________________________________ 
GReWeightINukeParams::MFP::~MFP()
{

}
//___________________________________________________________________________ 
double GReWeightINukeParams::MFP::ScaleFactor(void) const
{
  if(! this->IsTweaked()) return 1.;

  GSystUncertainty * uncert = GSystUncertainty::Instance();

  double mfp_sigma    = uncert->OneSigmaErr(fSyst);
  double mfp_twkdial  = this->TwkDial();
  double scale_factor = 1.0 + mfp_sigma * mfp_twkdial;

  return scale_factor;  
}
//___________________________________________________________________________ 
double GReWeightINukeParams::MFP::TwkDial(void) const
{
  return fTwkDial;
}
//___________________________________________________________________________ 
bool GReWeightINukeParams::MFP::IsIncluded(void) const
{
  return fIsIncluded;
}
//___________________________________________________________________________ 
bool GReWeightINukeParams::MFP::IsTweaked(void) const
{
  return (TMath::Abs(fTwkDial) >= controls::kASmallNum);
}
//___________________________________________________________________________ 
double GReWeightINukeParams::MFP::ChisqPenalty(void) const
{
  return TMath::Power(fTwkDial, 2.0);
}
//___________________________________________________________________________ 
void GReWeightINukeParams::MFP::Reset(void) 
{
  fTwkDial    = 0.0;
  fIsIncluded = false;
}
//___________________________________________________________________________ 
void GReWeightINukeParams::MFP::SetTwkDial(double val) 
{
  fTwkDial    = val;
  fIsIncluded = true;
}
//___________________________________________________________________________ 

