//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
         Imperial College London

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

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
*/
//____________________________________________________________________________

#include <cassert>

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
  fParmPionFates -> Reconfigure();
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
void GReWeightINukeParams::SetCurTwkDial(GSyst_t syst, double val)
{
  if(GSyst::IsINukePionFateSystematic(syst)) 
  {
    fParmPionFates->SetCurTwkDial(syst,val);
  }
  else
  if(GSyst::IsINukeNuclFateSystematic(syst))
  {
    fParmNuclFates->SetCurTwkDial(syst,val);
  }
  else
  if(GSyst::IsINukePionMeanFreePathSystematic(syst))
  {
    fParmPionMFP->SetCurTwkDial(val);
  }
  else
  if(GSyst::IsINukeNuclMeanFreePathSystematic(syst))
  {
    fParmNuclMFP->SetCurTwkDial(val);
  } 
}
//___________________________________________________________________________ 
//
// Fates nested class 
//
//___________________________________________________________________________ 
GReWeightINukeParams::Fates::Fates(GReWeightINukeParams::HadronType_t ht) 
{
  assert(ht != kRwINukeUndefined);

  fHadType       = ht;
  fNSysts        = 5;
  fNCushionTerms = 1;
  fCushAvgChisq  = 0;

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
  fSystListMap.clear();
  fIsCushionMap.clear();
}
//___________________________________________________________________________
void GReWeightINukeParams::Fates::Reconfigure()
{
  this->AddCushionTerms();
  this->CheckUnitarity();
}
//___________________________________________________________________________
void GReWeightINukeParams::Fates::SetCurTwkDial(GSyst_t syst, double val)
{
  // check type of systematic 
  if(!this->IsHandled(syst)) return;

  // update tweaking dial
  fSystListMap[syst]  = val;
  fIsCushionMap[syst] = false;
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

  double fractional_error  = uncert->OneSigmaErr(syst);
  double twk_dial          = this->CurTwkDial(syst, KE);

  double fate_fraction_scale = 1. + twk_dial * fractional_error;

  return fate_fraction_scale;
}
//___________________________________________________________________________
double GReWeightINukeParams::Fates::CurTwkDial(GSyst_t syst, double KE) const
{
  if(! this->IsIncluded(syst)) 
  {
    LOG("ReW", pWARN) 
      << "Systematic " << GSyst::AsString(syst) << " not included";
    return 0.;
  }
  bool is_cushion = this->IsCushionTerm(syst);

  // If it not acting as a cushion term just return the tweaking dial 
  // value that was set. 
  if(!is_cushion){
    map<GSyst_t, double>::const_iterator dial_iter = fSystListMap.find(syst);
    return dial_iter->second;
  }

  // For cushion terms calculate the required treak dial value that will
  // maintain unitarity.
  // JIMTODO - need to properly reference/explain this formulae

  if(KE < 0.) {
    LOG("ReW", pWARN)  << "The cushion term depends on kinetic energy!";
    LOG("ReW", pERROR) << "Unspecified kinetic energy!";
    return 0.;
  }

  GSystUncertainty * uncert = GSystUncertainty::Instance();

  double total_fraction_nocushion = 0.;
  double sum_non_cushion_fractions = 0.;
  double sum_error_non_cushion_fractions =0.;

  map<GSyst_t, double>::const_iterator iter = fSystListMap.begin();  
  for( ; iter != fSystListMap.end(); ++iter)
  {
     GSyst_t curr_syst       = iter->first;
     double  curr_twk_dial   = iter->second;
     bool    curr_is_cushion = this->IsCushionTerm(curr_syst);
     if(!curr_is_cushion){
       double fractional_error    = uncert->OneSigmaErr(curr_syst);
       double fate_fraction_scale = 1. + curr_twk_dial * fractional_error;

       total_fraction_nocushion += 
          genie::utils::rew::FateFraction(curr_syst, KE, fate_fraction_scale);
     }// no cushion
     else {
       double fractional_error    = uncert->OneSigmaErr(curr_syst);
       double fate_fraction_scale = 1.; //+ curr_twk_dial * fractional_error;
       double default_fate_fraction = genie::utils::rew::FateFraction(curr_syst, KE, fate_fraction_scale);
       double error_on_default_fate_fraction = fractional_error*default_fate_fraction;
       sum_non_cushion_fractions += default_fate_fraction;
       sum_error_non_cushion_fractions += error_on_default_fate_fraction;
     } // is cushion
  } // systs loop

  double twk_dial = ( (1.0- total_fraction_nocushion) - sum_non_cushion_fractions)/sum_error_non_cushion_fractions;

  // This is a temporary measure to ensure that we are not tweaking any of the parameters in such a way
  // that would require a negative cross section. This issue needs to be resolved.
  assert( (twk_dial*uncert->OneSigmaErr(syst)) > -1.0 );

  return twk_dial;
}
//___________________________________________________________________________
double GReWeightINukeParams::Fates::ChisqPenalty(void) const
{
  double chisq = 0.0;

  // loop over the non cushion terms & compute the chisq_{penalty}
  map<GSyst_t, bool>::const_iterator iter = fIsCushionMap.begin();
  for( ; iter != fIsCushionMap.end(); ++iter) 
  {
      GSyst_t syst = iter->first;
      if(this->IsCushionTerm(syst)) continue;

      double dial = this->CurTwkDial(syst);
      chisq += (TMath::Power(dial, 2.0));
  } 

  // add the contribution from the cushion term, averaged over kinetic energy
  chisq += fCushAvgChisq;

  return chisq;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsIncluded(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fSystListMap.find(syst);
  if(iter != fSystListMap.end()) return true;
  return false;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsCushionTerm(GSyst_t syst) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsCushionMap.find(syst);
  if(iter != fIsCushionMap.end()) return iter->second;

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

  map<GSyst_t, double>::const_iterator dial_iter = fSystListMap.find(syst);
  double dial = dial_iter->second;

  bool tweaked = (TMath::Abs(dial) >= controls::kASmallNum);
  return tweaked;
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::IsTweaked(void) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsCushionMap.begin();
  for ( ; iter != fIsCushionMap.end(); ++iter) {
    GSyst_t syst = iter->first;
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
// if have added all non cushion terms then fill rest of map with the 
// remaining cushion terms.

  //int ncursyst = (int)fSystListMap.size();
  bool add_cushion_terms = true;//(ncursyst == (fNSysts-4/*fNCushionTerms*/));

  if (!add_cushion_terms) {
     LOG("ReW", pDEBUG) << "Cannot determine the cushion terms just yet";
     return;
  }

  GSyst_t syst = kNullSystematic;
  int i=0;
  while( 
         (syst = (fHadType == kRwINukePion) ?
            GSyst::NextPionFateSystematic(i++) :
            GSyst::NextNuclFateSystematic(i++)
         ) != kNullSystematic
       )
  {      
      if(this->IsIncluded(syst)) {
	LOG("ReW", pDEBUG) 
          << "Systematic " << GSyst::AsString(syst) << " already included";
      }
      else {
	LOG("ReW", pDEBUG) 
          << "Adding "<< GSyst::AsString(syst) << " as cushion term.";

        fSystListMap[syst]  = 0.; /// JIMTODO ?
        fIsCushionMap[syst] = true;
      }
  } // end loop over relevant fate systematics
}
//___________________________________________________________________________
bool GReWeightINukeParams::Fates::CheckUnitarity(int n_points)
{
  this->AddCushionTerms();

  fCushAvgChisq = 0.0;

  double KEmin  = INukeHadroData::fMinKinEnergy;   // min kinetic energy
  double KEmax  = INukeHadroData::fMaxKinEnergyHA; // max kinetic energy
  double KEstep = (KEmax-KEmin)/(n_points-1.);

  // on the way work out the average chi-square for the cushion terms
  double total_chisq = 0.0;

  // loop over the hadro-data energy range and check that it is possible
  // to maintain unitarity using the cushion terms over this range.

  double KE = KEmin;
  for(int i = 0 ; i < n_points; i++) {

    // total fraction of non-cushion term
    double total_fraction_nocushion = 0;
  
    // loop over the fates 
    map<GSyst_t, double>::const_iterator iter = fSystListMap.begin();  
    for( ; iter != fSystListMap.end(); ++iter)
    {
	GSyst_t syst     = iter->first;
	double  twk_dial = iter->second;

	bool is_cushion = this->IsCushionTerm(syst);

        if(!is_cushion) {
          double fractional_error    = GSystUncertainty::Instance()->OneSigmaErr(syst);
          double fate_fraction_scale = 1 + twk_dial * fractional_error;

          total_fraction_nocushion += 
             genie::utils::rew::FateFraction(syst, KE, fate_fraction_scale);
        } // is cushion
        else {
          total_chisq += TMath::Power(twk_dial, 2.0);
        }

    } //syst loop

    // now check that can maintain unitarity for this energy 
    bool can_maintain_unitarity = (total_fraction_nocushion < 1.); 
    if (!can_maintain_unitarity)
    {
      LOG("ReW", pNOTICE) 
        << "Cannot maintain unitarity for current set of fate params. "
	<< "The test failed for hadron KE = " << KE;
      return false;
    }

    KE += KEstep;

  } // end kinetic energy loop

  // now set the current average
  fCushAvgChisq = total_chisq / ((double) n_points);

  return true;
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

  // corresponding GSyst_t
  if(ht == kRwINukePion) fSyst = kINukeTwkDial_MFP_pi;
  if(ht == kRwINukeNucl) fSyst = kINukeTwkDial_MFP_N;

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
  double mfp_twkdial  = this->CurTwkDial();
  double scale_factor = 1.0 + mfp_sigma * mfp_twkdial;

  return scale_factor;  
}
//___________________________________________________________________________ 
double GReWeightINukeParams::MFP::CurTwkDial(void) const
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
void GReWeightINukeParams::MFP::SetCurTwkDial(double val) 
{
  fTwkDial    = val;
  fIsIncluded = true;
}
//___________________________________________________________________________ 

