//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
         Imperial College London

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________
#include <iostream>

#include "ReWeight/GReWeightINukeFateParams.h"
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
using std::cout;
using std::endl;

//___________________________________________________________________________
GReWeightINukeFateParams::GReWeightINukeFateParams(int n_systs, int n_cushions)
{
  cout << "----- Constructing GReWeightINukeFateParams::GReWeightINukeFateParams("<< 
    n_systs <<","<< n_cushions <<") ! -----" << endl;
  this->Initialize(n_systs, n_cushions);
}
//___________________________________________________________________________
GReWeightINukeFateParams::~GReWeightINukeFateParams(void)
{

}
//___________________________________________________________________________
void GReWeightINukeFateParams::Initialize(int n_systs=5, int n_cushions=1)
{
  fNSysts = n_systs;
  fNCushionTerms = n_cushions;
  fArePionSysts = false;
  fAreNuclSysts = false;
}
//___________________________________________________________________________
void GReWeightINukeFateParams::Reset()
{
  fSystListMap.clear();
  fIsCushionMap.clear();
//  fIsCushionTermList.clear();
//  fSystTwkDials.clear();
  this->Initialize();

  return;
}
//___________________________________________________________________________
void GReWeightINukeFateParams::Reconfigure()
{

  fAverageChi2 = 0.0;

  this->AddCushionTerms();
  this->CheckFatesUnity();

  return;
}
//___________________________________________________________________________
void GReWeightINukeFateParams::SetSystematic(EGSyst syst, double val)
{
    this->SetFateSyst(syst, false, val);
}
//___________________________________________________________________________
void GReWeightINukeFateParams::SetFateSyst(EGSyst syst, bool is_cushion, double val)
{
  if(this->CheckSystType(syst)){
    fSystListMap[syst]=val;
    fIsCushionMap[syst]=is_cushion;
  }
  return;
}
//___________________________________________________________________________
bool GReWeightINukeFateParams::CheckSystType(EGSyst syst)
{
  bool is_inuke = GSyst::Type(syst) == kSystType_INuke; 
  if(!is_inuke){ 
    std::cout << "Cannot add this systematic to the " <<
      "inuke helper as they are not inuke parameters" << std::endl;
    return false;
  }

  // do not set if not a fate systematic
  if(this->GetHadronType(syst) == kNull){
    std::cout << "Cannot add this systematic to the inuke " <<
      "inuke helper as is not a pion or nucl fate parameter" << std::endl;
    return false;
  }

  // if hadron type is unset set it for the first time 
  if(fHadronType == kNull){ fHadronType = this->GetHadronType(syst);}

  // allow only type of pion or nucleon parameters per GReWeightINukeFateParams
  if(this->GetHadronType(syst) == fHadronType){
    return true;
  }
  
  std::cout << "Cannot add this systematic to the inuke " <<
    "helper as it controls the fate for a different hadron "<<
    "type to those already added to the helper"<< std::endl;  

  return false;
  
}
//___________________________________________________________________________
double GReWeightINukeFateParams::GetFateXSecTwk(EGSyst syst, double kinE, double twk_dial_val)
{
  double fractional_error = GSystUncertainty::Instance()->OneSigmaErr(syst);
  double fate_xsec_tweak = 
    (1.0 + twk_dial_val * fractional_error * genie::utils::rew::FateXSec(syst, kinE) );
	
  return fate_xsec_tweak;
}
//___________________________________________________________________________
double GReWeightINukeFateParams::TwkDialValue(EGSyst syst, double kinE)
{
  double twk_dial_val = 0.0;

  fSystListMap_iter  = fSystListMap.find(syst);
  fIsCushionMap_iter = fIsCushionMap.find(syst);

  if(fSystListMap_iter == fSystListMap.end() || fIsCushionMap_iter == fIsCushionMap.end()){ 
    //JIMTODO - need to add warning here and do something else than just setting twk_dial_val
    //return 0.0;
    std::cout << "Cannot get twkdial value for systematic as it cannot "<<
      "be found in the inuke helper internal maps" << std::endl;
  }

  bool is_cushion = fIsCushionMap_iter->second;

  ///< If it not acting as a cushion term just return the tweaking dial 
  ///< value that was set. 
  if(!is_cushion){
    return fSystListMap_iter->second;
  }
  ///< For cushion terms calculate the required treak dial value that will
  ///< maintain unitarity.
  else {

    ///< Now work out values needed to calculate cushion tweak value.
    ///< JIMTODO - need to properly reference/explain this formulae
    double numerator_val   = 0.0;
    double denominator_val = 0.0;
    for(fIsCushionMap_iter  = fIsCushionMap.begin();  
	fIsCushionMap_iter != fIsCushionMap.end();
	fIsCushionMap_iter ++)
      {
	bool loop_is_cushion = fIsCushionMap_iter->second;
	EGSyst loop_curr_syst = fIsCushionMap_iter->first;

	double fractional_error = GSystUncertainty::Instance()->OneSigmaErr(loop_curr_syst);
	///< contribbution to sum over syst_frac_error(ke)*syst_xsec(ke)  
	if(loop_is_cushion){
	  denominator_val += (fractional_error *
			      genie::utils::rew::FateXSec(loop_curr_syst, kinE) );
	}
	///< contribution to sum over syst_twk_val*syst_frac_error(ke)*syst_xsec(ke)  
	else{ 
	 numerator_val -=  (fSystListMap.find(loop_curr_syst)->second *
			    fractional_error *
			    genie::utils::rew::FateXSec(loop_curr_syst, kinE) );
	}
      } ///< end loop over systs


    ///< JIMTODO - need to correctly handle case when denominator is zero
    twk_dial_val = numerator_val/denominator_val; 

  } ///< end if(cushion term)
  return twk_dial_val;
}
//___________________________________________________________________________
EHadronType_t GReWeightINukeFateParams::GetHadronType(EGSyst syst)
{
  EHadronType_t hadron_type = kNull;

  if(syst == kSystINuke_CExTwk_pi  ||
     syst == kSystINuke_ElTwk_pi   ||
     syst == kSystINuke_InelTwk_pi ||
     syst == kSystINuke_AbsTwk_pi  ||
     syst == kSystINuke_PiProdTwk_pi){fHadronType = kPionInNucleus;}

  else if(syst == kSystINuke_CExTwk_N  ||
  	  syst == kSystINuke_ElTwk_N   ||
  	  syst == kSystINuke_InelTwk_N ||
  	  syst == kSystINuke_AbsTwk_N  ||
  	  syst == kSystINuke_PiProdTwk_N){fHadronType = kNuclInNucleus;}

  return hadron_type;
}
//___________________________________________________________________________
void GReWeightINukeFateParams::AddCushionTerms()
{
  // if have added all non cushion terms then fill rest
  // of map with the remaining cushion terms.
  if(fSystListMap.size() == fNSysts - fNCushionTerms){

    const int n_fates = 5;
    EGSyst all_nucl_syst[n_fates] = { kSystINuke_CExTwk_N   ,
                                      kSystINuke_ElTwk_N    ,
				      kSystINuke_InelTwk_N  ,
				      kSystINuke_AbsTwk_N   ,
				      kSystINuke_PiProdTwk_N };

    EGSyst all_pion_syst[n_fates] = { kSystINuke_CExTwk_pi   ,
                                      kSystINuke_ElTwk_pi    ,
				      kSystINuke_InelTwk_pi  ,
				      kSystINuke_AbsTwk_pi   ,
				      kSystINuke_PiProdTwk_pi };

    // now add as a cushion term if not already added 
    for(int ith_fate = 0; ith_fate < n_fates ; ith_fate ++ ){
      
      EGSyst syst_to_add = kSystNull;
      if(fArePionSysts){ syst_to_add = all_pion_syst[ith_fate];}
      if(fAreNuclSysts){ syst_to_add = all_nucl_syst[ith_fate];}

      fIsCushionMap_iter = fIsCushionMap.find(syst_to_add);
      bool is_included = fIsCushionMap_iter != fIsCushionMap.end();
  
      if(is_included){
	std::cout << "Fate "<< GSyst::AsString(syst_to_add) << " is already set as a tweaking." << std::endl;
      }
      else {
	std::cout << "Adding "<< GSyst::AsString(syst_to_add) << " as a cushion term." << std::endl;
	this->SetFateSyst(syst_to_add, true, 0.0);
      }
    } // end loop over fates
  } // end fill in cushion terms
  else if(fSystListMap.size() < (fNSysts - fNCushionTerms)){
    std::cout << "Cannot determine cushion terms as have" <<
      " not filled in all the non-cushion terms." << std::endl;
  }

  return;
}
//___________________________________________________________________________
bool GReWeightINukeFateParams::CheckFatesUnity(int n_points)
{
  this->AddCushionTerms();

  double KEmin = INukeHadroData::fMinKinEnergy;
  double KEmax = INukeHadroData::fMaxKinEnergyHA;

  double KE_hadron = KEmin;
  double E_step = (double) (KEmax-KEmin)/n_points;

  // work out the average chi-squared
  double total_chisq = 0.0;

  // loop over the hadro-data energy range and check that it is possible
  // to maintain unitarity using the cushion terms over this range.
  for(int i = 0 ; i < n_points ; i++){
    KE_hadron += E_step;
    // variables to store the sum of all the fates and also the sum of the tweaked 
    // non cushion terms.
    double total_xsec = 0.0;
    double total_tweaked_xsec_non_cushion = 0.0;
  
    // loop over the fates 
    for(fIsCushionMap_iter  = fIsCushionMap.begin();  
	fIsCushionMap_iter != fIsCushionMap.end();
	fIsCushionMap_iter ++)
      {
	bool loop_is_cushion = fIsCushionMap_iter->second;
	EGSyst loop_curr_syst = fIsCushionMap_iter->first;

	total_xsec += this->GetFateXSecTwk(loop_curr_syst, KE_hadron, 0.0);
	if(loop_is_cushion==false){ 
	  double curr_tweak_val = fSystListMap.find(loop_curr_syst)->second;
	  total_tweaked_xsec_non_cushion += this->GetFateXSecTwk(loop_curr_syst, KE_hadron, curr_tweak_val);
	  total_chisq += this->GetFatesChi2(KE_hadron);
	}
      } ///< end loop over systs

    // now check that can maintain unitarity for this energy 
    bool can_maintain_unity = total_tweaked_xsec_non_cushion < total_xsec; 
    if(can_maintain_unity == false){
      std::cout << "Warning cannot maintain unitarity for this particular "<<
	" set of fate parameters and tweak values. The test failed for hadron ke = "<< KE_hadron << "!" << std::endl;
      return false;
    }

    // now set the current average
    fAverageChi2 = (double) total_chisq / n_points;

  } // end energy loop
  
  return true;
}
//___________________________________________________________________________
double GReWeightINukeFateParams::GetFatesChi2(double kinE)
{
  double chi2_val = 0.0;

  // loop over the fates 
  for(fIsCushionMap_iter  = fIsCushionMap.begin();  
      fIsCushionMap_iter != fIsCushionMap.end();
      fIsCushionMap_iter ++)
    {
      EGSyst loop_curr_syst = fIsCushionMap_iter->first;
      double twk_dial_val = this->TwkDialValue(loop_curr_syst , kinE);
      chi2_val += (twk_dial_val*twk_dial_val);
    } ///< end loop over systs

  return chi2_val;
}
//___________________________________________________________________________
EGSyst GReWeightINukeFateParams::GetRelevantSyst(EINukeFateHA_t hadron_fate)
{

  if(fArePionSysts){ 
    switch(hadron_fate){
     case (kIHAFtCEx    ) : return kSystINuke_CExTwk_pi   ; break;
     case (kIHAFtElas   ) : return kSystINuke_ElTwk_pi    ; break;
     case (kIHAFtInelas ) : return kSystINuke_InelTwk_pi  ; break;
     case (kIHAFtAbsNP  ) : return kSystINuke_AbsTwk_pi   ; break;
     case (kIHAFtAbsPP  ) : return kSystINuke_AbsTwk_pi   ; break;
     case (kIHAFtAbsNPP ) : return kSystINuke_AbsTwk_pi   ; break;
     case (kIHAFtAbsNNP ) : return kSystINuke_AbsTwk_pi   ; break;
     case (kIHAFtAbs2N2P) : return kSystINuke_AbsTwk_pi   ; break;
     case (kIHAFtAbs2N3P) : return kSystINuke_AbsTwk_pi   ; break;
     case (kIHAFtNPip   ) : return kSystINuke_PiProdTwk_pi; break;
     case (kIHAFtNPipPi0) : return kSystINuke_PiProdTwk_pi; break;
     default: return kSystNull;
     } 
  }// endif pions
  else if(fAreNuclSysts){
    switch(hadron_fate){
     case (kIHAFtCEx    ) : return kSystINuke_CExTwk_N   ; break;
     case (kIHAFtElas   ) : return kSystINuke_ElTwk_N    ; break;
     case (kIHAFtInelas ) : return kSystINuke_InelTwk_N  ; break;
     case (kIHAFtAbsNP  ) : return kSystINuke_AbsTwk_N   ; break;
     case (kIHAFtAbsPP  ) : return kSystINuke_AbsTwk_N   ; break;
     case (kIHAFtAbsNPP ) : return kSystINuke_AbsTwk_N   ; break;
     case (kIHAFtAbsNNP ) : return kSystINuke_AbsTwk_N   ; break;
     case (kIHAFtAbs2N2P) : return kSystINuke_AbsTwk_N   ; break;
     case (kIHAFtAbs2N3P) : return kSystINuke_AbsTwk_N   ; break;
     case (kIHAFtNPip   ) : return kSystINuke_PiProdTwk_N; break;
     case (kIHAFtNPipPi0) : return kSystINuke_PiProdTwk_N; break;
     default: return kSystNull;
     }  
  }// endif nucleons
  return kSystNull;  
}  
//___________________________________________________________________________




