//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Controls.h"
#include "Registry/Registry.h"
#include "ReWeight/GReWeightNuXSecParams.h"
#include "ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//____________________________________________________________________________
GReWeightNuXSecParams::GReWeightNuXSecParams()
{
  this->Init();
}
//____________________________________________________________________________
GReWeightNuXSecParams::~GReWeightNuXSecParams()
{

}
//____________________________________________________________________________
double GReWeightNuXSecParams::DefValue(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fDefParams.find(syst);
  if(iter != fDefParams.end()) return iter->second;
  else return 0;
}
//____________________________________________________________________________
double GReWeightNuXSecParams::CurValue(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fCurParams.find(syst);
  if(iter != fCurParams.end()) return iter->second;
  else return 0;
}
//____________________________________________________________________________
double GReWeightNuXSecParams::CurTwkDial(GSyst_t syst) const
{
  map<GSyst_t, double>::const_iterator iter = fCurTwkDial.find(syst);
  if(iter != fCurTwkDial.end()) return iter->second;
  else return 0;
}
//____________________________________________________________________________
bool GReWeightNuXSecParams::IsIncluded(GSyst_t syst) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsTweaked.find(syst);
  if(iter != fIsTweaked.end()) return true;
  else return false;
}
//____________________________________________________________________________
bool GReWeightNuXSecParams::IsTweaked(GSyst_t syst) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsTweaked.find(syst);
  if(iter != fIsTweaked.end()) return iter->second;
  else return false;
}
//____________________________________________________________________________
bool GReWeightNuXSecParams::IsTweaked(void) const
{
  map<GSyst_t, bool>::const_iterator iter = fIsTweaked.begin();
  for( ; iter != fIsTweaked.end(); ++iter)
  {
    if(iter->second) return true;
  }
  return false;
}
//____________________________________________________________________________
double GReWeightNuXSecParams::ChisqPenalty(void) const
{
 double chisq = 0;

  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_MaQEL),         2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_MvQEL),         2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_MaRES),         2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_MvRES),         2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvpCC1pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvpCC2pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvpNC1pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvpNC2pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvnCC1pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvnCC2pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvnNC1pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvnNC2pi),      2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarpCC1pi),   2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarpCC2pi),   2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarpNC1pi),   2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarpNC2pi),   2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarnCC1pi),   2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarnCC2pi),   2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarnNC1pi),   2.0);
  chisq += TMath::Power(this->CurTwkDial(kSystNuXSec_RvbarnNC2pi),   2.0);

  return chisq;
}
//____________________________________________________________________________
void GReWeightNuXSecParams::Reconfigure(void) const
{
  if (! this->IsTweaked() ) return;

  // Use current reweighting params inputs to tweak the GENIE configuration
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * user_config = conf_pool->GlobalParameterList();

  if(this->IsIncluded(kSystNuXSec_MaQEL)) {
      user_config->Set("QEL-Ma", this->CurValue(kSystNuXSec_MaQEL));
  }
  if(this->IsIncluded(kSystNuXSec_MvQEL)) {
      user_config->Set("QEL-Mv", this->CurValue(kSystNuXSec_MvQEL));
  }
  if(this->IsIncluded(kSystNuXSec_MaRES)) {
      user_config->Set("RES-Ma", this->CurValue(kSystNuXSec_MaRES));
  }
  if(this->IsIncluded(kSystNuXSec_MvRES)) {
      user_config->Set("RES-Mv", this->CurValue(kSystNuXSec_MvRES));
  }

  //
  // note: *** need to add resonance background params ***
  //


  // Force reconfiguration of GENIE algorithms
  AlgFactory * algf = AlgFactory::Instance();
  algf->ForceReconfiguration();
}
//____________________________________________________________________________
void GReWeightNuXSecParams::LoadDefaults(void)
{
  // Get the default cross section parameters 
  // (supported by the current version of the ReWeight package)
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * user_config = conf_pool->GlobalParameterList();
  user_config->UnLock();

  double MaQEL_def       = user_config->GetDouble("QEL-Ma");
  double MvQEL_def       = user_config->GetDouble("QEL-Mv");
  double MaRES_def       = user_config->GetDouble("RES-Ma");
  double MvRES_def       = user_config->GetDouble("RES-Mv");
  double RvpCC1pi_def    = user_config->GetDouble("DIS-HMultWgt-vp-CC-m2" );
  double RvpCC2pi_def    = user_config->GetDouble("DIS-HMultWgt-vp-CC-m3" );
  double RvpNC1pi_def    = user_config->GetDouble("DIS-HMultWgt-vp-NC-m2" );
  double RvpNC2pi_def    = user_config->GetDouble("DIS-HMultWgt-vp-NC-m3" );
  double RvnCC1pi_def    = user_config->GetDouble("DIS-HMultWgt-vn-CC-m2" );
  double RvnCC2pi_def    = user_config->GetDouble("DIS-HMultWgt-vn-CC-m3" );
  double RvnNC1pi_def    = user_config->GetDouble("DIS-HMultWgt-vn-NC-m2" );
  double RvnNC2pi_def    = user_config->GetDouble("DIS-HMultWgt-vn-NC-m3" );
  double RvbarpCC1pi_def = user_config->GetDouble("DIS-HMultWgt-vbp-CC-m2");
  double RvbarpCC2pi_def = user_config->GetDouble("DIS-HMultWgt-vbp-CC-m3");
  double RvbarpNC1pi_def = user_config->GetDouble("DIS-HMultWgt-vbp-NC-m2");
  double RvbarpNC2pi_def = user_config->GetDouble("DIS-HMultWgt-vbp-NC-m3");
  double RvbarnCC1pi_def = user_config->GetDouble("DIS-HMultWgt-vbn-CC-m2");
  double RvbarnCC2pi_def = user_config->GetDouble("DIS-HMultWgt-vbn-CC-m3");
  double RvbarnNC1pi_def = user_config->GetDouble("DIS-HMultWgt-vbn-NC-m2");
  double RvbarnNC2pi_def = user_config->GetDouble("DIS-HMultWgt-vbn-NC-m3");

  // store local copies of the defaults
  this->SetDefValue(kSystNuXSec_MaQEL,       MaQEL_def       );
  this->SetDefValue(kSystNuXSec_MvQEL,       MvQEL_def       );
  this->SetDefValue(kSystNuXSec_MaRES,       MaRES_def       );
  this->SetDefValue(kSystNuXSec_MvRES,       MvRES_def       );
  this->SetDefValue(kSystNuXSec_RvpCC1pi,    RvpCC1pi_def    );
  this->SetDefValue(kSystNuXSec_RvpCC2pi,    RvpCC2pi_def    );
  this->SetDefValue(kSystNuXSec_RvpNC1pi,    RvpNC1pi_def    );
  this->SetDefValue(kSystNuXSec_RvpNC2pi,    RvpNC2pi_def    );
  this->SetDefValue(kSystNuXSec_RvnCC1pi,    RvnCC1pi_def    );
  this->SetDefValue(kSystNuXSec_RvnCC2pi,    RvnCC2pi_def    );
  this->SetDefValue(kSystNuXSec_RvnNC1pi,    RvnNC1pi_def    );
  this->SetDefValue(kSystNuXSec_RvnNC2pi,    RvnNC2pi_def    );
  this->SetDefValue(kSystNuXSec_RvbarpCC1pi, RvbarpCC1pi_def );
  this->SetDefValue(kSystNuXSec_RvbarpCC2pi, RvbarpCC2pi_def );
  this->SetDefValue(kSystNuXSec_RvbarpNC1pi, RvbarpNC1pi_def );
  this->SetDefValue(kSystNuXSec_RvbarpNC2pi, RvbarpNC2pi_def );
  this->SetDefValue(kSystNuXSec_RvbarnCC1pi, RvbarnCC1pi_def );
  this->SetDefValue(kSystNuXSec_RvbarnCC2pi, RvbarnCC2pi_def );
  this->SetDefValue(kSystNuXSec_RvbarnNC1pi, RvbarnNC1pi_def );
  this->SetDefValue(kSystNuXSec_RvbarnNC2pi, RvbarnNC2pi_def );
}
//____________________________________________________________________________
void GReWeightNuXSecParams::Reset(GSyst_t syst) 
{
  if(this->IsIncluded(syst)) {
    double def_param = this->DefValue(syst);

    this->SetCurValue    (syst, def_param);
    this->SetCurTwkDial  (syst, 0.);
    this->SetTweakedFlag (syst, false);
  }
}
//____________________________________________________________________________
void GReWeightNuXSecParams::Reset(void) 
{
  map<GSyst_t, bool>::const_iterator iter = fIsTweaked.begin();

  for( ; iter != fIsTweaked.end(); ++iter) {
     GSyst_t syst = iter->first;
     this->Reset(syst);
  }
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetDefValue(GSyst_t syst, double value)
{
  fDefParams[syst] = value;
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetCurValue(GSyst_t syst, double value)
{
  fCurParams[syst] = value;
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetCurTwkDial(GSyst_t syst, double value)
{
  // check size of tweaking dial
  if(TMath::Abs(value) < controls::kASmallNum){ 
    this->SetTweakedFlag (syst, false); 
  }
  else {
    this->SetTweakedFlag (syst, true); 
  }

  // update tweaking dial
  fCurTwkDial[syst] = value;

  // update current value of corresponding physics parameter
  GSystUncertainty * syst_uncertainties = GSystUncertainty::Instance();
  double sigma = syst_uncertainties->OneSigmaErr(syst);
          
  double def_param = this->DefValue(syst);
  double cur_param = def_param * (1 + value * sigma);
          
  this->SetCurValue(syst, cur_param);
}
//____________________________________________________________________________
void GReWeightNuXSecParams::SetTweakedFlag(GSyst_t syst, bool value)
{
  fIsTweaked[syst] = value;
}
//____________________________________________________________________________
void GReWeightNuXSecParams::Init(void)
{
  fDefParams.clear ();
  fCurParams.clear ();
  fCurTwkDial.clear();
  fIsTweaked.clear ();
}      
//____________________________________________________________________________
