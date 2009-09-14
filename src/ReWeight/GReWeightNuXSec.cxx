//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code. 
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "BaryonResonance/BaryonResonance.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightNuXSec.h"
#include "ReWeight/GSystUncertainty.h"
#include "ReWeight/GSystSet.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSec::GReWeightNuXSec() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSec::~GReWeightNuXSec()
{

}
//_______________________________________________________________________________________
void GReWeightNuXSec::Init(void)
{
  // Get the default cross section parameters (supported by the ReWeight package)
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

  // Keep local copies of the defaults
  fXSecRwParams.SetDefValue(kSystNuXSec_MaQEL,       MaQEL_def       );
  fXSecRwParams.SetDefValue(kSystNuXSec_MvQEL,       MvQEL_def       );
  fXSecRwParams.SetDefValue(kSystNuXSec_MaRES,       MaRES_def       );
  fXSecRwParams.SetDefValue(kSystNuXSec_MvRES,       MvRES_def       );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvpCC1pi,    RvpCC1pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvpCC2pi,    RvpCC2pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvpNC1pi,    RvpNC1pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvpNC2pi,    RvpNC2pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvnCC1pi,    RvnCC1pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvnCC2pi,    RvnCC2pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvnNC1pi,    RvnNC1pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvnNC2pi,    RvnNC2pi_def    );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarpCC1pi, RvbarpCC1pi_def );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarpCC2pi, RvbarpCC2pi_def );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarpNC1pi, RvbarpNC1pi_def );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarpNC2pi, RvbarpNC2pi_def );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarnCC1pi, RvbarnCC1pi_def );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarnCC2pi, RvbarnCC2pi_def );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarnNC1pi, RvbarnNC1pi_def );
  fXSecRwParams.SetDefValue(kSystNuXSec_RvbarnNC2pi, RvbarnNC2pi_def );
}
//_______________________________________________________________________________________
void GReWeightNuXSec::SetSystematic(GSyst_t syst, double twk_dial)
{
  GSystUncertainty * syst_uncertainties = GSystUncertainty::Instance();
  double sigma = syst_uncertainties->OneSigmaErr(syst);

  double def_param = fXSecRwParams.DefValue(syst);
  double cur_param = def_param * (1 + twk_dial * sigma);

  fXSecRwParams.SetCurValue   (syst, cur_param);
  fXSecRwParams.SetCurTwkDial (syst, twk_dial );

  double eta  = 1E-07;
  if(TMath::Abs(twk_dial) > eta) {
     fXSecRwParams.SetTweakedFlag (syst, true);
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reset(void)
{
// restore 'factory' defaults
//
  fXSecRwParams.Reset ( kSystNuXSec_MaQEL       );
  fXSecRwParams.Reset ( kSystNuXSec_MvQEL       );
  fXSecRwParams.Reset ( kSystNuXSec_MaRES       );
  fXSecRwParams.Reset ( kSystNuXSec_MvRES       );
  fXSecRwParams.Reset ( kSystNuXSec_RvpCC1pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvpCC2pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvpNC1pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvpNC2pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvnCC1pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvnCC2pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvnNC1pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvnNC2pi    );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarpCC1pi );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarpCC2pi );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarpNC1pi );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarpNC2pi );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarnCC1pi );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarnCC2pi );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarnNC1pi );
  fXSecRwParams.Reset ( kSystNuXSec_RvbarnNC2pi );

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reconfigure(void)
{
  if (! fXSecRwParams.IsTweaked() ) return;

  // Use current reweighting params inputs to tweak the GENIE configuration
  AlgConfigPool * conf_pool = AlgConfigPool::Instance();
  Registry * user_config = conf_pool->GlobalParameterList();

  if(fXSecRwParams.IsIncluded(kSystNuXSec_MaQEL)) {
  	user_config->Set("QEL-Ma", fXSecRwParams.CurValue(kSystNuXSec_MaQEL));
  }
  if(fXSecRwParams.IsIncluded(kSystNuXSec_MvQEL)) {
  	user_config->Set("QEL-Mv", fXSecRwParams.CurValue(kSystNuXSec_MvQEL));
  }
  if(fXSecRwParams.IsIncluded(kSystNuXSec_MaRES)) {
  	user_config->Set("RES-Ma", fXSecRwParams.CurValue(kSystNuXSec_MaRES));
  }
  if(fXSecRwParams.IsIncluded(kSystNuXSec_MvRES)) {
  	user_config->Set("RES-Mv", fXSecRwParams.CurValue(kSystNuXSec_MvRES));
  }


  //
  // note: *** need to add resonance background params ***
  //


  // Force reconfiguration of GENIE algorithms
  AlgFactory * algf = AlgFactory::Instance();
  algf->ForceReconfiguration();
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcWeight(const genie::EventRecord & event) 
{
  if (! fXSecRwParams.IsTweaked() ) return 1.;

  double wght = fXSecRwHelper.NewWeight(event);
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcChisq()
{
  double chisq = 0;

  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_MaQEL),         2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_MvQEL),         2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_MaRES),         2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_MvRES),         2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvpCC1pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvpCC2pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvpNC1pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvpNC2pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvnCC1pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvnCC2pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvnNC1pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvnNC2pi),      2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarpCC1pi),   2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarpCC2pi),   2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarpNC1pi),   2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarpNC2pi),   2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarnCC1pi),   2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarnCC2pi),   2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarnNC1pi),   2.0);    
  chisq += TMath::Power( fXSecRwParams.CurTwkDial(kSystNuXSec_RvbarnNC2pi),   2.0);    

  return chisq;
}
//_______________________________________________________________________________________
