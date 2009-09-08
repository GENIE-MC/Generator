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
#include "ReWeight/GSystErrors.h"
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
  // Access the algorithm factory and the configuration pool

  fAlgFactory = genie::AlgFactory::Instance();
  fConfigPool = genie::AlgConfigPool::Instance();

  // Get the default cross section parameters (supported by the ReWeight package)

  fUserPhysicsConfig = fConfigPool->GlobalParameterList();
  fUserPhysicsConfig->UnLock();

  fUserPhysicsConfig->Get("QEL-Ma",                 fMaQEL_def);
  fUserPhysicsConfig->Get("QEL-Mv",                 fMvQEL_def);
  fUserPhysicsConfig->Get("RES-Ma",                 fMaRES_def);
  fUserPhysicsConfig->Get("RES-Mv",                 fMvRES_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vp-CC-m2",  fRvpCC1pi_def);    
  fUserPhysicsConfig->Get("DIS-HMultWgt-vp-CC-m3",  fRvpCC2pi_def);   
  fUserPhysicsConfig->Get("DIS-HMultWgt-vp-NC-m2",  fRvpNC1pi_def);   
  fUserPhysicsConfig->Get("DIS-HMultWgt-vp-NC-m3",  fRvpNC2pi_def);   
  fUserPhysicsConfig->Get("DIS-HMultWgt-vn-CC-m2",  fRvnCC1pi_def);   
  fUserPhysicsConfig->Get("DIS-HMultWgt-vn-CC-m3",  fRvnCC2pi_def);   
  fUserPhysicsConfig->Get("DIS-HMultWgt-vn-NC-m2",  fRvnNC1pi_def);   
  fUserPhysicsConfig->Get("DIS-HMultWgt-vn-NC-m3",  fRvnNC2pi_def);   
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbp-CC-m2", fRvbarpCC1pi_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbp-CC-m3", fRvbarpCC2pi_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbp-NC-m2", fRvbarpNC1pi_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbp-NC-m3", fRvbarpNC2pi_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbn-CC-m2", fRvbarnCC1pi_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbn-CC-m3", fRvbarnCC2pi_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbn-NC-m2", fRvbarnNC1pi_def);
  fUserPhysicsConfig->Get("DIS-HMultWgt-vbn-NC-m3", fRvbarnNC2pi_def);

  fMaQEL       = fMaQEL_def;
  fMvQEL       = fMvQEL_def;
  fMaRES       = fMaRES_def;
  fMvRES       = fMvRES_def;
  fRvpCC1pi    = fRvpCC1pi_def;   
  fRvpCC2pi    = fRvpCC2pi_def;   
  fRvpNC1pi    = fRvpNC1pi_def;   
  fRvpNC2pi    = fRvpNC2pi_def;   
  fRvnCC1pi    = fRvnCC1pi_def;   
  fRvnCC2pi    = fRvnCC2pi_def;   
  fRvnNC1pi    = fRvnNC1pi_def;   
  fRvnNC2pi    = fRvnNC2pi_def;   
  fRvbarpCC1pi = fRvbarpCC1pi_def;
  fRvbarpCC2pi = fRvbarpCC2pi_def;
  fRvbarpNC1pi = fRvbarpNC1pi_def;
  fRvbarpNC2pi = fRvbarpNC2pi_def;
  fRvbarnCC1pi = fRvbarnCC1pi_def;
  fRvbarnCC2pi = fRvbarnCC2pi_def;
  fRvbarnNC1pi = fRvbarnNC1pi_def;
  fRvbarnNC2pi = fRvbarnNC2pi_def;
  
  fTweaked              = false;
  fMaQEL_included       = false;
  fMvQEL_included       = false;
  fMaRES_included       = false;
  fMvRES_included       = false;
  fRvpCC1pi_included    = false;    
  fRvpCC2pi_included    = false;   
  fRvpNC1pi_included    = false;   
  fRvpNC2pi_included    = false;   
  fRvnCC1pi_included    = false;   
  fRvnCC2pi_included    = false;   
  fRvnNC1pi_included    = false;   
  fRvnNC2pi_included    = false;   
  fRvbarpCC1pi_included = false;
  fRvbarpCC2pi_included = false;
  fRvbarpNC1pi_included = false;
  fRvbarpNC2pi_included = false;
  fRvbarnCC1pi_included = false;
  fRvbarnCC2pi_included = false;
  fRvbarnNC1pi_included = false;
  fRvbarnNC2pi_included = false;
}
//_______________________________________________________________________________________
void GReWeightNuXSec::SetSystematic(GSyst_t syst, double val)
{
  double eta  = 1E-15;
  double val2 = val*val;

  switch(syst) {

  case ( kSystNuXSec_MaQEL ) : 
    fMaQEL = fMaQEL_def*(1.0+val*kNuXSec_MaQELFracError);
    fNuXSec_MaQELTwkDial = val; 
    fMaQEL_included = false; 
    if(val2 > eta) {
      fMaQEL_included = true; 
      fTweaked = true;
    }
    break;

  case ( kSystNuXSec_MvQEL ) : 
    fMvQEL = fMvQEL_def*(1.0+val*kNuXSec_MvQELFracError);
    fNuXSec_MvQELTwkDial = val;
    fMvQEL_included = false; 
    if(val2 > eta) {
      fMvQEL_included = true; 
      fTweaked = true;
    }
    break;

  case ( kSystNuXSec_MaRES ) :
    fMaRES = fMaRES_def*(1.0+val*kNuXSec_MaRESFracError);  
    fNuXSec_MaRESTwkDial = val;
    fMaRES_included = false; 
    if(val2 > eta) {
      fMaRES_included = true; 
      fTweaked = true;
    }
    break;

  case ( kSystNuXSec_MvRES ) : 
    fMvRES = fMvRES_def*(1.0+val*kNuXSec_MvRESFracError); 
    fNuXSec_MvRESTwkDial = val;
    fMvRES_included = false; 
    if(val2 > eta) {
      fMvRES_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvpCC1pi):  
    fRvpCC1pi = fRvpCC1pi_def*(1.0+val*kNuXSec_RvpCC1piFracError);  
    fTweaked = true;
    fNuXSec_RvpCC1piTwkDial = val;
    fRvpCC1pi_included = false; 
    if(val2 > eta) {
      fRvpCC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvpCC2pi):  
    fRvpCC2pi = fRvpCC2pi_def*(1.0+val*kNuXSec_RvpCC2piFracError);
    fNuXSec_RvpCC2piTwkDial = val;
    fRvpCC2pi_included = false; 
    if(val2 > eta) {
      fRvpCC2pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvpNC1pi): 
    fRvpNC1pi = fRvpNC1pi_def*(1.0+val*kNuXSec_RvpNC1piFracError); 
    fNuXSec_RvpNC1piTwkDial = val; 
    fRvpNC1pi_included = false; 
    if(val2 > eta) {
      fRvpNC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvpNC2pi): 
    fRvpNC2pi = fRvpNC2pi_def*(1.0+val*kNuXSec_RvpNC2piFracError); 
    fNuXSec_RvpNC2piTwkDial = val;
    fRvpNC2pi_included = false; 
    if(val2 > eta) {
      fRvpNC2pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvnCC1pi): 
    fRvnCC1pi = fRvnCC1pi_def*(1.0+val*kNuXSec_RvnCC1piFracError);  
    fNuXSec_RvnCC1piTwkDial = val; 
    fRvnCC1pi_included = false; 
    if(val2 > eta) {
      fRvnCC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvnCC2pi): 
    fRvnCC2pi = fRvnCC2pi_def*(1.0+val*kNuXSec_RvnCC2piFracError); 
    fNuXSec_RvnCC2piTwkDial = val;
    fRvnCC2pi_included = false; 
    if(val2 > eta) {
      fRvnCC2pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvnNC1pi):
    fRvnNC1pi = fRvnNC1pi_def*(1.0+val*kNuXSec_RvnNC1piFracError);
    fNuXSec_RvnNC1piTwkDial = val;
    fRvnNC1pi_included = false; 
    if(val2 > eta) {
      fRvnNC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvnNC2pi):
    fRvnNC2pi = fRvnNC2pi_def*(1.0+val*kNuXSec_RvnNC2piFracError);  
    fNuXSec_RvnNC2piTwkDial = val; 
      fRvnNC2pi_included = false; 
    if(val2 > eta) {
      fRvnNC2pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarpCC1pi): 
    fRvbarpCC1pi = fRvbarpCC1pi_def*(1.0+val*kNuXSec_RvbarpCC1piFracError); 
    fNuXSec_RvbarpCC1piTwkDial = val;
    fRvbarpCC1pi_included = false;
    if(val2 > eta) {
      fRvbarpCC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarpCC2pi): 
    fRvbarpCC2pi = fRvbarpCC2pi_def*(1.0+val*kNuXSec_RvbarpCC2piFracError); 
    fNuXSec_RvbarpCC2piTwkDial = val;
    fRvbarpCC2pi_included = false;
    if(val2 > eta) {
      fRvbarpCC2pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarpNC1pi): 
    fRvbarpNC1pi = fRvbarpNC1pi_def*(1.0+val*kNuXSec_RvbarpNC1piFracError);  
    fNuXSec_RvbarpNC1piTwkDial = val; 
    fRvbarpNC1pi_included = false;
    if(val2 > eta) {
      fRvbarpNC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarpNC2pi):
    fRvbarpNC2pi = fRvbarpNC2pi_def*(1.0+val*kNuXSec_RvbarpNC2piFracError); 
    fNuXSec_RvbarpNC2piTwkDial = val;
    fRvbarpNC2pi_included = false;
    if(val2 > eta) {
      fRvbarpNC2pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarnCC1pi): 
    fRvbarnCC1pi = fRvbarnCC1pi_def*(1.0+val*kNuXSec_RvbarnCC1piFracError); 
    fNuXSec_RvbarnCC1piTwkDial = val;
    fRvbarnCC1pi_included = false;
    if(val2 > eta) {
      fRvbarnCC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarnCC2pi): 
    fRvbarnCC2pi = fRvbarnCC2pi_def*(1.0+val*kNuXSec_RvbarnCC2piFracError); 
    fNuXSec_RvbarnCC2piTwkDial = val;
    fRvbarnCC2pi_included = false;
    if(val2 > eta) {
      fRvbarnCC2pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarnNC1pi): 
    fRvbarnNC1pi = fRvbarnNC1pi_def*(1.0+val*kNuXSec_RvbarnNC1piFracError);  
    fNuXSec_RvbarnNC1piTwkDial = val;
    fRvbarnNC1pi_included = false;
    if(val2 > eta) {
      fRvbarnNC1pi_included = true; 
      fTweaked = true;
    }
    break;

  case (kSystNuXSec_RvbarnNC2pi): 
    fRvbarnNC2pi = fRvbarnNC2pi_def*(1.0+val*kNuXSec_RvbarnNC2piFracError); 
    fNuXSec_RvbarnNC2piTwkDial = val;
    fRvbarnNC2pi_included = false;
    if(val2 > eta) {
      fRvbarnNC2pi_included = true;
      fTweaked = true;
    }
    break;

  default:
    break;
  }
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reset(void)
{
// restore 'factory' defaults

  fMaQEL = fMaQEL_def;
  fMvQEL = fMvQEL_def;
  fMaRES = fMaRES_def;
  fMvRES = fMvRES_def;

  // need to add resonance background defaults



  this->Reconfigure();

  fTweaked = false;
  fMaQEL_included       = false;
  fMvQEL_included       = false;
  fMaRES_included       = false;
  fMvRES_included       = false;
  fRvpCC1pi_included    = false;    
  fRvpCC2pi_included    = false;   
  fRvpNC1pi_included    = false;   
  fRvpNC2pi_included    = false;   
  fRvnCC1pi_included    = false;   
  fRvnCC2pi_included    = false;   
  fRvnNC1pi_included    = false;   
  fRvnNC2pi_included    = false;   
  fRvbarpCC1pi_included = false;
  fRvbarpCC2pi_included = false;
  fRvbarpNC1pi_included = false;
  fRvbarpNC2pi_included = false;
  fRvbarnCC1pi_included = false;
  fRvbarnCC2pi_included = false;
  fRvbarnNC1pi_included = false;
  fRvbarnNC2pi_included = false;
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reconfigure(void)
{
  if(!fTweaked) return;

  fUserPhysicsConfig->Set("QEL-Ma", fMaQEL);
  fUserPhysicsConfig->Set("QEL-Mv", fMvQEL);
  fUserPhysicsConfig->Set("RES-Ma", fMaRES);
  fUserPhysicsConfig->Set("RES-Mv", fMvRES);

  fAlgFactory->ForceReconfiguration();
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcWeight(const genie::EventRecord & event) 
{
  if (!fTweaked) return 1.;

  double wght = fXSecRwHelper.NewWeight(event);
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcChisq()
{
  double chisq = 0;

  chisq += TMath::Power(fNuXSec_MaQELTwkDial,       2.0);    
  chisq += TMath::Power(fNuXSec_MvQELTwkDial,       2.0);    
  chisq += TMath::Power(fNuXSec_MaRESTwkDial,       2.0);
  chisq += TMath::Power(fNuXSec_MvRESTwkDial,       2.0);
  chisq += TMath::Power(fNuXSec_RvpCC1piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvpCC2piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvpNC1piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvpNC2piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvnCC1piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvnCC2piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvnNC1piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvnNC2piTwkDial,    2.0);
  chisq += TMath::Power(fNuXSec_RvbarpCC1piTwkDial, 2.0);
  chisq += TMath::Power(fNuXSec_RvbarpCC2piTwkDial, 2.0);
  chisq += TMath::Power(fNuXSec_RvbarpNC1piTwkDial, 2.0);
  chisq += TMath::Power(fNuXSec_RvbarpNC2piTwkDial, 2.0);
  chisq += TMath::Power(fNuXSec_RvbarnCC1piTwkDial, 2.0);
  chisq += TMath::Power(fNuXSec_RvbarnCC2piTwkDial, 2.0);
  chisq += TMath::Power(fNuXSec_RvbarnNC1piTwkDial, 2.0);
  chisq += TMath::Power(fNuXSec_RvbarnNC2piTwkDial, 2.0);

  return chisq;
}
//_______________________________________________________________________________________
