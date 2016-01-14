//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 25, 2010 - CA
   First included in v2.7.1.

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightNuXSecNC.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightNuXSecNC::GReWeightNuXSecNC() 
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightNuXSecNC::~GReWeightNuXSecNC()
{

}
//_______________________________________________________________________________________
bool GReWeightNuXSecNC::IsHandled(GSyst_t syst)
{
   switch(syst) {
     case ( kXSecTwkDial_NC ) : 
       return true;
       break;
     default:
       return false;
       break;
   }
   return false;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::SetSystematic(GSyst_t syst, double twk_dial)
{
  if(!this->IsHandled(syst)) return;

   switch(syst) {
     case ( kXSecTwkDial_NC ) : 
       fNCTwkDial = twk_dial;
       break;
     default:
       return;
       break;
   }
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::Reset(void)
{
  fNCTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::Reconfigure(void)
{
  fNCTwkDial = TMath::Max(0., fNCTwkDial);
}
//_______________________________________________________________________________________
double GReWeightNuXSecNC::CalcWeight(const genie::EventRecord & event) 
{
  bool tweaked = (TMath::Abs(fNCTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  Interaction * interaction = event.Summary();
  bool is_nc  = interaction->ProcInfo().IsWeakNC();
  if(!is_nc) return 1.;

  bool is_qel = event.Summary()->ProcInfo().IsQuasiElastic();
  if(is_qel && !fRewQE) return 1.;

  bool is_res = event.Summary()->ProcInfo().IsResonant();
  if(is_res && !fRewRES) return 1.;

  bool is_dis = event.Summary()->ProcInfo().IsDeepInelastic();
  if(is_dis && !fRewDIS) return 1.;

  int nupdg = interaction->InitState().ProbePdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  return fNCTwkDial;
}
//_______________________________________________________________________________________
double GReWeightNuXSecNC::CalcChisq(void)
{
  double chisq = TMath::Power(fNCTwkDial,  2.);
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightNuXSecNC::Init(void)
{
  fNCTwkDial = 0.;

  this->RewNue       (true);
  this->RewNuebar    (true);
  this->RewNumu      (true);
  this->RewNumubar   (true);

  this->RewQE  (true );
  this->RewRES (false); // assume GReWeightNuXSecNCRES is going to be your 1st choice
  this->RewDIS (false); // assume GReWeightNuXSecDIS   is going to be your 1st choice

}
//_______________________________________________________________________________________


