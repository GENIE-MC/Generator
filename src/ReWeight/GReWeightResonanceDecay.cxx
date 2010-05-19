//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 10, 2010 - CA
   First included in v2.7.1.

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightResonanceDecay::GReWeightResonanceDecay() :
GReWeightI()
{
  this->Reset();
}
//_______________________________________________________________________________________
GReWeightResonanceDecay::~GReWeightResonanceDecay()
{

}
//_______________________________________________________________________________________
bool GReWeightResonanceDecay::IsHandled(GSyst_t syst)
{
  switch(syst) {
     case ( kRDcyTwkDial_BR_Delta2Ngamma  ) :
     case ( kRDcyTwkDial_Theta_Delta2Npi  ) :
        return true;
        break;
     default:
        return false;
        break;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightResonanceDecay::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
     case ( kRDcyTwkDial_BR_Delta2Ngamma  ) :
        fBRDelta2NgammaTwkDial = twk_dial;
        break;
     case ( kRDcyTwkDial_Theta_Delta2Npi  ) :
        fThetaDelta2NpiTwkDial = twk_dial;
        break;
     default:
        return;
        break;
  }
}
//_______________________________________________________________________________________
void GReWeightResonanceDecay::Reset(void)
{
  fBRDelta2NgammaTwkDial = 0.0;
  fThetaDelta2NpiTwkDial = 0.0;
}
//_______________________________________________________________________________________
void GReWeightResonanceDecay::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::CalcWeight(const EventRecord & event) 
{  
  double wght = 
    this->RewBRDelta2Ngamma(event) * 
    this->RewThetaDelta2Npi(event);

  return wght;
}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::CalcChisq(void)
{
  return 0.;
}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::RewBRDelta2Ngamma(const EventRecord & event)
{
  bool tweaked = (TMath::Abs(fBRDelta2NgammaTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  bool is_res = event.Summary()->ProcInfo().IsResonant();
  if(!is_res) return 1.;

  bool is_Ngamma = false;
  GHepParticle * p = 0;
  TIter iter(&event);
  while((p=(GHepParticle*)iter.Next())) {
    if(utils::res::IsBaryonResonance(p->Pdg())) {
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if(nd==2) {
        int fpdg = event.Particle(fd)->Pdg();
        int lpdg = event.Particle(ld)->Pdg();
        if((fpdg==kPdgProton||fpdg==kPdgNeutron) && lpdg==kPdgGamma) is_Ngamma = true;
        if((lpdg==kPdgProton||lpdg==kPdgNeutron) && fpdg==kPdgGamma) is_Ngamma = true;
      }
      if(is_Ngamma) break;
    }
  }

  if(!is_Ngamma) return 1.;

  if(is_Ngamma) {
    LOG("ReW", pDEBUG) << "A Resonance -> N + gamma event:";
    LOG("ReW", pDEBUG) << event;
  }

//GSystUncertainty * uncertainty = GSystUncertainty::Instance();

  double err  = 1;
  double dial = fBRDelta2NgammaTwkDial;

  double wght = (1. + dial * err);

  return wght;
}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::RewThetaDelta2Npi(const EventRecord & event)
{
  bool tweaked = (TMath::Abs(fThetaDelta2NpiTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  bool is_res = event.Summary()->ProcInfo().IsResonant();
  if(!is_res) return 1.;

  bool is_Deltapp_1pi = false;
  GHepParticle * p = 0;
  TIter iter(&event);
  while((p=(GHepParticle*)iter.Next())) {
    bool is_Deltapp = (p->Pdg()==kPdgP33m1232_DeltaPP);
    if(is_Deltapp) {
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if(nd==2) {
        int fpdg = event.Particle(fd)->Pdg();
        int lpdg = event.Particle(ld)->Pdg();
        if(fpdg==kPdgProton && lpdg==kPdgPiP) is_Deltapp_1pi = true;
        if(lpdg==kPdgProton && fpdg==kPdgPiP) is_Deltapp_1pi = true;
      }
      if(is_Deltapp_1pi) break;
    }
  }
  if(!is_Deltapp_1pi) return 1.;

  if(is_Deltapp_1pi) {
    LOG("ReW", pDEBUG) << "A Delta++ -> p pi+ event:";
    LOG("ReW", pDEBUG) << event;
  }

  return 0.;
}
//_______________________________________________________________________________________
