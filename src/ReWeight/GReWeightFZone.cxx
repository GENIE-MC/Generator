//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 20, 2009 - CA
   Skeleton first included in v2.5.1.
 @ Dec 17, 2010 - JD
   First implementation of FZone reweighting.
 @ Feb 08, 2013 - CA
   Mean free path is function of Z too. Make sure it is passed on.
   Use new INTRANUKE FSI flags.	Use formation zone code from PhysUtils.
*/
//____________________________________________________________________________

#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightUtils.h"
#include "ReWeight/GSystUncertainty.h"
#include "Utils/PrintUtils.h"
#include "Utils/PhysUtils.h"

using namespace genie;
using namespace genie::rew;
using namespace genie::utils;

//_______________________________________________________________________________________
GReWeightFZone::GReWeightFZone() :
GReWeightI()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightFZone::~GReWeightFZone()
{

}
//_______________________________________________________________________________________
bool GReWeightFZone::IsHandled(GSyst_t syst)
{
  switch(syst) {
    case (kHadrNuclTwkDial_FormZone) :
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightFZone::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
    case (kHadrNuclTwkDial_FormZone) :
      fFZoneTwkDial = twk_dial;
      break;
    default:
      return;
  }
}
//_______________________________________________________________________________________
void GReWeightFZone::Reset(void)
{
  fFZoneTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightFZone::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightFZone::CalcWeight(const EventRecord & event) 
{  
  // Physics parameter tweaked?
  bool tweaked = (TMath::Abs(fFZoneTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  // Skip events not involving nuclear targets.
  GHepParticle * tgt = event.TargetNucleus();
  if (!tgt) return 1.;
  double A = tgt->A();
  double Z = tgt->Z();
  if (A<=1) return 1.;
  if (Z<=1) return 1.0;


  // Skip event if was not DIS scattering.
  bool is_dis = event.Summary()->ProcInfo().IsDeepInelastic();
  if(!is_dis) {
    LOG("ReW", pDEBUG) << "Not a DIS event";
    return 1.;
  }

  //
  // Calculate event weight.
  //

  double event_weight = 1.;

  // hadronic system 4-/3-momentum
  TLorentzVector p4hadr = utils::rew::Hadronic4pLAB(event);
  TVector3 p3hadr = p4hadr.Vect();

  // vertex
  assert(event.HitNucleon());
  const TLorentzVector & vtx = *(event.HitNucleon()->X4());

  // formation zone: fractional 1sigma err
  GSystUncertainty * uncertainty = GSystUncertainty::Instance();
  double fracerr = uncertainty->OneSigmaErr(kHadrNuclTwkDial_FormZone);

  // Loop over particles calculate weights for all primary hadrons inside the nucleus.
  int ip=-1;
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
      ip++;

     // Skip particles with code other than 'hadron in the nucleus'
     GHepStatus_t ist  = p->Status();
     if(ist != kIStHadronInTheNucleus)
     {
        continue;
     }
     // JIMTODO - Skip if is not a hadron 
     // I am not sure why the above conditional would not catch this out - need to investigate.
     int pdgc = p->Pdg();   // hadron pdg code
     if (pdg::IsHadron(pdgc) == false) 
     {
        continue;
     } 

     // Determine whether it interacted or not
     int fsi_code = p->RescatterCode();
     if(fsi_code == -1 || fsi_code == (int)kIHAFtUndefined) {
       LOG("ReW", pFATAL) << "INTRANUKE didn't set a valid rescattering code for event in position: " << ip;
       LOG("ReW", pFATAL) << "Here is the problematic event:";
       LOG("ReW", pFATAL) << event;
//       exit(1);
      fsi_code = kIHAFtNoInteraction;
     }
     bool escaped    = (fsi_code == (int)kIHAFtNoInteraction);
     bool interacted = !escaped;

     LOG("ReW", pDEBUG)
        << "Attempting to reweight hadron at position = " << ip
        << " with PDG code = " << pdgc
        << ". The hadron "
        << ((interacted) ? "re-interacted" : "did not re-interact");
     
     // Default formation zone
     double m = p->Mass();
     TLorentzVector * p4  = p->P4();
 
     double ct0=0.;
     pdg::IsNucleon(pdgc) ? ct0=fct0nucleon : ct0=fct0pion;
     double fz_def = phys::FormationZone(m,*p4,p3hadr,ct0,fK);

     double fz_scale_factor  = (1 + fFZoneTwkDial * fracerr);
     double fz_twk  = fz_def * fz_scale_factor;
     fz_twk = TMath::Max(0.,fz_twk);

     LOG("ReW", pDEBUG)  
        << "Formation zone = " << fz_def << " fm (nominal), "
        << fz_twk << " fm (tweaked) - scale factor = " << fz_scale_factor;

     // Calculate hadron's position at the end of the formation zone step
     TVector3 step3v = p4->Vect(); 
     step3v.SetMag(fz_def);
     TLorentzVector step4v(step3v, 0.);

     TLorentzVector x4 = vtx + step4v;  
     LOG("ReW", pDEBUG)  << "Vtx position: "<< print::X4AsString(&vtx);
     LOG("ReW", pDEBUG)  << "Hadron position: "<< print::X4AsString(&x4);
 
     // Calculate particle weight
     double hadron_weight = genie::utils::rew::FZoneWeight(
        pdgc, vtx, x4, *p4, A, Z, fz_scale_factor, interacted);

     // Update event weight
     event_weight *= hadron_weight;
  }

  return event_weight;
}
//_______________________________________________________________________________________
double GReWeightFZone::CalcChisq(void)
{
  double chisq = TMath::Power(fFZoneTwkDial, 2.);
  return chisq;
}
//_______________________________________________________________________________________
void GReWeightFZone::Init(void)
{
  fFZoneTwkDial = 0.;

  this->SetR0         (1.3);//fm
  this->SetNR         (3.);
  this->SetCT0Pion    (0.342);//fm
  this->SetCT0Nucleon (2.300);//fm
  this->SetK          (0.);
}
//_______________________________________________________________________________________

