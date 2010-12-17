//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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

using namespace genie;
using namespace genie::rew;

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
  if (A<=1) return 1.;

  // Skip event if was not DIS scattering.
  bool is_dis = event.Summary()->ProcInfo().IsDeepInelastic();
  if(!is_dis) return 1.;

  //
  // Calculate event weight.
  //

  double event_weight = 1.;

  // nuclear radius for particle tracking purposes
  //double R = (fNR * fR0 * TMath::Power(A, 1./3.))/units::fm;

  // Get the event vtx needed to infer the formation zone.
  GHepParticle * hitnucl = event.HitNucleon();
  if(!hitnucl){ 
    LOG("ReW", pINFO)  << "Cannot calculate formation zone as target not a hit nucleon!"; 
    return 1.;
  }

  // hadronic system 4-/3-momentum
  TLorentzVector p4hadr = utils::rew::Hadronic4pLAB(event);
  TVector3 p3hadr = p4hadr.Vect();

  TLorentzVector * vtx = hitnucl->X4();

  // formation zone: fractional 1sigma err
  GSystUncertainty * uncertainty = GSystUncertainty::Instance();
  double fracerr = uncertainty->OneSigmaErr(kHadrNuclTwkDial_FormZone);

  // Loop over particles calculate weights for all primary hadrons inside the nucleus.
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     // Skip particles with code other than 'hadron in the nucleus'
     GHepStatus_t ist  = p->Status();
     if(ist != kIStHadronInTheNucleus)
     {
        continue;
     }
     // JIMTODO - Skip if is not a hadron - I am not sure why the above conditional would not catch this out - need to investigate.
     int            pdgc  = p->Pdg();       // hadron pdg code
     if (pdg::IsHadron(pdgc) == false) 
     {
        continue;
     } 
     
     // Default formation zone
     // Calculation below copied from the actual generation code.
     // Could also get default formation zone as the distance between the actual
     // vertex and the current position.
     TVector3 p3  = p->P4()->Vect();      // hadron's: p (px,py,pz)
     double   m   = p->Mass();            //           m
     double   m2  = m*m;                  //           m^2
     double   P   = p->P4()->P();         //           |p|
     double   Pt  = p3.Pt(p3hadr);        //           pT
     double   Pt2 = Pt*Pt;                //           pT^2

     LOG("ReW", pDEBUG)  << "FZone: m2 = "<< m2 << ", fK = "<< fK << ", Pt2 = "<< Pt2;

     double fz_def = P*fct0*m/(m2+fK*Pt2) / units::fm; 

     // Tweaked formation zone
     bool interacted = (p->RescatterCode() != -1);

     double fz_scale_factor  = (1 + fFZoneTwkDial * fracerr);
     LOG("ReW", pDEBUG)  << "Using fzone scale factor = "<< fz_scale_factor;
     double fz_twk  = fz_def * fz_scale_factor;
     fz_twk = TMath::Max(0.,fz_twk);
     LOG("ReW", pDEBUG)  << "FZone: "<< fz_def << ", FZone tweaked: "<< fz_twk;

     // hadron's: (px,py,pz,E)
     TLorentzVector * p4  = p->P4(); 
     TVector3 fz_vect = p4->Vect(); fz_vect.SetMag(fz_def);
     TLorentzVector fz_4vect(fz_vect, 0.0);
     TLorentzVector x4 = (*vtx) + fz_4vect;  
     LOG("ReW", pDEBUG)  << "Vtx positions (X,Y,Z,T): "<< 
         vtx->X() << ", "<< vtx->Y() << ", "<<vtx->Z() << ", "<< vtx->T();
     LOG("ReW", pDEBUG)  << "Hadron position (X,Y,Z,T): "<< 
           x4.X() << ", "<< x4.Y() << ", "<<x4.Z() << ", "<< x4.T();
 
     // Calculate particle weight
     double hadron_weight = genie::utils::rew::FZoneWeight(
        pdgc, *vtx, x4, *p4, A, fz_scale_factor, interacted);

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

  this->SetR0  (1.3 * units::fm);
  this->SetNR  (3.);
  this->SetCT0 (0.342 * units::fm);
  this->SetK   (0.);
}
//_______________________________________________________________________________________

