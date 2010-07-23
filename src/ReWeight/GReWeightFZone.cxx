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

*/
//____________________________________________________________________________

#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
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

  //
  // Calculate event weight.
  //

  double event_weight = 1.;

  // nuclear radius for particle tracking purposes
  //double R = (fNR * fR0 * TMath::Power(A, 1./3.))/units::fm;

  // hadronic system 4-/3-momentum
  TLorentzVector p4hadr = utils::rew::Hadronic4pLAB(event);
  TVector3 p3hadr = p4hadr.Vect();

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

     double fz_def = P*fct0*m/(m2+fK*Pt2) / units::fm; 

     // Tweaked formation zone
     double fz_twk  = fz_def * (1 + fFZoneTwkDial * fracerr);
     fz_twk = TMath::Max(0.,fz_twk);

     // Calculate particle weight
     double hadron_weight = 1.;

     //     
     // ...
     // ...
     // ...
     //     


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

