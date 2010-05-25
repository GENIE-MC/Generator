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
     
     // Calculate particle weight
     double hadron_weight = 1.;
     //     
     // ...
     // ...
     // ...
     //     


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
}
//_______________________________________________________________________________________

