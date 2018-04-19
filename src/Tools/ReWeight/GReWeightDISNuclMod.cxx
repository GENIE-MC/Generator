//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Apr 27, 2010 - CA
   First included in v2.7.1.

*/
//____________________________________________________________________________

#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Tools/ReWeight/GReWeightDISNuclMod.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightDISNuclMod::GReWeightDISNuclMod() :
GReWeightModel("DISNuclMod")
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightDISNuclMod::~GReWeightDISNuclMod()
{

}
//_______________________________________________________________________________________
bool GReWeightDISNuclMod::IsHandled(GSyst_t syst)
{
  switch(syst) {
    case(kXSecTwkDial_DISNuclMod) :
      return true;
      break;
    default:
      return false;
      break;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
    case(kXSecTwkDial_DISNuclMod) :
      fNuclModTwkDial = twk_dial;
      break;
    default:
      break;
  }
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::Reset(void)
{
  fNuclModTwkDial = 0.;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::Reconfigure(void)
{
  
}
//_______________________________________________________________________________________
double GReWeightDISNuclMod::CalcWeight(const EventRecord & /*event*/) 
{ 
  LOG("ReW",pFATAL) << "Not implemented.";
  exit(-1);
  return 1.;
}
//_______________________________________________________________________________________
void GReWeightDISNuclMod::Init(void)
{
  fNuclModTwkDial = 0.;

}
//_______________________________________________________________________________________
