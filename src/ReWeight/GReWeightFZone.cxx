//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
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

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightFZone.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightFZone::GReWeightFZone() :
GReWeightI()
{

}
//_______________________________________________________________________________________
GReWeightFZone::~GReWeightFZone()
{

}
//_______________________________________________________________________________________
void GReWeightFZone::SetSystematic(GSyst_t /*syst*/, double /*val*/)
{

}
//_______________________________________________________________________________________
void GReWeightFZone::Reset(void)
{

}
//_______________________________________________________________________________________
void GReWeightFZone::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightFZone::CalcWeight(const EventRecord & /*event*/) 
{  
  return 1.0;
}
//_______________________________________________________________________________________
double GReWeightFZone::CalcChisq(void)
{
  return 0.;
}
//_______________________________________________________________________________________
