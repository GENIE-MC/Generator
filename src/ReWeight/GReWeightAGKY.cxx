//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Skeleton first included in v2.5.1.

*/
//____________________________________________________________________________

#include "EVGCore/EventRecord.h"
#include "Messenger/Messenger.h"
#include "ReWeight/GReWeightAGKY.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightAGKY::GReWeightAGKY() :
GReWeightI()
{

}
//_______________________________________________________________________________________
GReWeightAGKY::~GReWeightAGKY()
{

}
//_______________________________________________________________________________________
void GReWeightAGKY::SetSystematic(GSyst_t /*syst*/, double /*val*/)
{

}
//_______________________________________________________________________________________
void GReWeightAGKY::Reset(void)
{

}
//_______________________________________________________________________________________
void GReWeightAGKY::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightAGKY::CalcWeight(const EventRecord & /*event*/) 
{  
  return 1.;
}
//_______________________________________________________________________________________
double GReWeightAGKY::CalcChisq(void)
{
  return 0.;
}
//_______________________________________________________________________________________
