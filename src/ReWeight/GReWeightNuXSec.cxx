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
  // Get the default cross section parameters 
  fXSecRwParams.LoadDefaults();
}
//_______________________________________________________________________________________
void GReWeightNuXSec::SetSystematic(GSyst_t syst, double twk_dial)
{
  fXSecRwParams.SetCurTwkDial (syst, twk_dial );
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reset(void)
{
  fXSecRwParams.Reset();
  fXSecRwParams.Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightNuXSec::Reconfigure(void)
{
  fXSecRwParams.Reconfigure();
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcWeight(const genie::EventRecord & event) 
{
  if (! fXSecRwParams.IsTweaked() ) return 1.;

  double wght = fXSecRwHelper.NewWeight(event);
  return wght;
}
//_______________________________________________________________________________________
double GReWeightNuXSec::CalcChisq()
{
  return fXSecRwParams.ChisqPenalty();
}
//_______________________________________________________________________________________
