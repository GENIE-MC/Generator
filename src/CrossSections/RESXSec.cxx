//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "CrossSections/RESXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESXSec::RESXSec() :
XSecIntegratorI("genie::RESXSec")
{

}
//____________________________________________________________________________
RESXSec::RESXSec(string config) :
XSecIntegratorI("genie::RESXSec", config)
{

}
//____________________________________________________________________________
RESXSec::~RESXSec()
{

}
//____________________________________________________________________________
double RESXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t Q2l = kps.Limits(kKVQ2);
  Range1D_t Wl  = kps.Limits(kKVW);

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  GXSecFunc * func = new Integrand_D2XSec_DWDQ2_E(model, interaction);
  func->SetParam(0,"W",  Wl);
  func->SetParam(1,"Q2", Q2l);
  double xsec = fIntegrator->Integrate(*func);

  //SLOG("RESXSec", pINFO)  << "XSec[RES] (Ev = " << Ev << " GeV) = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
void RESXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESXSec::LoadConfig(void)
{
  //-- get the specified integration algorithm
  fIntegrator = 
       dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
  assert(fIntegrator);
}
//____________________________________________________________________________
