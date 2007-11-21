//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 21, 2007 - CA
   Was renamed to COHPiXSec (from COHXSec)

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "CrossSections/COHPiXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
COHPiXSec::COHPiXSec() :
XSecIntegratorI("genie::COHPiXSec")
{

}
//____________________________________________________________________________
COHPiXSec::COHPiXSec(string config) :
XSecIntegratorI("genie::COHPiXSec", config)
{

}
//____________________________________________________________________________
COHPiXSec::~COHPiXSec()
{

}
//____________________________________________________________________________
double COHPiXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("COHPiXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  LOG("COHPiXSec", pINFO)
            << "x integration range = [" << xl.min << ", " << xl.max << "]";
  LOG("COHPiXSec", pINFO)
            << "y integration range = [" << yl.min << ", " << yl.max << "]";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  GXSecFunc * func = new Integrand_D2XSec_DxDy_E(model, interaction);
  func->SetParam(0,"x",xl);
  func->SetParam(1,"y",yl);
  double xsec = fIntegrator->Integrate(*func);

  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
  LOG("COHPiXSec", pINFO)  << "XSec[COHPi] (E = " << Ev << " GeV) = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
void COHPiXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHPiXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHPiXSec::LoadConfig(void)
{
  //-- get specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
  assert(fIntegrator);
}
//____________________________________________________________________________
