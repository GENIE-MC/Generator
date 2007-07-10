//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "CrossSections/DISXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISXSec::DISXSec() :
XSecIntegratorI("genie::DISXSec")
{

}
//____________________________________________________________________________
DISXSec::DISXSec(string config) :
XSecIntegratorI("genie::DISXSec", config)
{

}
//____________________________________________________________________________
DISXSec::~DISXSec()
{

}
//____________________________________________________________________________
double DISXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("DISXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  Interaction * interaction = new Interaction(*in);

  interaction->SetBit(kISkipProcessChk);
  //interaction->SetBit(kISkipKinematicChk);

  interaction->SetBit(kINoNuclearCorrection);

  Range1D_t Wl  = kps.WLim();
  Range1D_t Q2l = kps.Q2Lim();
  LOG("DISXSec", pINFO)  
            << "W integration range = [" << Wl.min << ", " << Wl.max << "]";
  LOG("DISXSec", pINFO)  
         << "Q2 integration range = [" << Q2l.min << ", " << Q2l.max << "]";

  GXSecFunc * func = new Integrand_D2XSec_DWDQ2_E(model, interaction);
  func->SetParam(0,"W", Wl);
  func->SetParam(1,"Q2",Q2l);
  double xsec = fIntegrator->Integrate(*func);

/*
  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);
  LOG("DISXSec", pINFO)  
            << "x integration range = [" << xl.min << ", " << xl.max << "]";
  LOG("DISXSec", pINFO)  
            << "y integration range = [" << yl.min << ", " << yl.max << "]";

  GXSecFunc * func = new Integrand_D2XSec_DxDy_E(model, interaction);
  func->SetParam(0,"x",xl);
  func->SetParam(1,"y",yl);
  double xsec = fIntegrator->Integrate(*func);
*/
  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfHitNucRest);
  LOG("DISXSec", pINFO)  << "XSec[DIS] (E = " << Ev << " GeV) = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
void DISXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISXSec::LoadConfig(void)
{
  //-- get specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
  assert(fIntegrator);
}
//____________________________________________________________________________

