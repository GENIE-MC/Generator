//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
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
#include "Conventions/KineVar.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/QELXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
QELXSec::QELXSec() :
XSecIntegratorI("genie::QELXSec")
{

}
//____________________________________________________________________________
QELXSec::QELXSec(string config) :
XSecIntegratorI("genie::QELXSec", config)
{

}
//____________________________________________________________________________
QELXSec::~QELXSec()
{

}
//____________________________________________________________________________
double QELXSec::Integrate(
                  const XSecAlgorithmI * model, const Interaction * in) const
{
  if(! model->ValidProcess(in)) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("QELXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }
  Range1D_t rQ2 = kps.Limits(kKVQ2);
  LOG("QELXSec", pDEBUG) 
          << "Q2 integration range = (" << rQ2.min << ", " << rQ2.max << ")";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  GXSecFunc * func = new Integrand_DXSec_DQ2_E(model, interaction);
  func->SetParam(0,"Q2",rQ2);
  double xsec = fIntegrator->Integrate(*func);

  //LOG("QELXSec", pDEBUG) << "XSec[QEL] (E = " << E << ") = " << xsec;

  delete interaction;
  delete func;
  return xsec;
}
//____________________________________________________________________________
void QELXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELXSec::LoadConfig(void)
{
  //-- get the specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *>(this->SubAlg("Integrator"));
  assert(fIntegrator);
}
//____________________________________________________________________________

