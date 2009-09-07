//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 18, 2008 - CA
   Add protection against unphysical Q2 limits
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/QELXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/Range1.h"
#include "Utils/GSLUtils.h"

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
  if(rQ2.min<0 || rQ2.max<0) return 0;
  LOG("QELXSec", pDEBUG) 
          << "Q2 integration range = (" << rQ2.min << ", " << rQ2.max << ")";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

#ifdef __GENIE_GSL_ENABLED__
  ROOT::Math::IBaseFunctionOneDim * func = new 
      utils::gsl::wrap::dXSec_dQ2_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type = 
      utils::gsl::Integration1DimTypeFromString(fGSLIntgType);     
  ROOT::Math::Integrator ig(ig_type);
  ig.SetFunction(*func);
  ig.SetRelTolerance(fGSLRelTol);
  double xsec = ig.Integral(rQ2.min, rQ2.max) * (1E-38 * units::cm2);
     
#else
  GXSecFunc * func = new Integrand_DXSec_DQ2_E(model, interaction);
  func->SetParam(0,"Q2",rQ2);
  double xsec = fIntegrator->Integrate(*func);

#endif

  //LOG("QELXSec", pDEBUG) << "XSec[QEL] (E = " << E << ") = " << xsec;

  delete func;
  delete interaction;

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
  // Get the specified GENIE integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *>(this->SubAlg("Integrator"));
  assert(fIntegrator);

  // Get GSL integration type & relative tolerance
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type",  "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance", 0.01);
}
//____________________________________________________________________________

