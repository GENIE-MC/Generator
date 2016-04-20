//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 04, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 18, 2008 - CA
   Add protection against unphysical Q2 limits
 @ Sep 07, 2009 - CA
   Integrated with GNU Numerical Library (GSL) via ROOT's MathMore library.
 @ Nov 30, 2009 - CA
   Added option to integrate over the hit nucleon momentum distribution.
 @ Mar 18, 2016 - JJ (SD)
   Moved generation of the struck nucleon position and momentum to the
   QEL implementations of the XSecAlgorithmI class. Integrate() now contains
   the code that was previously in IntegrateOnce().
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/Integrator.h>

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Conventions/KineVar.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/QELXSec.h"
#include "CrossSections/GSLXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"
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
  LOG("QELXSec",pDEBUG) << "Beginning integrate";
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

  ROOT::Math::IBaseFunctionOneDim * func = new 
      utils::gsl::dXSec_dQ2_E(model, interaction);
  ROOT::Math::IntegrationOneDim::Type ig_type = 
      utils::gsl::Integration1DimTypeFromString(fGSLIntgType);
  
  double abstol = 1; //We mostly care about relative tolerance
  ROOT::Math::Integrator ig(*func,ig_type,abstol,fGSLRelTol,fGSLMaxEval);
  double xsec = ig.Integral(rQ2.min, rQ2.max) * (1E-38 * units::cm2);
     
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
  // Get GSL integration type & relative tolerance
  fGSLIntgType = fConfig->GetStringDef("gsl-integration-type", "adaptive");
  fGSLRelTol   = fConfig->GetDoubleDef("gsl-relative-tolerance", 0.01);
  fGSLMaxEval  = (unsigned int) fConfig->GetIntDef ("gsl-max-eval", 100000);
}
//____________________________________________________________________________

