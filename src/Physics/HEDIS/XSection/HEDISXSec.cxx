//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/XSection/HEDISXSec.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Numerical/GSLUtils.h"

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
HEDISXSec::HEDISXSec() :
XSecIntegratorI("genie::HEDISXSec")
{

}
//____________________________________________________________________________
HEDISXSec::HEDISXSec(string config) :
XSecIntegratorI("genie::HEDISXSec", config)
{

}
//____________________________________________________________________________
HEDISXSec::~HEDISXSec()
{

}
//____________________________________________________________________________
double HEDISXSec::Integrate(
                 const XSecAlgorithmI * model, const Interaction * in) const
{

  if(! model->ValidProcess(in) ) return 0.;

  const KPhaseSpace & kps = in->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("DISXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  const InitialState & init_state = in->InitState();
  double Ev = init_state.ProbeE(kRfLab);
 
  // If the input interaction is off a nuclear target, then chek whether 
  // the corresponding free nucleon cross section already exists at the 
  // cross section spline list. 
  // If yes, calculate the nuclear cross section based on that value.
  //
  XSecSplineList * xsl = XSecSplineList::Instance();
  if(init_state.Tgt().IsNucleus() && !xsl->IsEmpty() ) {

    int nucpdgc = init_state.Tgt().HitNucPdg();
    int NNucl   = (pdg::IsProton(nucpdgc)) ? init_state.Tgt().Z() : init_state.Tgt().N();

    Interaction * interaction = new Interaction(*in);
    Target * target = interaction->InitStatePtr()->TgtPtr();
    if(pdg::IsProton(nucpdgc)) { target->SetId(kPdgTgtFreeP); }
    else                       { target->SetId(kPdgTgtFreeN); }
    if(xsl->SplineExists(model,interaction)) {
      const Spline * spl = xsl->GetSpline(model, interaction);
      double xsec = spl->Evaluate(Ev);
      LOG("HEDISXSec", pINFO) << "From XSecSplineList: XSec[HEDIS,free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if( !interaction->TestBit(kIAssumeFreeNucleon) ) { 
          xsec *= NNucl; 
          LOG("HEDISXSec", pINFO)  << "XSec[HEDIS] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
  }

  LOG("HEDISXSec", pINFO)  << in->AsString();

  Range1D_t xl  = kps.XLim();
  Range1D_t Q2l = kps.Q2Lim();
  LOG("HEDISXSec", pDEBUG) << "X only kinematic range = [" << xl.min << ", " << xl.max << "]";
  LOG("HEDISXSec", pDEBUG) << "Q2 only kinematic range = [" << Q2l.min << ", " << Q2l.max << "]";

  if (xl.min < fSFXmin)  xl.min=fSFXmin;
  if (Q2l.min < fSFQ2min) Q2l.min=fSFQ2min;
  if (Q2l.max > fSFQ2max) Q2l.max=fSFQ2max;

  LOG("HEDISXSec", pDEBUG) << "X kinematic+PDF range = [" << xl.min << ", " << xl.max << "]";
  LOG("HEDISXSec", pDEBUG) << "Q2 kinematic+PDF range = [" << Q2l.min << ", " << Q2l.max << "]";


  bool phsp_ok = 
      (Q2l.min >= 0. && Q2l.max >= 0. && Q2l.max >= Q2l.min &&
        xl.min >= 0. &&  xl.max <= 1. &&  xl.max >=  xl.min);

  if (!phsp_ok) return 0.;


  // Just go ahead and integrate the input differential cross section for the 
  // specified interaction.
  //
  double xsec = 0.;

  Interaction * interaction = new Interaction(*in);
  
  // If a GSL option has been chosen, then the total xsec is recomptued
  ROOT::Math::IBaseFunctionMultiDim * func = new utils::gsl::d2XSec_dlog10xdlog10Q2_E(model, interaction);
  ROOT::Math::IntegrationMultiDim::Type ig_type = utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
  double abstol = 0; //In vegas we care about number of interactions
  double reltol = 0; //In vegas we care about number of interactions
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, reltol, fGSLMaxEval);
  double kine_min[2] = { TMath::Log10(xl.min), TMath::Log10(Q2l.min) };
  double kine_max[2] = {TMath::Log10(xl.max), TMath::Log10(Q2l.max) };
  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
  delete func;

  delete interaction;

  LOG("HEDISXSec", pINFO)  << "XSec[HEDIS] (E = " << Ev << " GeV) = " << xsec * (1E+38/units::cm2) << " x 1E-38 cm^2";

  return xsec;

}
//____________________________________________________________________________
void HEDISXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISXSec::LoadConfig(void)
{

  // Get GSL integration type & relative tolerance
  GetParamDef( "gsl-integration-type", fGSLIntgType, string("vegas") ) ;

  int max_eval ;
  GetParamDef( "gsl-max-eval", max_eval, 500000 ) ;
  fGSLMaxEval  = (unsigned int) max_eval ;

  // Limits from the SF tables that are useful to reduce computation 
  // time of the total cross section
  GetParam("XGrid-Min",  fSFXmin ) ;
  GetParam("Q2Grid-Min", fSFQ2min ) ;
  GetParam("Q2Grid-Max", fSFQ2max ) ;

}
