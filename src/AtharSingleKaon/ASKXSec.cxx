//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/IntegratorMultiDim.h>
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "AtharSingleKaon/ASKXSec.h"
#include "AtharSingleKaon/ASKGSL.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
//#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"
#include "Utils/GSLUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
ASKXSec::ASKXSec() :
XSecIntegratorI("genie::ASKXSec")
{

}
//____________________________________________________________________________
ASKXSec::ASKXSec(string config) :
XSecIntegratorI("genie::ASKXSec", config)
{

}
//____________________________________________________________________________
ASKXSec::~ASKXSec()
{

}
//____________________________________________________________________________
double ASKXSec::Integrate(
      const XSecAlgorithmI * model, const Interaction * in) const
{

  LOG( "ASKXSec", pDEBUG ) << "Thank you for calling ASKXSec::Integrate";

  #ifndef __GENIE_GSL_ENABLED__
    SLOG("ASKXSec",pFATAL)<<"Trying to call AtharSingleKaon without GSL";
    return 0.;
  #endif
  
  const InitialState & init_state = in -> InitState();
  
  if(! model->ValidProcess(in) ) return 0.;
  
  const KPhaseSpace & kps = in->PhaseSpace(); // only OK phase space for this
  if(!kps.IsAboveThreshold()) {
     LOG("ASKXSec", pDEBUG)  << "*** Below energy threshold";
     return 0;
  }

  // If the input interaction is off a nuclear target, then chek whether 
  // the corresponding free nucleon cross section already exists at the 
  // cross section spline list. 
  // Cross section for PP scales with number of protons, NP and NN scale 
  // with number of neutrons
  int nucpdgc = init_state.Tgt().HitNucPdg();
  int NNucl   = (pdg::IsProton(nucpdgc)) ? 
                   init_state.Tgt().Z() : init_state.Tgt().N();
  double Ev = init_state.ProbeE(kRfHitNucRest);

  XSecSplineList * xsl = XSecSplineList::Instance();
  if(init_state.Tgt().IsNucleus() && !xsl->IsEmpty() ) {
    Interaction * interaction = new Interaction(*in);
    Target * target = interaction->InitStatePtr()->TgtPtr();
    if(pdg::IsProton(nucpdgc)) { target->SetId(kPdgTgtFreeP); }
    else                       { target->SetId(kPdgTgtFreeN); }
    if(xsl->SplineExists(model,interaction)) {
      const Spline * spl = xsl->GetSpline(model, interaction);
      double xsec = spl->Evaluate(Ev);
      LOG("ASKXSec", pINFO)  
        << "From XSecSplineList: XSec[ASK,free nucleon] (E = " << Ev << " GeV) = " << xsec;
      if(! interaction->TestBit(kIAssumeFreeNucleon) ) { 
          xsec *= NNucl; 
          LOG("ASKXSec", pINFO)  << "XSec[ASK] (E = " << Ev << " GeV) = " << xsec;
      }
      delete interaction;
      return xsec;
    }
    delete interaction;
  }

  // no free nucelon spline exists -- do the integration
 
  // Check this
  double Enu = init_state.ProbeE(kRfLab);
  int kpdg = in->ExclTag().StrangeHadronPdg();
  double mk   = PDGLibrary::Instance()->Find(kpdg)->Mass();
  double ml   = PDGLibrary::Instance()->Find(in->FSPrimLeptonPdg())->Mass();

  // integration bounds for T (kinetic energy)
  double zero    = 0.0;
  double tmax = Enu - mk - ml;

  LOG("ASKXSec", pDEBUG)
       << "Lepton/Kaon KE integration range = [" << 0.0 << ", " << tmax << "]";

  Interaction * interaction = new Interaction(*in);
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);
  
  double xsec = 0;

  // do the integration over log(1-costheta) so it's not so sharply peaked
    
  ROOT::Math::IBaseFunctionMultiDim * func = 
        new utils::gsl::d3Xsec_dTldTkdCosThetal(model, interaction);
  double kine_min[3] = { zero, zero, -20 }; // Tlep, Tkaon, cosine theta lep
  double kine_max[3] = { tmax, tmax,  0.69314718056 }; // Tlep, Tkaon, cosine theta lep
    
  ROOT::Math::IntegrationMultiDim::Type ig_type = 
    utils::gsl::IntegrationNDimTypeFromString(fGSLIntgType);
        
  double abstol = 1; //We mostly care about relative tolerance.
  ROOT::Math::IntegratorMultiDim ig(*func, ig_type, abstol, fGSLRelTol, fGSLMaxEval);
  
  xsec = ig.Integral(kine_min, kine_max) * (1E-38 * units::cm2);
  delete func;

  delete interaction;

  return xsec;
}
//____________________________________________________________________________
void ASKXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ASKXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ASKXSec::LoadConfig(void)
{
  // Get specified GENIE integration algorithm
//  fIntegrator = dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
//  assert(fIntegrator);

  // Get GSL integration type & relative tolerance
  fGSLIntgType   = fConfig->GetStringDef("gsl-integration-type" ,    "vegas");
  fGSLMaxEval    = (unsigned int) fConfig->GetIntDef("gsl-max-evals", 20000);
  fGSLRelTol     = fConfig->GetDoubleDef("gsl-relative-tolerance",    0.01);
  fSplitIntegral = fConfig->GetBoolDef("split-integral",              true);
}
//_____________________________________________________________________________


