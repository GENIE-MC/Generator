//____________________________________________________________________________
/*!

\class    genie::NuElectronXSec

\brief    nu/nubar + e- scattering cross section. Integrates the loaded
          differential cross section model. An analytical cross section
          model also exists, so you cal also use that if you do not apply
          any kinematical cuts.

          The cross section algorithm handles:
             - nue/nuebar + e- -> nue/nuebar + e- [CC + NC + interference]
             - numu/nutau + e- -> numu/nutau + e- [NC]
             - numubar/nutaubar + e- -> numubar/nutaubar + e- [NC]
             - numu/nutau + e- -> l- + nu_e [CC]

          NuElectronXSec is a concrete implementation of the XSecAlgorithmI
          interface. \n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  February 10, 2006

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/NuElectronXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
NuElectronXSec::NuElectronXSec() :
XSecAlgorithmI("genie::NuElectronXSec")
{

}
//____________________________________________________________________________
NuElectronXSec::NuElectronXSec(string config) :
XSecAlgorithmI("genie::NuElectronXSec", config)
{

}
//____________________________________________________________________________
NuElectronXSec::~NuElectronXSec()
{

}
//____________________________________________________________________________
double NuElectronXSec::XSec(const Interaction * interaction) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get initial & final state information
  const InitialState & init_state = interaction->GetInitialState();
  double E  = init_state.GetProbeE(kRfLab);

  double e    = 1E-6;
  double ymin = e;
  double ymax = 1-e;

  GXSecFunc * func = new Integrand_DXSec_Dy_E(fDiffXSecModel, interaction);
  func->SetParam(0,"y",ymin,ymax);
  double xsec = fIntegrator->Integrate(*func);

  LOG("Elastic", pDEBUG) << "*** XSec[ve-] (E=" << E << ") = " << xsec;

  delete func;
  return xsec;
}
//____________________________________________________________________________
bool NuElectronXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  return fDiffXSecModel->ValidProcess(interaction);
}
//____________________________________________________________________________
bool NuElectronXSec::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  return true;
}
//____________________________________________________________________________
void NuElectronXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NuElectronXSec::LoadConfig(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed and sets configuration variables to avoid looking up at the Registry
// all the time.

  //-- get an algorithm to calculate differential cross sections dxsec/dy
  fDiffXSecModel =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fDiffXSecModel);

  //-- get the specified integration algorithm
  fIntegrator = dynamic_cast<const IntegratorI *> (
                 this->SubAlg("integrator-alg-name", "integrator-param-set"));
  assert(fIntegrator);
}
//____________________________________________________________________________

