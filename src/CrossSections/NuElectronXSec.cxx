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

#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "CrossSections/NuElectronXSec.h"
#include "Messenger/Messenger.h"
#include "Numerical/UnifGrid.h"
#include "Numerical/FunctionMap.h"
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

  // Integration grid
  double e    = 1E-6;
  double ymin = e;
  double ymax = 1-e;

  const double logymin = TMath::Log(ymin);
  const double logymax = TMath::Log(ymax);
  const double dlogy   = (logymax-logymin)/(fNBins-1);

  UnifGrid grid;
  grid.AddDimension(fNBins, logymin, logymax);

  FunctionMap fmap(grid);

  // Compute the differential cross section over the allowed phase space
  for(int i = 0; i < fNBins; i++) {
    double y = TMath::Exp(logymin + i * dlogy);

    //-- update the scattering parameters
    interaction->GetKinematicsPtr()->Sety(y);
    //-- compute dsec/dy
    double dsig_dy  = fDiffXSecModel->XSec(interaction);
    //-- push y*(dxsec/dy) to the FunctionMap
    fmap.AddPoint(y*dsig_dy, i);
  }

  // Do the numerical integration
  double xsec = fIntegrator->Integrate(fmap);
  LOG("Elastic", pDEBUG) << "*** XSec[ve-] (E=" << E << ") = " << xsec;
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

  //----- Get an algorithm to calculate differential cross sections dxsec/dy
  fDiffXSecModel =
          dynamic_cast<const XSecAlgorithmI *> (this->SubAlg(
                         "partial-xsec-alg-name", "partial-xsec-param-set"));
  assert(fDiffXSecModel);

  //----- Get an integrator
  string intgr = fConfig->GetStringDef("integrator", "genie::Simpson1D");

  AlgFactory * algf = AlgFactory::Instance();
  fIntegrator = dynamic_cast<const IntegratorI *>(algf->GetAlgorithm(intgr));

  assert(fIntegrator);

  //----- get the input number of dxsec/dlogQ2 points for num. integration
  //      or use a default if no number is specified
  //      (must be odd number for the Simpson rule)
  fNBins = fConfig->GetIntDef("N-logy-bins", 131);
  LOG("Elastic", pDEBUG) << "Number of integration (logy) bins = " << fNBins;
  assert(fNBins>2);
}
//____________________________________________________________________________

