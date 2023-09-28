//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Brinden Carlson bcarlson1@ufl.edu
 University of Florida
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/NuElectron/XSection/PXSecOnElectron.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/ElectronVelocity.h"
#include "Framework/Numerical/MathUtils.h"

#include "TFile.h"
#include "TGraph.h"
#include <fstream>
#include <iterator>

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
PXSecOnElectron::PXSecOnElectron(string name, string config) :
XSecAlgorithmI(name, config)
{

}
//____________________________________________________________________________
PXSecOnElectron::PXSecOnElectron() :
XSecAlgorithmI("genie::PXSecOnElectron")
{

}
//____________________________________________________________________________
PXSecOnElectron::~PXSecOnElectron()
{

}
//____________________________________________________________________________
double PXSecOnElectron::Integral(const Interaction * interaction) const
{
  double xsec_sum = 0; 
  double xsec_sum2 = 0; 
  int NInt = 0; //Count number of integrations
  do{
    NInt++;
    Interaction in_curr(*interaction); //Copy interaction object
    fElectronVelocity->InitializeVelocity(in_curr); //Modify interaction to give electron random velocity from selected distribution
    double xsec = fXSecIntegrator->Integrate(this,&in_curr); // In ele at rest
    //get beta - comps orthogonal to beam x,y
    //scale = sqrt(1-b_t^2)
    //TVector3 beta = in_curr.InitState().Tgt().HitPartP4().BoostVector(); // beta
    //double beta_tangent = sqrt(TMath::Power(beta[0],2)+TMath::Power(beta[1],2)); //Component tangential to beam
    
    auto probe_direction = in_curr.InitState().GetProbeP4(kRfHitElRest)->Vect().Unit();
    auto electron_mom = in_curr.InitState().Tgt().HitPartP4();
    auto transverse_mom = utils::math::GetOrthogonal(electron_mom, probe_direction);
    auto beta_transverse = transverse_mom.BoostVector().Mag();
    //double beta_transverse = beta_transverse.Mag();

    xsec *= sqrt(1-TMath::Power(beta_transverse,2)); //Correct for lorentz factor
    xsec_sum+=xsec;
    xsec_sum2+=TMath::Power(xsec,2);
    double xsec_mean = xsec_sum/NInt;
    //var = (sum(xi^2)/N-xsec_mean^2)
    //rel_err = sigma/sqrt(n)*mean
    double xsec_err = sqrt((xsec_sum2/NInt-TMath::Power(xsec_mean,2))/NInt)/xsec_mean;
    if (NInt > 1 && xsec_err < fErrTolerance){break;} //Break condition for dipping below set tolerance
  }
  while ( NInt < fNIntegration); 
  double xsec_avg = xsec_sum/NInt;
  return xsec_avg;
}
//____________________________________________________________________________
bool PXSecOnElectron::ValidProcess(const Interaction * interaction) const
{
  //if(pdg::IsElectron(interaction->InitState().ProbePdg())) return true;
  if(interaction->TestBit(kISkipProcessChk)) return true;
  return true;
}
//____________________________________________________________________________
bool PXSecOnElectron::ValidKinematics(const Interaction* interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;
  return true;
}
//____________________________________________________________________________
void PXSecOnElectron::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PXSecOnElectron::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PXSecOnElectron::LoadConfig(void)
{
  //Integration parameters
  GetParam( "N-Integration-Samples", fNIntegration ) ; //
  GetParam( "NuE-XSecRelError" , fErrTolerance ) ; //

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator")); //
  fElectronVelocity =
      dynamic_cast<const ElectronVelocity *> (this->SubAlg("Electron-Velocity")); //
  if (!fElectronVelocity) {
    LOG("PXSecOnElectron", pDEBUG)
    << "fElectronVelocity is not initialized correctly.";
    exit(78);
  }
}
//____________________________________________________________________________
