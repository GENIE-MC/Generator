//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 19, 2009 - CA
   First included in v2.5.1.
 @ Nov 26, 2009 - CA
   Fix mistake in convertion from dsigma/dOmega --> dsigma/dQ2 uncovered at
   the first comparison against electron QE data.

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/QuasiElastic/XSection/RosenbluthPXSec.h"
#include "Physics/QuasiElastic/XSection/ELFormFactorsModelI.h"
#include "Physics/QuasiElastic/XSection/TransverseEnhancementFFModel.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
RosenbluthPXSec::RosenbluthPXSec() :
XSecAlgorithmI("genie::RosenbluthPXSec")
{

}
//____________________________________________________________________________
RosenbluthPXSec::RosenbluthPXSec(string config) :
XSecAlgorithmI("genie::RosenbluthPXSec", config)
{

}
//____________________________________________________________________________
RosenbluthPXSec::~RosenbluthPXSec()
{
  if (fCleanUpfElFFModel) {
    delete fElFFModel;
  }
}
//____________________________________________________________________________
double RosenbluthPXSec::XSec(
    const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get interaction information
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  int nucpdgc = target.HitNucPdg();
  double E  = init_state.ProbeE(kRfHitNucRest);
  double Q2 = kinematics.Q2();
  double M  = target.HitNucMass();

  double E2 = E*E;
  double E3 = E*E2;
  double M2 = M*M;

  // Calculate scattering angle
  //
  // Q^2 = 4 * E^2 * sin^2 (theta/2) / ( 1 + 2 * (E/M) * sin^2(theta/2) ) =>
  // sin^2 (theta/2) = MQ^2 / (4ME^2 - 2EQ^2)

  double sin2_halftheta = M*Q2 / (4*M*E2 - 2*E*Q2);
  double sin4_halftheta = TMath::Power(sin2_halftheta, 2.);
  double cos2_halftheta = 1.-sin2_halftheta;
  //unused double cos_halftheta  = TMath::Sqrt(cos2_halftheta);
  double tan2_halftheta = sin2_halftheta/cos2_halftheta;

  // Calculate the elastic nucleon form factors
  fELFF.Calculate(interaction);
  double Gm  = pdg::IsProton(nucpdgc) ? fELFF.Gmp() : fELFF.Gmn();
  double Ge  = pdg::IsProton(nucpdgc) ? fELFF.Gep() : fELFF.Gen();
  double Ge2 = Ge*Ge;
  double Gm2 = Gm*Gm;

  // Calculate tau and the virtual photon polarization (epsilon)
  double tau     = Q2/(4*M2);
  double epsilon = 1. / (1. + 2.*(1.+tau)*tan2_halftheta);

  // Calculate the scattered lepton energy
  double Ep  = E / (1. + 2.*(E/M)*sin2_halftheta);
  //double Ep2 = Ep*Ep; // unused variable

  // Calculate the Mott cross section dsigma/dOmega
  double xsec_mott = (0.25 * kAem2 * Ep / E3) * (cos2_halftheta/sin4_halftheta);

  // Calculate the electron-nucleon elastic cross section dsigma/dOmega
  double xsec = xsec_mott * (Ge2 + (tau/epsilon)*Gm2) / (1+tau);

  // Convert dsigma/dOmega --> dsigma/dQ2
  xsec *= (kPi/(Ep*E));

  // The algorithm computes dxsec/dQ2
  // Check whether variable tranformation is needed
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("Rosenbluth", pDEBUG)
       << "Jacobian for transformation to: "
       << KinePhaseSpace::AsString(kps) << ", J = " << J;
#endif
    xsec *= J;
  }

  // If requested, return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Take into account the number of nucleons/tgt
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();
  xsec *= NNucl;

  // Compute & apply nuclear suppression factor
  // (R(Q2) is adapted from NeuGEN - see comments therein)
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);
  xsec *= R;

  return xsec;
}
//____________________________________________________________________________
double RosenbluthPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool RosenbluthPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsQuasiElastic()) return false;
  if(!proc_info.IsEM()) return false;

  const InitialState & init_state = interaction->InitState();

  int  hitnuc = init_state.Tgt().HitNucPdg();
  bool is_pn = (pdg::IsProton(hitnuc) || pdg::IsNeutron(hitnuc));
  if (!is_pn) return false;

  int  probe   = init_state.ProbePdg();
  bool is_chgl = pdg::IsChargedLepton(probe);
  if (!is_chgl) return false;

  return true;
}
//____________________________________________________________________________
void RosenbluthPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RosenbluthPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RosenbluthPXSec::LoadConfig(void)
{
  fElFFModel = 0;

  // load elastic form factors model
  fElFFModel = dynamic_cast<const ELFormFactorsModelI *> ( this -> SubAlg("ElasticFormFactorsModel" ) ) ;

  assert(fElFFModel);

  fCleanUpfElFFModel = false;
  bool useFFTE = false ;
  GetParam( "UseElFFTransverseEnhancement", useFFTE ) ;
  if( useFFTE ) {
    const ELFormFactorsModelI* sub_alg = fElFFModel;
    fElFFModel = dynamic_cast<const ELFormFactorsModelI *> ( this -> SubAlg("TransverseEnhancement") ) ;
    dynamic_cast<const TransverseEnhancementFFModel*>(fElFFModel)->SetElFFBaseModel( sub_alg );
    fCleanUpfElFFModel = true;
  }
  fELFF.SetModel(fElFFModel);

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
