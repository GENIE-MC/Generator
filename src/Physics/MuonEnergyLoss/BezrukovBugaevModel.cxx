//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - December 10, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/Integrator.h>

#include "Physics/MuonEnergyLoss/BezrukovBugaevModel.h"
//#include "Numerical/IntegratorI.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::mueloss;
using namespace genie::constants;

//____________________________________________________________________________
BezrukovBugaevModel::BezrukovBugaevModel() :
MuELossI("genie::mueloss::BezrukovBugaevModel")
{

}
//____________________________________________________________________________
BezrukovBugaevModel::BezrukovBugaevModel(string config) :
MuELossI("genie::mueloss::BezrukovBugaevModel", config)
{

}
//____________________________________________________________________________
BezrukovBugaevModel::~BezrukovBugaevModel()
{

}
//____________________________________________________________________________
double BezrukovBugaevModel::dE_dx(double E, MuELMaterial_t material) const
{
// Calculate the muon -dE/dx due to muon nuclear interaction (in GeV^-2).
// To convert the result to more handly units, eg MeV/(gr/cm^2), just write:
// dE_dx /= (units::MeV/(units::g/units::cm2));

  if(material == eMuUndefined) return 0;
  if(E<=MuELProcess::Threshold(this->Process()) || E>=kMaxMuE) return 0;

  // material Z,E
  double Z = MuELMaterial::Z(material);
  double A = MuELMaterial::A(material);

  // calculate (the min,max) fraction of energy, v,  carried to the photon
  double Vmin = 0.;
  double Vmax = 1. - 0.75*kSqrtNapierConst* (kMuonMass/E) * TMath::Power(Z,1/3.);

  // integrate the Bezrukov-Bugaev differential cross section v*ds/dv for
  // muon nuclear interaction over v

  ROOT::Math::IBaseFunctionOneDim * integrand = 
          new gsl::BezrukovBugaevIntegrand(E,A);
  ROOT::Math::IntegrationOneDim::Type ig_type =
          utils::gsl::Integration1DimTypeFromString("adaptive");

  double abstol   = 1;    // We mostly care about relative tolerance
  double reltol   = 1E-4;
  int    nmaxeval = 100000;
  ROOT::Math::Integrator ig(*integrand,ig_type,abstol,reltol,nmaxeval);

  // calculate the b factor (bE = -dE/dx) in GeV^-3
  A *= units::g;
  double bnucl = (kNA/A) * ig.Integral(Vmin, Vmax); 

  delete integrand;

  // calculate the dE/dx due to muon nuclear interaction in GeV^-2
  double de_dx = bnucl*E;
  return de_dx;
}
//____________________________________________________________________________
//void BezrukovBugaevModel::Configure(const Registry & config)
//{
//  Algorithm::Configure(config);
//  this->LoadConfig();
//}
////____________________________________________________________________________
//void BezrukovBugaevModel::Configure(string config)
//{
//  Algorithm::Configure(config);
//  this->LoadConfig();
//}
////____________________________________________________________________________
//void BezrukovBugaevModel::LoadConfig(void)
//{
//  //fIntegrator = 
////       dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
////  assert(fIntegrator);
//}
//____________________________________________________________________________
//
// BezrukovBugaevIntegrand
//
//____________________________________________________________________________
gsl::BezrukovBugaevIntegrand::BezrukovBugaevIntegrand(double E, double A) :
ROOT::Math::IBaseFunctionOneDim()
{
  fE = E;
  fA = A;
}
//____________________________________________________________________________
gsl::BezrukovBugaevIntegrand::~BezrukovBugaevIntegrand()
{

}
//____________________________________________________________________________
unsigned int gsl::BezrukovBugaevIntegrand::NDim(void) const
{
  return 1;
}
//____________________________________________________________________________
double gsl::BezrukovBugaevIntegrand::DoEval(double xin) const
{
// Returns v*(ds/dv)

  double v = xin; // v, the fraction of energy transfered to the photon

  if (! (v >0)) return 0;
  if (   v >1)  return 0;
  if (! (fE>0)) return 0;

  double a    = kAem;
  double pi   = kPi;
  double mmu2 = kMuonMass2;
  double v2   = TMath::Power(v,2.);
  double t    = mmu2 *v2/(1-v);
  double k    = 1. - 2./v + 2./v2;
  double A13  = TMath::Power(fA,1./3.);
  double M1_2 = 0.54; // m1^2 in photonuclear diff. xsec formula (in GeV^2)
  double M2_2 = 1.80; // m2^2 in photonuclear diff. xsec formula (in GeV^2)
  double M1_2_t = M1_2 / t;
  double M2_2_t = M2_2 / t;
  double mmu2_t = mmu2 / t;
  double d      = M1_2 / (t + M1_2);

  // Calculate the cross section (in ub) for photonuclear interaction
  double Ep   = v*fE; // photon energy (GeV)
  double loge = TMath::Log(0.0213*Ep); // factor 0.0213 has units of GeV^-1
  double sig  = 114.3 + 1.647 * loge*loge; // in ub

  // Calculate the factor G (dimensionless) appearing in the differential
  // photonuclear interaction cross section
  double x   = 0.00282*A13*sig;  // factor 0.00282 has units of ub^-1
  double x2  = x*x;
  double x3  = x2*x;
  double G   = 3*(0.5*x2 - 1. + (1.+x)*TMath::Exp(-x)) /x3;

  // Calculate the differential cross section ds/dv for muon nuclear
  // interaction based on the Bezrukov-Bugaev formula.
  double bbA = 0.5*(a/pi) * fA * sig * v;
  double bbB = 0.75*G     * ( k*TMath::Log(1.+M1_2_t) - k*d - 2.*mmu2_t );
  double bbC = 0.25       * ( k*TMath::Log(1.+M2_2_t) - 2.*mmu2_t );
  double bbD = 0.5*mmu2_t * ( 0.75*G*d + 0.25*M2_2_t*TMath::Log(1.+1./M2_2_t) );

  double ds_dv  = bbA*(bbB+bbC+bbD); // in um (microbarns)
  double vds_dv = v*ds_dv;

  vds_dv *= units::ub; // ub -> GeV^-2
  return vds_dv;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
  gsl::BezrukovBugaevIntegrand::Clone(void) const
{
  return new gsl::BezrukovBugaevIntegrand(fE, fA);
}
//____________________________________________________________________________
