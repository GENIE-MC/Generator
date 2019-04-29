//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/IntegratorMultiDim.h>

#include "Physics/MuonEnergyLoss/KokoulinPetrukhinModel.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;
using namespace genie::mueloss;
using namespace genie::constants;

//____________________________________________________________________________
KokoulinPetrukhinModel::KokoulinPetrukhinModel() :
MuELossI("genie::mueloss::KokoulinPetrukhinModel")
{

}
//____________________________________________________________________________
KokoulinPetrukhinModel::KokoulinPetrukhinModel(string config) :
MuELossI("genie::mueloss::KokoulinPetrukhinModel", config)
{

}
//____________________________________________________________________________
KokoulinPetrukhinModel::~KokoulinPetrukhinModel()
{

}
//____________________________________________________________________________
double KokoulinPetrukhinModel::dE_dx(double E, MuELMaterial_t material) const
{
// Calculate the muon -dE/dx due to e+e- pair production (in GeV^-2).
// To convert the result to more handly units, eg MeV/(gr/cm^2), just write:
// dE_dx /= (units::MeV/(units::g/units::cm2));

  if(material == eMuUndefined) return 0;
  if(E<=MuELProcess::Threshold(this->Process()) || E>=kMaxMuE) return 0;

  // material Z,E
  double Z = MuELMaterial::Z(material);
  double A = MuELMaterial::A(material);

  // calculate (the min,max) fraction of energy, v, carried to the photon
  double Vmin = 4.*kElectronMass/E;
  double Vmax = 1. - 0.75*kSqrtNapierConst* (kMuonMass/E) * TMath::Power(Z,1/3.);

  // claculate the limits of the asymmetry parameter of the e+e- pair
  // p = (E(+) - E(-)) / (E(+) + E(-))
  double Pmin = 0.;
  double Pmax = 1.;

  // integrate the Kokulin-Petrukhin differential cross section v*ds/dv for
  // muon e+e- pair production over v and p

  ROOT::Math::IBaseFunctionMultiDim * integrand =
      new gsl::KokoulinPetrukhinIntegrand(E, Z);
  ROOT::Math::IntegrationMultiDim::Type ig_type =
      utils::gsl::IntegrationNDimTypeFromString("adaptive");
  
  double abstol   = 1; // We mostly care about relative tolerance.
  double reltol   = 1E-4;
  int    nmaxeval = 500000;

  ROOT::Math::IntegratorMultiDim ig(*integrand, ig_type, abstol, reltol, nmaxeval);

  double min[2] = { Vmin, Pmin };
  double max[2] = { Vmax, Pmax };

  // calculate the b factor (bE = -dE/dx) in GeV^-3
  A *= units::g;
  double bpair = (2*kNA/A) * ig.Integral(min, max);

  delete integrand;

  // calculate the dE/dx due to muon bremsstrahlung in GeV^-2
  double de_dx = bpair*E;
  return de_dx;
}
//____________________________________________________________________________
//void KokoulinPetrukhinModel::Configure(const Registry & config)
//{
//  Algorithm::Configure(config);
//  this->LoadConfig();
//}
////____________________________________________________________________________
//void KokoulinPetrukhinModel::Configure(string config)
//{
//  Algorithm::Configure(config);
//  this->LoadConfig();
//}
////____________________________________________________________________________
//void KokoulinPetrukhinModel::LoadConfig(void)
//{
////  fIntegrator = dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));
////  assert(fIntegrator);
//}
//____________________________________________________________________________
//
// KokoulinPetrukhinIntegrand
//
//____________________________________________________________________________
gsl::KokoulinPetrukhinIntegrand::KokoulinPetrukhinIntegrand(double E, double Z) :
ROOT::Math::IBaseFunctionMultiDim()
{
  fE = E;
  fZ = Z;
}
//____________________________________________________________________________
gsl::KokoulinPetrukhinIntegrand::~KokoulinPetrukhinIntegrand()
{

}
//____________________________________________________________________________
unsigned int gsl::KokoulinPetrukhinIntegrand::NDim(void) const
{
  return 2;
} 
//____________________________________________________________________________
double gsl::KokoulinPetrukhinIntegrand::DoEval (const double * xin) const
{
// Returns v*(d^2s/dvdp)

  double v  = xin[0]; // v, the fraction of energy transfered to the photon
  double p  = xin[1]; //

  if (! (v >0)) return 0;
  if (   v >1)  return 0;
  if (! (fE>0)) return 0;

  double pmax_v = (1. - 6.*kMuonMass2 / (fE*fE*(1.-v)) ) * 
                                 TMath::Sqrt(1.-4.*kElectronMass/(fE*v));
  if(p>pmax_v) return 0;

  double v2      = TMath::Power(v,2.);
  double p2      = TMath::Power(p,2.);
  double R       = 189.;
  double a4      = TMath::Power(kAem,4.);
  double pi      = kPi;
  double ZLe2    = TMath::Power(fZ*kLe,2.);
  double Zm13    = TMath::Power(fZ,-1./3.);
  double Zm23    = TMath::Power(fZ,-2./3.);
  double Z13     = TMath::Power(fZ,1./3.);
  double me      = kElectronMass;
  double me2     = kElectronMass2;
  double mmu     = kMuonMass;
  double mmu2    = kMuonMass2;
  double memu2   = me2/mmu2;
  double memu    = me/mmu;
  double mume    = mmu/me;

  double b    = 0.5*v2/(1.-v);
  double xi   = (1.-p2) * TMath::Power(0.5*v*mume, 2.) / (1.-v);

  // Approximate the FIm factor (dimensionless) for the Kokoulin Petrukhin
  // pair production cross section

  double FImA = ( (1.+p2)*(1.+3.*b/2.) - (1.+2.*b)*(1.-p2)/xi ) * TMath::Log(1.+xi);
  double FImB = xi*(1.-p2-b) / (1.+xi);
  double FImC = (1.+2.*b) * (1.-p2);
  double Ym   = (4. + p2 + 3.*b*(1.+p2)) /
                      ((1.+p2)*(1.5+2.*b)*TMath::Log(3.+xi) + 1. - 1.5*p2);
  double LmA  = (2./3.) * mume * R * Zm23;
  double LmB  = 1. + (2.*me*R * kSqrtNapierConst * Zm13 * (1+xi) * (1+Ym)) / (fE*v*(1-p2) );
  double Lm   = TMath::Log(LmA/LmB);
  double FIm  = (FImA+FImB+FImC)*Lm;

  // Approximate the FIe factor (dimensionless) for the Kokoulin Petrukhin
  // pair production cross section

  double FIeA = ( (2.+p2) * (1.+b) + xi*(3.+p2) ) * log(1.+1./xi);
  double FIeB = (1.-p2-b) / (1.+xi);
  double FIeC = 3. + p2;
  double Ye   = (5. - p2 + 4*b*(1+p2)) /
                    (2.*(1.+3.*b)*TMath::Log(3.+1./xi) - p2 - 2.*b*(2.-p2));
  double x_Y  = (1+xi)*(1+Ye);
  double LeA  = R*Zm13*TMath::Sqrt(x_Y);
  double LeB  = 1. + (2.*me*R*kSqrtNapierConst*Zm13*x_Y) / (fE*v*(1-p2));
  double LeC  = 1.5 * memu * Z13;
  double LeD  = 1 + TMath::Power(LeC,2.)*x_Y;
  double Le   = TMath::Log(LeA/LeB) - 0.5*TMath::Log(LeD);
  double FIe  = (FIeA+FIeB-FIeC)*Le;

  double d2s_dvdp = a4 * (2./(3.*pi)) * ZLe2 * ((1.-v)/v) * (FIe + FIm*memu2);

  double vd2s_dvdp = v*d2s_dvdp;
  return vd2s_dvdp;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim * 
  gsl::KokoulinPetrukhinIntegrand::Clone(void) const
{
  return new gsl::KokoulinPetrukhinIntegrand(fE,fZ);
}
//____________________________________________________________________________

