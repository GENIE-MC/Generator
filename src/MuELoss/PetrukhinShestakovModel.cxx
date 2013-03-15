//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - December 10, 2003

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "MuELoss/PetrukhinShestakovModel.h"
#include "Numerical/IntegratorI.h"

using namespace genie;
using namespace genie::mueloss;
using namespace genie::constants;

//____________________________________________________________________________
PetrukhinShestakovModel::PetrukhinShestakovModel() :
MuELossI("genie::mueloss::PetrukhinShestakovModel")
{

}
//____________________________________________________________________________
PetrukhinShestakovModel::PetrukhinShestakovModel(string config) :
MuELossI("genie::mueloss::PetrukhinShestakovModel", config)
{

}
//____________________________________________________________________________
PetrukhinShestakovModel::~PetrukhinShestakovModel()
{

}
//____________________________________________________________________________
double PetrukhinShestakovModel::dE_dx(double E, MuELMaterial_t material) const
{
// Calculate the muon -dE/dx due to muon bremsstrahlung (in GeV^-2).
// To convert the result to more handly units, eg MeV/(gr/cm^2), just write:
// dE_dx /= (units::MeV/(units::g/units::cm2));

  if(material == eMuUndefined) return 0;
  if(E<=MuELProcess::Threshold(this->Process()) || E>=kMaxMuE) return 0;

  // material Z,E
  double Z = MuELMaterial::Z(material);
  double A = MuELMaterial::A(material);

  // calculate (the min,max) fraction of energy, v,  carried to the photon
  double Vmin = 0.;
  double Vmax = 1. - 0.75*kSqrte* (kMuonMass/E) * TMath::Power(Z,1/3.);

  // integrate the Bethe-Heitler muon bremsstrahlung cross section v*ds/dv
  // over v
  PetrukhinShestakovIntegrand vds_dv(E,Z);
  vds_dv.SetParam(0,"v",Vmin,Vmax);

  // calculate the b factor (bE = -dE/dx) in GeV^-3
  A *= units::g;
  double bbrem = (kNA/A) * fIntegrator->Integrate(vds_dv);

  // calculate the dE/dx due to muon bremsstrahlung in GeV^-2
  double de_dx = bbrem*E;
  return de_dx;
}
//____________________________________________________________________________
void PetrukhinShestakovModel::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PetrukhinShestakovModel::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PetrukhinShestakovModel::LoadConfig(void)
{
  fIntegrator = dynamic_cast<const IntegratorI *>(this->SubAlg("Integrator"));
  assert(fIntegrator);
}
//____________________________________________________________________________
PetrukhinShestakovIntegrand::PetrukhinShestakovIntegrand(double E, double Z) :
GSFunc(1)
{
  fE = E;
  fZ = Z;
}
//____________________________________________________________________________
PetrukhinShestakovIntegrand::~PetrukhinShestakovIntegrand()
{

}
//____________________________________________________________________________
double PetrukhinShestakovIntegrand::operator () (const vector<double> & x)
{
// Calculate the Bethe-Heitler cross section ds/dv for muon bremsstrahlung
// Returns v*(ds/dv)

  double v  = x[0]; // v, the fraction of energy transfered to the photon
  double v2 = TMath::Power(v,2.);

  if (! v >0) return 0;
  if (  v >1) return 0;
  if (! fE>0) return 0;

  // Some constants...
  double Z2    = TMath::Power(fZ,2.);
  double Zm13  = TMath::Power(fZ,-1./3.);
  double Zm23  = TMath::Power(fZ,-2./3.);
  double a3    = TMath::Power(kAem,3.); // (e/m coupling const)^3
  double me    = kElectronMass;
  double mmu   = kMuonMass;
  double mmu2  = kMuonMass2;
  double mmue  = mmu/me;
  double memu  = me/mmu;
  double memu2 = TMath::Power(memu,2);

  // Calculate the minimum monentum transfer to the nucleus (in GeV)
  double delta = (mmu2/fE) * 0.5*v/(1.-v);

  // Calculate the fi(delta) factor for the bremsstrahlung cross section
  // ds/dv according to the Petrukhin/Shestakov formula (dimensionless)
  double a  = ( (fZ<10) ? 189.*mmue * Zm13 : 189.*mmue * (2./3.)*Zm23 );
  double b  = 1. + (189./me) * Zm13 * kSqrte * delta;
  double fi = TMath::Log(a/b);

  // Calculate the Bethe-Heitler cross section ds/dv for muon bremsstrahlung
  double ds_dv  = (a3*memu2*kLe2) * (4*Z2) * (fi) * (4/3.-4*v/3.+v2)/v;
  double vds_dv = v*ds_dv;
  return vds_dv; // in GeV^-2
}
//____________________________________________________________________________

