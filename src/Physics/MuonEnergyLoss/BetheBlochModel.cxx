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

#include "Physics/MuonEnergyLoss/BetheBlochModel.h"
#include "Physics/MuonEnergyLoss/BetheBlochMaterialParams.h"
#include "Framework/Messenger/Messenger.h"

using namespace genie;
using namespace genie::mueloss;
using namespace genie::constants;

//____________________________________________________________________________
BetheBlochModel::BetheBlochModel() :
MuELossI("genie::mueloss::BetheBlochModel")
{

}
//____________________________________________________________________________
BetheBlochModel::BetheBlochModel(string config) :
MuELossI("genie::mueloss::BetheBlochModel", config)
{

}
//____________________________________________________________________________
BetheBlochModel::~BetheBlochModel()
{

}
//____________________________________________________________________________
double BetheBlochModel::dE_dx(double E, MuELMaterial_t mt) const
{
// Calculates ionization dE/dx for muons via Bethe-Bloch formula (in GeV^-2)
// To convert the result to more handly units, eg MeV/(gr/cm^2), just write:
// dE_dx /= (units::MeV/(units::g/units::cm2));

   if(mt == eMuUndefined) return 0;
   if(E<=MuELProcess::Threshold(this->Process()) || E>=kMaxMuE) return 0;

   double Z       = MuELMaterial::Z(mt);
   double A       = MuELMaterial::A(mt);
   double Z_A     = Z/A;              // in mol/gr
   double a2      = kAem2;            // (em coupling const)^2
   double Na      = kNA;              // Avogadro's number
   double lamda2 =  kLe2/units::cm2;  // (e compton wavelength)^2 in cm^2
   double me      = kElectronMass;    // in GeV
   double me2     = kElectronMass2; 
   double mmu     = kMuonMass;        // in GeV
   double mmu2    = kMuonMass2;
   double E2      = TMath::Power(E,2);
   double beta    = TMath::Sqrt(E2-mmu2)/E;
   double beta2   = TMath::Power(beta,2);
   double gamma   = E/mmu;
   double gamma2  = TMath::Power(gamma,2);
   double I       = BetheBlochMaterialParams::IonizationPotential(mt) * units::eV; 
   double I2      = TMath::Power(I,2); // in GeV^2

   // Calculate the maximum energy transfer to the electron (in GeV)

   double p2      = E2-mmu2;
   double Emaxt   = 2*me*p2 / (me2 + mmu2 + 2*me*E);
   double Emaxt2  = TMath::Power(Emaxt,2);

   // Calculate the density correction factor delta

   double X0 =  BetheBlochMaterialParams::DensityCorrection_X0(mt);
   double X1 =  BetheBlochMaterialParams::DensityCorrection_X1(mt);
   double a  =  BetheBlochMaterialParams::DensityCorrection_a(mt);
   double m  =  BetheBlochMaterialParams::DensityCorrection_m(mt);
   double C  =  BetheBlochMaterialParams::DensityCorrection_C(mt);
   double X  =  TMath::Log10(beta*gamma);

   double delta = 0;
   if(X0<X && X<X1) delta = 4.6052*X + a*TMath::Power(X1-X,m) + C;
   if(X>X1)         delta = 4.6052*X + C;

   LOG("MuELoss", pDEBUG) << "density correction factor = " << delta;
   LOG("MuELoss", pDEBUG) << "max energy transfer (GeV) = " << Emaxt;
   LOG("MuELoss", pDEBUG) << "ionization potential (GeV)= " << I;
   LOG("MuELoss", pDEBUG) << "E = " << E << ", p2 = " << p2;
   LOG("MuELoss", pDEBUG) << "beta = " << beta << ", gamma = " << gamma;

   // Calculate the -dE/dx
   double de_dx =  a2 * (2*kPi*Na*lamda2) * Z_A * (me/beta2) *
                    (TMath::Log( 2*me*beta2*gamma2*Emaxt/I2 ) - 
                                      2*beta2 + 0.25*(Emaxt2/E2) - delta);

   de_dx *= (units::GeV/(units::g/units::cm2));
   return de_dx; // in GeV^-2
}
//____________________________________________________________________________


