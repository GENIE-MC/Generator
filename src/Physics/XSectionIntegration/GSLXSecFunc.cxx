//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>
#include <algorithm>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Numerical/GSLUtils.h"

using namespace genie;

//____________________________________________________________________________
genie::utils::gsl::dXSec_dQ2_E::dXSec_dQ2_E(
    const XSecAlgorithmI * m, const Interaction * i, double scale) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fScale(scale)
{

}
genie::utils::gsl::dXSec_dQ2_E::~dXSec_dQ2_E()
{

}
unsigned int genie::utils::gsl::dXSec_dQ2_E::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::dXSec_dQ2_E::DoEval(double xin) const
{
// inputs:  
//    Q2 [GeV^2]
// outputs: 
//   differential cross section [10^-38 cm^2 / GeV^2]
//
  double Q2 = xin;
  fInteraction->KinePtr()->SetQ2(Q2);
  double xsec = fModel->XSec(fInteraction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GSLXSecFunc", pDEBUG) << "xsec(Q2 = " << Q2 << ") = " << xsec;
#endif
  return fScale*xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::dXSec_dQ2_E::Clone() const
{
  return
    new genie::utils::gsl::dXSec_dQ2_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::dXSec_dy_E::dXSec_dy_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i)
{

}
genie::utils::gsl::dXSec_dy_E::~dXSec_dy_E()
{

}
unsigned int genie::utils::gsl::dXSec_dy_E::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::dXSec_dy_E::DoEval(double xin) const
{
// inputs:  
//    y [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double y = xin;
  fInteraction->KinePtr()->Sety(y);
  double xsec = fModel->XSec(fInteraction, kPSyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GXSecFunc", pDEBUG) << "xsec(y = " << y << ") = " << xsec;
#endif
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::dXSec_dy_E::Clone() const
{
  return
    new genie::utils::gsl::dXSec_dy_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dxdy_E::d2XSec_dxdy_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
  
}
genie::utils::gsl::d2XSec_dxdy_E::~d2XSec_dxdy_E()
{

}
unsigned int genie::utils::gsl::d2XSec_dxdy_E::NDim(void) const
{
  return 2;
}
double genie::utils::gsl::d2XSec_dxdy_E::DoEval(const double * xin) const
{
// inputs:  
//    x [-]
//    y [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double x = xin[0];
  double y = xin[1];
  fInteraction->KinePtr()->Setx(x);
  fInteraction->KinePtr()->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim * 
   genie::utils::gsl::d2XSec_dxdy_E::Clone() const
{
  return 
    new genie::utils::gsl::d2XSec_dxdy_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dQ2dy_E::d2XSec_dQ2dy_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
  
}
genie::utils::gsl::d2XSec_dQ2dy_E::~d2XSec_dQ2dy_E()
{

}
unsigned int genie::utils::gsl::d2XSec_dQ2dy_E::NDim(void) const
{
  return 2;
}
double genie::utils::gsl::d2XSec_dQ2dy_E::DoEval(const double * xin) const
{
// inputs:  
//   Q2 [-]
//    y [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double Q2 = xin[0];
  double  y = xin[1];
  fInteraction->KinePtr()->SetQ2(Q2);
  fInteraction->KinePtr()->Sety(y);
  kinematics::UpdateXFromQ2Y(fInteraction);
  double xsec = fModel->XSec(fInteraction, kPSQ2yfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim * 
   genie::utils::gsl::d2XSec_dQ2dy_E::Clone() const
{
  return 
    new genie::utils::gsl::d2XSec_dQ2dy_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dQ2dydt_E::d2XSec_dQ2dydt_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{
  
}
genie::utils::gsl::d2XSec_dQ2dydt_E::~d2XSec_dQ2dydt_E()
{

}
unsigned int genie::utils::gsl::d2XSec_dQ2dydt_E::NDim(void) const
{
  return 3;
}
double genie::utils::gsl::d2XSec_dQ2dydt_E::DoEval(const double * xin) const
{
// inputs:  
//   Q2 [-]
//    y [-]
//    t [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  //double  E = fInteraction->InitState().ProbeE(kRfLab);
  double Q2 = xin[0];
  double  y = xin[1];
  double  t = xin[2];
  fInteraction->KinePtr()->SetQ2(Q2);
  fInteraction->KinePtr()->Sety(y);
  fInteraction->KinePtr()->Sett(t);
  kinematics::UpdateXFromQ2Y(fInteraction);
  double xsec = fModel->XSec(fInteraction, kPSQ2yfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim * 
   genie::utils::gsl::d2XSec_dQ2dydt_E::Clone() const
{
  return 
    new genie::utils::gsl::d2XSec_dQ2dydt_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::d3XSec_dxdydt_E::d3XSec_dxdydt_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{

}
genie::utils::gsl::d3XSec_dxdydt_E::~d3XSec_dxdydt_E()
{

}
unsigned int genie::utils::gsl::d3XSec_dxdydt_E::NDim(void) const
{
  return 3;
}
double genie::utils::gsl::d3XSec_dxdydt_E::DoEval(const double * xin) const
{
// inputs:
//    x [-]
//    y [-]
//    t [-]
// outputs:
//   differential cross section [10^-38 cm^2]
//
  //double  E = fInteraction->InitState().ProbeE(kRfLab);
  double  x = xin[0];
  double  y = xin[1];
  double  t = xin[2];
  fInteraction->KinePtr()->Setx(x);
  fInteraction->KinePtr()->Sety(y);
  fInteraction->KinePtr()->Sett(t);
  double xsec = fModel->XSec(fInteraction, kPSxytfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d3XSec_dxdydt_E::Clone() const
{
  return
    new genie::utils::gsl::d3XSec_dxdydt_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dWdQ2_E::d2XSec_dWdQ2_E(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{

}
genie::utils::gsl::d2XSec_dWdQ2_E::~d2XSec_dWdQ2_E()
{

}
unsigned int genie::utils::gsl::d2XSec_dWdQ2_E::NDim(void) const
{
  return 2;
}
double genie::utils::gsl::d2XSec_dWdQ2_E::DoEval(const double * xin) const
{
// inputs:  
//    W  [GeV]
//    Q2 [GeV^2]
// outputs: 
//   differential cross section [10^-38 cm^2/GeV^3]
//
  double W  = xin[0];
  double Q2 = xin[1];
  fInteraction->KinePtr()->SetW(W);
  fInteraction->KinePtr()->SetQ2(Q2);
  if(fInteraction->ProcInfo().IsDeepInelastic() ||
     fInteraction->ProcInfo().IsDarkMatterDeepInelastic()) {
    double x=0,y=0;
    double E = fInteraction->InitState().ProbeE(kRfHitNucRest);
    double M = fInteraction->InitState().Tgt().HitNucP4Ptr()->M();

    kinematics::WQ2toXY(E,M,W,Q2,x,y);
    fInteraction->KinePtr()->Setx(x);
    fInteraction->KinePtr()->Sety(y);
  }
  double xsec = fModel->XSec(fInteraction, kPSWQ2fE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d2XSec_dWdQ2_E::Clone() const
{
  return 
    new genie::utils::gsl::d2XSec_dWdQ2_E(fModel,fInteraction);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dxdy_Ex::d2XSec_dxdy_Ex(
     const XSecAlgorithmI * m, const Interaction * i, double x) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fx(x)
{

}
genie::utils::gsl::d2XSec_dxdy_Ex::~d2XSec_dxdy_Ex()
{

}
unsigned int genie::utils::gsl::d2XSec_dxdy_Ex::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::d2XSec_dxdy_Ex::DoEval(double xin) const
{
// inputs:  
//    y [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double y = xin;
  fInteraction->KinePtr()->Setx(fx);
  fInteraction->KinePtr()->Sety(y);
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::d2XSec_dxdy_Ex::Clone() const
{
  return
    new genie::utils::gsl::d2XSec_dxdy_Ex(fModel,fInteraction,fx);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dxdy_Ey::d2XSec_dxdy_Ey(
     const XSecAlgorithmI * m, const Interaction * i, double y) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fy(y)
{

}
genie::utils::gsl::d2XSec_dxdy_Ey::~d2XSec_dxdy_Ey()
{

}
unsigned int genie::utils::gsl::d2XSec_dxdy_Ey::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::d2XSec_dxdy_Ey::DoEval(double xin) const
{
// inputs:  
//    x [-]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  double x = xin;
  fInteraction->KinePtr()->Setx(x);
  fInteraction->KinePtr()->Sety(fy);
  double xsec = fModel->XSec(fInteraction, kPSxyfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::d2XSec_dxdy_Ey::Clone() const
{
  return
    new genie::utils::gsl::d2XSec_dxdy_Ey(fModel,fInteraction,fy);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dWdQ2_EW::d2XSec_dWdQ2_EW(
     const XSecAlgorithmI * m, const Interaction * i, double W) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fW(W)
{

}
genie::utils::gsl::d2XSec_dWdQ2_EW::~d2XSec_dWdQ2_EW()
{

}
unsigned int genie::utils::gsl::d2XSec_dWdQ2_EW::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::d2XSec_dWdQ2_EW::DoEval(double xin) const
{
// inputs:  
//   Q2 [GeV^2]
// outputs: 
//   differential cross section [10^-38 cm^2/GeV^3]
//
  double Q2 = xin;
  fInteraction->KinePtr()->SetW(fW);
  fInteraction->KinePtr()->SetQ2(Q2);
  double xsec = fModel->XSec(fInteraction, kPSWQ2fE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::d2XSec_dWdQ2_EW::Clone(void) const
{
  return
    new genie::utils::gsl::d2XSec_dWdQ2_EW(fModel,fInteraction,fW);
}
//____________________________________________________________________________
genie::utils::gsl::d2XSec_dWdQ2_EQ2::d2XSec_dWdQ2_EQ2(
     const XSecAlgorithmI * m, const Interaction * i, double Q2) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
fQ2(Q2)
{

}
genie::utils::gsl::d2XSec_dWdQ2_EQ2::~d2XSec_dWdQ2_EQ2()
{

}
unsigned int genie::utils::gsl::d2XSec_dWdQ2_EQ2::NDim(void) const
{
  return 1;
}
double genie::utils::gsl::d2XSec_dWdQ2_EQ2::DoEval(double xin) const
{
// inputs:  
//   W [GeV]
// outputs: 
//   differential cross section [10^-38 cm^2/GeV^3]
//
  double W = xin;
  fInteraction->KinePtr()->SetW(W);
  fInteraction->KinePtr()->SetQ2(fQ2);
  double xsec = fModel->XSec(fInteraction,kPSWQ2fE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionOneDim *
   genie::utils::gsl::d2XSec_dWdQ2_EQ2::Clone() const
{
  return 
   new genie::utils::gsl::d2XSec_dWdQ2_EQ2(fModel,fInteraction,fQ2);
}
//____________________________________________________________________________
//
//  This just returns the 5-D differential cross section
//
genie::utils::gsl::d5XSecAR::d5XSecAR(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i),
flip(false)
{

}
genie::utils::gsl::d5XSecAR::~d5XSecAR()
{
}
 
unsigned int genie::utils::gsl::d5XSecAR::NDim(void) const
{
  return 5;
}
double genie::utils::gsl::d5XSecAR::DoEval(const double * xin) const
{
// inputs:
//    x [-]
//    y [-]
// outputs:
//   differential cross section [10^-38 cm^2]
//

  Kinematics * kinematics = fInteraction->KinePtr();
  const TLorentzVector * P4_nu = fInteraction->InitStatePtr()->GetProbeP4(kRfLab);
  double E_nu       = P4_nu->E();
  
  double E_l       = xin[0];
  double theta_l   = xin[1];
  double phi_l     = xin[2];
  double theta_pi  = xin[3];
  double phi_pi    = xin[4];
  
  double E_pi= E_nu-E_l;
  
  double y = E_pi/E_nu;
   
  double m_l = fInteraction->FSPrimLepton()->Mass();
  if (E_l < m_l) {
    return 0.;
  }
  
  double m_pi;
  if ( fInteraction->ProcInfo().IsWeakCC() ) {
    m_pi = constants::kPionMass;
  }
  else {
    m_pi = constants::kPi0Mass;
  }

   
  double p_l = TMath::Sqrt(E_l*E_l - m_l*m_l);
  TVector3 lepton_3vector = TVector3(0,0,0);
  lepton_3vector.SetMagThetaPhi(p_l,theta_l,phi_l);
  TLorentzVector P4_lep    = TLorentzVector(lepton_3vector , E_l );

  double p_pi = TMath::Sqrt(E_pi*E_pi - m_pi*m_pi);
  TVector3 pion_3vector = TVector3(0,0,0);
  pion_3vector.SetMagThetaPhi(p_pi,theta_pi,phi_pi);
  TLorentzVector P4_pion   = TLorentzVector(pion_3vector   , E_pi);
     
  double Q2 = -(*P4_nu-P4_lep).Mag2();

  double x = Q2/(2*E_pi*constants::kNucleonMass);

  Range1D_t xlim = fInteraction->PhaseSpace().XLim();

  if ( x <  xlim.min || x > xlim.max ) {
    return 0.;
  }
 
  kinematics->Setx(x);
  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);
  
  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_pion); // use Hadronic System variable to store pion momentum
 
  double xsec = fModel->XSec(fInteraction);
  if (xsec>0 && flip) {
    xsec = xsec*-1.0;
  }
  delete P4_nu;
  //return xsec/(1E-38 * units::cm2);
  return xsec;
}

ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d5XSecAR::Clone() const
{    
  return
    new genie::utils::gsl::d5XSecAR(fModel,fInteraction);
}

//____________________________________________________________________________
//
//  This is the original 5-D cross-section that Steve D coded
//
genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi::d5Xsec_dEldOmegaldOmegapi(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i)
{ 
  
}
genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi::~d5Xsec_dEldOmegaldOmegapi()
{  
  
}
unsigned int genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi::NDim(void) const
{
  return 5;
}
double genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi::DoEval(const double * xin) const
{  
// inputs:
//    x [-]
//    y [-]
// outputs:
//   differential cross section [10^-38 cm^2]
//
  Kinematics * kinematics = fInteraction->KinePtr();
  const TLorentzVector * P4_nu = fInteraction->InitStatePtr()->GetProbeP4(kRfLab);
  double E_nu       = P4_nu->E();

  double E_l       = xin[0];
  double theta_l   = xin[1];
  double phi_l     = xin[2];
  double theta_pi  = xin[3];
  double phi_pi    = xin[4];

  double E_pi= E_nu-E_l;

  double y = E_pi/E_nu;
    
  double m_l = fInteraction->FSPrimLepton()->Mass();
  if (E_l < m_l) {
    return 0.;
  }
 
  double m_pi;
  if ( fInteraction->ProcInfo().IsWeakCC() ) {
    m_pi = constants::kPionMass;
  }
  else {
    m_pi = constants::kPi0Mass;
  }

  double p_l = TMath::Sqrt(E_l*E_l - m_l*m_l);
  TVector3 lepton_3vector = TVector3(0,0,0);
  lepton_3vector.SetMagThetaPhi(p_l,theta_l,phi_l);
  TLorentzVector P4_lep    = TLorentzVector(lepton_3vector , E_l );
  
  double p_pi = TMath::Sqrt(E_pi*E_pi - m_pi*m_pi);
  TVector3 pion_3vector = TVector3(0,0,0);
  pion_3vector.SetMagThetaPhi(p_pi,theta_pi,phi_pi);
  TLorentzVector P4_pion   = TLorentzVector(pion_3vector   , E_pi);

  double Q2 = -(*P4_nu-P4_lep).Mag2();

  double x = Q2/(2*E_pi*constants::kNucleonMass);

  Range1D_t xlim = fInteraction->PhaseSpace().XLim();
    
  if ( x <  xlim.min || x > xlim.max ) {
    return 0.;
  }
 
  kinematics->Setx(x);
  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);
  
  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_pion); // use Hadronic System variable to store pion momentum
 
  delete P4_nu;

  double xsec = fModel->XSec(fInteraction)*TMath::Sin(theta_l)*TMath::Sin(theta_pi);
  return xsec/(1E-38 * units::cm2);
}

ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi::Clone() const
{
  return
    new genie::utils::gsl::d5Xsec_dEldOmegaldOmegapi(fModel,fInteraction);
}

//____________________________________________________________________________
//
// This is the same as the 5d space that Steve D coded,   
// But we remove the integration of phi_l
//

genie::utils::gsl::d4Xsec_dEldThetaldOmegapi::d4Xsec_dEldThetaldOmegapi(
									const XSecAlgorithmI * m, 
									const Interaction * i) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fModel(m),
  fInteraction(i)
{
  
}
genie::utils::gsl::d4Xsec_dEldThetaldOmegapi::~d4Xsec_dEldThetaldOmegapi()
{

}
unsigned int genie::utils::gsl::d4Xsec_dEldThetaldOmegapi::NDim(void) const
{
  return 4;
}
double genie::utils::gsl::d4Xsec_dEldThetaldOmegapi::DoEval(const double * xin) const
{
// inputs:  
//    El [GeV]
//    theta l [rad]
//    theta pi [rad]
//    phi pi [rad]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  Kinematics * kinematics = fInteraction->KinePtr();
  const TLorentzVector * P4_nu = fInteraction->InitStatePtr()->GetProbeP4(kRfLab);
  double E_nu       = P4_nu->E();
  
  double E_l       = xin[0];
  double theta_l   = xin[1];
  double phi_l     = 0.0;
  double theta_pi  = xin[2];
  double phi_pi    = xin[3];
  
  double sin_theta_l  = TMath::Sin(theta_l);
  double sin_theta_pi = TMath::Sin(theta_pi);
  
  double E_pi= E_nu-E_l;
    
  double y = E_pi/E_nu;
  
  double m_l = fInteraction->FSPrimLepton()->Mass();
  if (E_l < m_l) {
    return 0.;
  }
  
  double m_pi;
  if ( fInteraction->ProcInfo().IsWeakCC() ) {
    m_pi = constants::kPionMass;
  }
  else {
    m_pi = constants::kPi0Mass;
  }
  
  double p_l = TMath::Sqrt(E_l*E_l - m_l*m_l);
  TVector3 lepton_3vector = TVector3(0,0,0);
  lepton_3vector.SetMagThetaPhi(p_l,theta_l,phi_l);
  TLorentzVector P4_lep    = TLorentzVector(lepton_3vector , E_l );
  
  double p_pi = TMath::Sqrt(E_pi*E_pi - m_pi*m_pi);
  TVector3 pion_3vector = TVector3(0,0,0);
  pion_3vector.SetMagThetaPhi(p_pi,theta_pi,phi_pi);
  TLorentzVector P4_pion   = TLorentzVector(pion_3vector   , E_pi);
  
  double Q2 = -(*P4_nu-P4_lep).Mag2();
  
  double x = Q2/(2*E_pi*constants::kNucleonMass);
  
  Range1D_t xlim = fInteraction->PhaseSpace().XLim();
  
  if ( x <  xlim.min || x > xlim.max ) {
    return 0.;
  }

  kinematics->Setx(x);
  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);
  
  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_pion); // use Hadronic System variable to store pion momentum
  
  delete P4_nu;
  
  double xsec = sin_theta_l * sin_theta_pi * fModel->XSec(fInteraction,kPSElOlTpifE);
  return fFactor * xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim * 
   genie::utils::gsl::d4Xsec_dEldThetaldOmegapi::Clone() const
{
  return 
    new genie::utils::gsl::d4Xsec_dEldThetaldOmegapi(fModel,fInteraction);
}
void genie::utils::gsl::d4Xsec_dEldThetaldOmegapi::SetFactor(double factor)
{
  fFactor = factor;
}
double genie::utils::gsl::d4Xsec_dEldThetaldOmegapi::GetFactor(void) const
{
  return fFactor;
}
//____________________________________________________________________________
genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig::d4Xsec_dEgdThetaldThetagdPhig(
										const XSecAlgorithmI * m, 
										const Interaction * i) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fModel(m),
  fInteraction(i),
  fFactor(1.) 
{
  
}
genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig::~d4Xsec_dEgdThetaldThetagdPhig() 
{

}
unsigned int genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig::NDim(void) const
{
  return 4;
}
//____________________________________________________________________________
double genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig::DoEval(const double * xin) const
{
  // inputs:  
  //    E gamma [GeV]
  //    theta l [rad]
  //    theta gamma [rad]
  //    phi gamma [rad]
  //
  // outputs: 
  //   differential cross section [10^-38 cm^2]
  //
 
  Kinematics * kinematics = fInteraction->KinePtr();
  double E_nu = fInteraction->InitState().ProbeE(kRfLab);
  const TLorentzVector P4_nu( 0., 0., E_nu, E_nu ) ;
  
  //double E_g       = xin[0];  // No need to define this, left here for comments

  // Complex correct formula for E_l that takes into account the target mass
  double m_t = fInteraction->InitState().Tgt().Mass() ;
  double cos_theta_g = cos( xin[2] ) ;
  double cos_theta_l = cos( xin[1] ) ; 
  double E_l = ( m_t*(E_nu - xin[0]) - E_nu*xin[0]*(1. - cos_theta_g) ) 
    / ( m_t + E_nu*(1.-cos_theta_l) 
	- xin[0]*(1.-cos_theta_l*cos_theta_g  - sin(xin[1])*sin(xin[2])*cos(xin[3]) ) ) ;

  //double theta_l   = xin[1];
  //double phi_l     = 0.0;
  
  //double theta_g   = xin[2];
  //double phi_g     = xin[3];
  
  TVector3 lepton_3vector = TVector3(0,0,0);
  lepton_3vector.SetMagThetaPhi( E_l, xin[1], 0. ) ;
  TLorentzVector P4_lep = TLorentzVector( lepton_3vector , E_l );
  
  TVector3 photon_3vector = TVector3(0,0,0);
  photon_3vector.SetMagThetaPhi( xin[0], xin[2], xin[3] ) ;
  TLorentzVector P4_photon   = TLorentzVector(photon_3vector, xin[0] );

  TLorentzVector q = P4_nu-P4_lep ;
 
  double Q2 = -q.Mag2();

  double x = Q2/(2 * xin[0] * constants::kNucleonMass );
  
  // Range1D_t xlim = fInteraction->PhaseSpace().XLim();
  
  // if ( x <  xlim.min || x > xlim.max ) {
  //   return 0.;
  // }
  
  kinematics->Setx(x);

  double y = xin[0]/E_nu;

  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);

  double t = TMath::Abs( (q - P4_photon).Mag2() );
  kinematics -> Sett( t ) ;

  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_photon); // use Hadronic System variable to store photon momentum
  
  double xsec = fModel->XSec(fInteraction,kPSEgTlTgPgfE);
  
  return fFactor * xsec/(1E-38 * units::cm2);
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim * 
genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig::Clone() const
{
  return 
    new genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig(fModel,fInteraction);
}
void genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig::SetFactor(double factor)
{
  fFactor = factor;
}
double genie::utils::gsl::d4Xsec_dEgdThetaldThetagdPhig::GetFactor(void) const
{
  return fFactor;
}
//____________________________________________________________________________
genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig::d4Xsec_dEgdtdThetagdPhig(
										const XSecAlgorithmI * m, 
										const Interaction * i) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fModel(m),
  fInteraction(i),
  fFactor(1.) 
{
  
}
genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig::~d4Xsec_dEgdtdThetagdPhig() 
{

}
unsigned int genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig::NDim(void) const
{
  return 4;
}
//____________________________________________________________________________
double genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig::DoEval(const double * xin) const
{
  // inputs:  
  //    E gamma [GeV]     = xin[0]
  //    t [GeV^2]         = xin[1]
  //    theta gamma [rad] = xin[2]
  //    phi gamma [rad]   = xin[3]
  //
  // outputs: 
  //   differential cross section [10^-38 cm^2]
  //
 
  Kinematics * kinematics = fInteraction->KinePtr();
  double E_nu = fInteraction->InitState().ProbeE(kRfLab);
  const TLorentzVector P4_nu( 0., 0., E_nu, E_nu ) ;
  
  double m_t = fInteraction->InitState().Tgt().Mass() ;

  // Compute values so there are no repeated computations
  double cos_theta_g = cos( xin[2] ) ;
  double A = xin[0] * cos_theta_g - E_nu  ;
  double B = E_nu - xin[0] ;
  double C = xin[0] * sin( xin[2] ) * cos( xin[3] ) ;

  // E_l and theta_l as a function of t
  double E_l = B - (0.5*xin[1])/m_t ;
  double alpha = atan2( C, A ) ;

  double arcosin = acos( ( (0.5*xin[1] - E_nu * xin[0]*( 1-cos_theta_g ) )/E_l - B  ) / 
			 sqrt( pow( A, 2 ) + pow( C, 2 ) ) 
			 ) ;

  std::array<double,2> solutions = { alpha + arcosin, alpha - arcosin } ;
  for ( auto & s : solutions ) {
    if ( s < 0. ) s+= 2 * constants::kPi ;
  }
  
  double theta_l = std::min( solutions[0], solutions[1] ) ;

  TVector3 lepton_3vector = TVector3(0,0,0);
  lepton_3vector.SetMagThetaPhi( E_l, theta_l, 0. ) ;
  TLorentzVector P4_lep = TLorentzVector( lepton_3vector , E_l );
  
  TVector3 photon_3vector = TVector3(0,0,0);
  photon_3vector.SetMagThetaPhi( xin[0], xin[2], xin[3] ) ;
  TLorentzVector P4_photon = TLorentzVector(photon_3vector, xin[0] );

  TLorentzVector q = P4_nu-P4_lep ;
 
  double Q2 = -q.Mag2();

  double x = Q2/( 2 * xin[0] * constants::kNucleonMass );
  
  // Range1D_t xlim = fInteraction->PhaseSpace().XLim();
  
  // if ( x <  xlim.min || x > xlim.max ) {
  //   return 0.;
  // }
  
  kinematics->Setx(x);

  double y = xin[0]/E_nu;

  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);

  kinematics -> Sett( xin[1] ) ;

  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_photon); // use Hadronic System variable to store photon momentum
  
  double xsec = fModel->XSec(fInteraction,kPSEgtTgPgfE);
  
  return fFactor * xsec/(1E-38 * units::cm2);
}

//____________________________________________________________________________
ROOT::Math::IBaseFunctionMultiDim * 
genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig::Clone() const
{
  return 
    new genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig(fModel,fInteraction);
}
void genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig::SetFactor(double factor)
{
  fFactor = factor;
}
double genie::utils::gsl::d4Xsec_dEgdtdThetagdPhig::GetFactor(void) const
{
  return fFactor;
}
//____________________________________________________________________________
genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag::d5Xsec_dEgdOmegaldOmegag(
										const XSecAlgorithmI * m, 
										const Interaction * i) :
  ROOT::Math::IBaseFunctionMultiDim(),
  fModel(m),
  fInteraction(i),
  fFactor(1.) 
{
  
}
genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag::~d5Xsec_dEgdOmegaldOmegag() 
{

}
unsigned int genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag::NDim(void) const
{
  return 4;
}
double genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag::DoEval(const double * xin) const
{
  // inputs:  
  //    E gamma [GeV]
  //    cos theta l 
  //    cos theta gamma
  //    phi gamma [rad]
  //
  // outputs: 
  //   differential cross section [10^-38 cm^2]
  //
 
  Kinematics * kinematics = fInteraction->KinePtr();
  double E_nu = fInteraction->InitState().ProbeE(kRfLab);
  const TLorentzVector P4_nu( 0., 0., E_nu, E_nu ) ;
  
  double sin_theta_l = sqrt( 1. - xin[1]*xin[1] ) ;
  double sin_theta_g = sqrt( 1. - xin[2]*xin[2] ) ;
  
  //double E_g       = xin[0];  // No need to define this, left here for comments
  
  // Complex correct formula for E_l that takes into account the target mass
  double m_t = fInteraction->InitState().Tgt().Mass() ;
  //double cos_theta_g = xin[2] ;
  //double cos_theta_l = xin[1]  ; 
  double E_l = ( m_t*(E_nu - xin[0]) - E_nu*xin[0]*(1. - xin[2]) ) 
    / ( m_t + E_nu*(1.-xin[1]) 
	- xin[0]*(1.-xin[1]*xin[2] - sin_theta_l*sin_theta_g*cos(xin[3]) ) ) ;

  //double theta_l   = xin[1];
  //double phi_l     = 0.0;
  
  //double theta_g   = xin[2];
  //double phi_g     = xin[3];
  
  TVector3 lepton_3vector( E_l * sin_theta_l, 0., E_l * xin[1] );
  TLorentzVector P4_lep = TLorentzVector( lepton_3vector , E_l );
  
  TVector3 photon_3vector(xin[0] * sin_theta_g * cos(xin[3] ), 
			  xin[0] * sin_theta_g * sin(xin[3] ), 
			  xin[0] * xin[2] ) ;  
  TLorentzVector P4_photon   = TLorentzVector(photon_3vector, xin[0] );

  TLorentzVector q = P4_nu-P4_lep ;
  
  double Q2 = -q.Mag2();
  
  double x = Q2/(2 * xin[0] * constants::kNucleonMass );
  
  // Range1D_t xlim = fInteraction->PhaseSpace().XLim();
  
  // if ( x <  xlim.min || x > xlim.max ) {
  //   return 0.;
  // }
  
  kinematics->Setx(x);

  double y = xin[0]/E_nu;

  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);
  
  double t = TMath::Abs( (q - P4_photon).Mag2() );
  kinematics -> Sett( t ) ;
  
  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_photon); // use Hadronic System variable to store photon momentum
  
  double xsec = fModel->XSec(fInteraction,kPSEgOlOgfE);
  
  return fFactor * xsec/(1E-38 * units::cm2);
}

ROOT::Math::IBaseFunctionMultiDim * 
genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag::Clone() const
{
  return 
    new genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag(fModel,fInteraction);
}
void genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag::SetFactor(double factor)
{
  fFactor = factor;
}
double genie::utils::gsl::d5Xsec_dEgdOmegaldOmegag::GetFactor(void) const
{
  return fFactor;
}
//____________________________________________________________________________
genie::utils::gsl::d3Xsec_dOmegaldThetapi::d3Xsec_dOmegaldThetapi(
     const XSecAlgorithmI * m, const Interaction * i) :
ROOT::Math::IBaseFunctionMultiDim(),
fModel(m),
fInteraction(i),
fElep(-1)
{
  
}
genie::utils::gsl::d3Xsec_dOmegaldThetapi::~d3Xsec_dOmegaldThetapi()
{

}
unsigned int genie::utils::gsl::d3Xsec_dOmegaldThetapi::NDim(void) const
{
  return 3;
}
double genie::utils::gsl::d3Xsec_dOmegaldThetapi::DoEval(const double * xin) const
{
// inputs:  
//    theta l [rad]
//    phi_l   [rad]
//    phi pi  [rad]
// outputs: 
//   differential cross section [10^-38 cm^2]
//
  Kinematics * kinematics = fInteraction->KinePtr();
  const TLorentzVector * P4_nu = fInteraction->InitStatePtr()->GetProbeP4(kRfLab);
  double E_nu = P4_nu->E();
  
  double E_l = fElep;
  
  double theta_l   = xin[0];
  double phi_l     = xin[1];
  double theta_pi  = xin[2];
  double phi_pi    = 0;
  
  double sin_theta_l  = TMath::Sin(theta_l);
  double sin_theta_pi = TMath::Sin(theta_pi);
  
  double E_pi= E_nu-E_l;
    
  double y = E_pi/E_nu;
  
  double m_l = fInteraction->FSPrimLepton()->Mass();
  if (E_l < m_l) {
    return 0.;
  }
  
  double m_pi;
  if ( fInteraction->ProcInfo().IsWeakCC() ) {
    m_pi = constants::kPionMass;
  }
  else {
    m_pi = constants::kPi0Mass;
  }
  
  double p_l = TMath::Sqrt(E_l*E_l - m_l*m_l);
  TVector3 lepton_3vector = TVector3(0,0,0);
  lepton_3vector.SetMagThetaPhi(p_l,theta_l,phi_l);
  TLorentzVector P4_lep    = TLorentzVector(lepton_3vector , E_l );
  
  double p_pi = TMath::Sqrt(E_pi*E_pi - m_pi*m_pi);
  TVector3 pion_3vector = TVector3(0,0,0);
  pion_3vector.SetMagThetaPhi(p_pi,theta_pi,phi_pi);
  TLorentzVector P4_pion   = TLorentzVector(pion_3vector   , E_pi);
  
  double Q2 = -(*P4_nu-P4_lep).Mag2();
  
  double x = Q2/(2*E_pi*constants::kNucleonMass);
  
  Range1D_t xlim = fInteraction->PhaseSpace().XLim();
  
  if ( x <  xlim.min || x > xlim.max ) {
    return 0.;
  }
  
  kinematics->Setx(x);
  kinematics->Sety(y);
  kinematics::UpdateWQ2FromXY(fInteraction);
  
  kinematics->SetFSLeptonP4(P4_lep );
  kinematics->SetHadSystP4 (P4_pion); // use Hadronic System variable to store pion momentum
  
  delete P4_nu;
  
  double xsec = (sin_theta_l * sin_theta_pi) * fModel->XSec(fInteraction,kPSElOlTpifE);
  return xsec/(1E-38 * units::cm2);
}
genie::utils::gsl::d3Xsec_dOmegaldThetapi * 
   genie::utils::gsl::d3Xsec_dOmegaldThetapi::Clone() const
{
  d3Xsec_dOmegaldThetapi * out = new genie::utils::gsl::d3Xsec_dOmegaldThetapi(fModel,fInteraction);
  out->SetE_lep(fElep);
  return out;
}
//____________________________________________________________________________
void genie::utils::gsl::d3Xsec_dOmegaldThetapi::SetE_lep(double E_lepton) const
{
  fElep = E_lepton;
}
//____________________________________________________________________________
genie::utils::gsl::dXSec_dElep_AR_pion::dXSec_dElep_AR_pion(
    const XSecAlgorithmI * m, const Interaction * i,
    string gsl_nd_integrator_type, double gsl_relative_tolerance,
    unsigned int max_n_calls) :
ROOT::Math::IBaseFunctionOneDim(),
fModel(m),
fInteraction(i),
integrator(utils::gsl::IntegrationNDimTypeFromString(gsl_nd_integrator_type),1, gsl_relative_tolerance, max_n_calls),
fGSLIntegratorType(gsl_nd_integrator_type),
fGSLRelTol(gsl_relative_tolerance),
fGSLMaxCalls(max_n_calls)
{
  func = new utils::gsl::d3Xsec_dOmegaldThetapi(fModel, fInteraction);
  
  integrator.SetRelTolerance(gsl_relative_tolerance);
  integrator.SetFunction(*func);
  
  kine_min[0] = kine_min[1] = kine_min[2] = controls::kASmallNum;
  
  kine_max[0] = kine_max[2] = constants::kPi-controls::kASmallNum;
  kine_max[1] = 2 * constants::kPi-controls::kASmallNum;
}
genie::utils::gsl::dXSec_dElep_AR_pion::~dXSec_dElep_AR_pion()
{
  delete func;
}
double genie::utils::gsl::dXSec_dElep_AR_pion::DoEval(double xin) const
{
  double Elep = xin;
  func->SetE_lep(Elep);
  double xsec = integrator.Integral(&kine_min[0], &kine_max[0]) ;
  LOG("GSLXSecFunc",pINFO) << "dXSec_dElep_AR_pion >> "<<func->NDim()<<"d integral done. (Elep = " <<Elep<< " , dxsec/dElep = "<<xsec << ")";
  return xsec;
}
genie::utils::gsl::dXSec_dElep_AR_pion *
   genie::utils::gsl::dXSec_dElep_AR_pion::Clone() const
{
  return
    new genie::utils::gsl::dXSec_dElep_AR_pion(fModel,fInteraction, fGSLIntegratorType, fGSLRelTol, fGSLMaxCalls);
}
//____________________________________________________________________________
genie::utils::gsl::dXSec_Log_Wrapper::dXSec_Log_Wrapper(
      const ROOT::Math::IBaseFunctionMultiDim * fn,
      bool * ifLog, double * mins, double * maxes) :
  fFn(fn),
  fIfLog(ifLog),
  fMins(mins),
  fMaxes(maxes)
{
}
genie::utils::gsl::dXSec_Log_Wrapper::~dXSec_Log_Wrapper()
{
} 
 
// ROOT::Math::IBaseFunctionMultiDim interface
unsigned int genie::utils::gsl::dXSec_Log_Wrapper::NDim   (void) const  
{
  return fFn->NDim();
}
double genie::utils::gsl::dXSec_Log_Wrapper::DoEval (const double * xin) const
{
  double * toEval = new double[this->NDim()];
  double a,b,x;
  for (unsigned int i = 0 ; i < this->NDim() ; i++ )
  {
    if (fIfLog[i]) {
      a = fMins[i];
      b = fMaxes[i];
      x = xin[i];
      toEval[i] = a + (b-a)/(constants::kNapierConst-1.) * (exp(x/(b-a)) - 1.);
    }
    else {
      toEval[i] = xin[i];
    }
  }
  double val = (*fFn)(toEval);
  delete[] toEval;
  return val; 
}
ROOT::Math::IBaseFunctionMultiDim * genie::utils::gsl::dXSec_Log_Wrapper::Clone (void) const
{
  return new dXSec_Log_Wrapper(fFn,fIfLog,fMins,fMaxes);
}
    
//____________________________________________________________________________

