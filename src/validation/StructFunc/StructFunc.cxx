//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Jelena Ilic <jelena.ilic \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ 25/04/12 - CA
   First added in v2.7.1

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TMatrixD.h>

#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "validation/StructFunc/StructFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::mc_vs_data;

//____________________________________________________________________________
StructFunc::StructFunc()
{
  fRESXSecModel      = 0;
  fDISXSecModel      = 0;
  fDISCharmXSecModel = 0;
}
//____________________________________________________________________________
bool StructFunc::ExtractF1F2xF3(
   double x, double Q2, int lepton_pdg, int nucleon_pdg,
   double & F1, double & F2, double & xF3) const
{
  // Extract F1, F2, xF3 from a "black box" cross-section model 

  // initialise structure functions
  F1  = 0;
  F2  = 0;
  xF3 = 0;

  const double sign = 1;
  const double M  = kNucleonMass;

  const int N = 3;

  double y[N]    = {0.25,  0.50,  0.75};
  double E[N]    = {0.,  0.,  0.};
  double xsec[N] = {0.,  0.,  0.};

  // E,y pairs for constant x,Q2 ( Q^2/2Mx = y*E)
  for(int i=0; i<N; i++) {
    E[i] = Q2/(2*M*x*y[i]);
  }
  // Calculate the total d2sigma/dxdy cross section for each {x,Q2,y[i],E[i]}
  for(int i=0; i<N; i++) {
    xsec[i] = this->d2SigTot_dxdy(x,y[i],E[i],lepton_pdg,nucleon_pdg);
  }
  //
  // Solve the system of equations:
  //   (Sigma) = (A) x (SF) => SF = (A)^-1 x (Sigma)
  //
  TMatrixD A     (3,3); /* (row,col) */
  TMatrixD Sigma (3,1);

  A(0,0) = x * y[0]*y[0]; //A(0,0) = x * TMath::Power(y[0],2.);  
  A(1,0) = x * y[1]*y[1]; //A(1,0) = x * TMath::Power(y[1],2.);  
  A(2,0) = x * y[2]*y[2]; //A(2,0) = x * TMath::Power(y[2],2.);  
  A(0,1) = 1. - y[0] - 0.5*M*x*y[0]/E[0];
  A(1,1) = 1. - y[1] - 0.5*M*x*y[1]/E[1];
  A(2,1) = 1. - y[2] - 0.5*M*x*y[2]/E[2];
  A(0,2) = sign * y[0] * (1.-0.5*y[0]);
  A(1,2) = sign * y[1] * (1.-0.5*y[1]);
  A(2,2) = sign * y[2] * (1.-0.5*y[2]);

  Sigma(0,0) = xsec[0] / (kGF2*M*E[0]/kPi);
  Sigma(1,0) = xsec[1] / (kGF2*M*E[1]/kPi);
  Sigma(2,0) = xsec[2] / (kGF2*M*E[2]/kPi);

  A.Invert();

  TMatrixD SF = A*Sigma;
  F1  = SF(0,0); 
  F2  = SF(1,0); 
  xF3 = SF(2,0); 

  return true;
}
//____________________________________________________________________________
double StructFunc::d2SigTot_dxdy(
    double x, double y, double E, int lepton_pdg, int nucleon_pdg) const
{
// Return inclusive d^2sigma/dxdy cross-section (only neutrino CC for the time-being)

  // If nucleon pdg is deuterium PDG then assume user asks cross-section for
  // an isoscalar target. Return avegare of free-proton and free-neutron
  // cross-sections
  if(nucleon_pdg == 1000010020) {
     double d2Sig_dxdy_p = this->d2SigTot_dxdy(x,y,E,lepton_pdg,kPdgProton);
     double d2Sig_dxdy_n = this->d2SigTot_dxdy(x,y,E,lepton_pdg,kPdgNeutron);
     return 0.5*(d2Sig_dxdy_p + d2Sig_dxdy_n);
  }

  // To proceed from this point on, ask that the input nucleon_pdg is either
  // a proton or a neutron
  assert(nucleon_pdg == kPdgProton || nucleon_pdg == kPdgNeutron);

  // Have free-nucleon target, not a nuclear target
  int tgt_pdg = 0;
  if(nucleon_pdg ==kPdgProton ) tgt_pdg = 1000010010;
  if(nucleon_pdg ==kPdgNeutron) tgt_pdg = 1000000010;

  // number of resonances and GENIE resonance IDs
  const int kNRes=18;  
  Resonance_t kResId[kNRes] = {
    kP33_1232, kS11_1535, kD13_1520, kS11_1650,
    kD13_1700, kD15_1675, kS31_1620, kD33_1700,
    kP11_1440, kP33_1600, kP13_1720, kF15_1680,
    kP31_1910, kP33_1920, kF35_1905, kF37_1950,
    kP11_1710, kF17_1970 
  };

  // Get W,Q2 from x,y
  double M  = kNucleonMass;
  double W  = 0;
  double Q2 = 0;
  utils::kinematics::XYtoWQ2 (E,M,W,Q2,x,y);

  //
  // Calculate resonance cross-section
  //
  double d2SigRES_dxdy  = 0;
  if(fRESXSecModel != 0) {
    Interaction * res_int = 
       Interaction::RESCC(tgt_pdg, nucleon_pdg, lepton_pdg, E);
    res_int->KinePtr()->SetW (W);
    res_int->KinePtr()->SetQ2(Q2);
    res_int->KinePtr()->Setx (x);
    res_int->KinePtr()->Sety (y);
    // loop over resonances
    for(int ires=0; ires<kNRes; ires++) {
       res_int->ExclTagPtr()->SetResonance(kResId[ires]);
       // resonance cross-section model typically computes d^2sigma/dWdQ2 but
       // Jacobean transformation must be implemented within the model to convert
       // cross-section to d^2sigma/dxdy.
       double xsec = fRESXSecModel->XSec(res_int,kPSxyfE); 
       xsec = TMath::Max(0., xsec);
       LOG("gvldtest", pINFO) 
         << "d2xsec_dxdy(CCRES; " << utils::res::AsString(kResId[ires])
	 << "; E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
	 << ", x = " << x << ", y = " << y << ") = "
	 << xsec/(1E-38*units::cm2) << " 1E-38 cm^2";
       d2SigRES_dxdy += xsec;
    }//res
    delete res_int;
  }

  //
  // Calculate DIS cross-section
  //
  double d2SigDIS_dxdy = 0;
  if(fDISXSecModel != 0) {
    // Note: Not setting quark ID. 
    // If the quark ID is set, the code returns (eg for neutrino CC) the vq->lq' cross-section.
    // But if a quark ID is not specified then the code loops over all relevant
    // valence and sea quark species and returns (eg for neutrino CC) the vN->lX cross-section.
    Interaction * dis_int = 
        Interaction::DISCC(tgt_pdg, nucleon_pdg, lepton_pdg, E);
    dis_int->KinePtr()->SetW (W);
    dis_int->KinePtr()->SetQ2(Q2);
    dis_int->KinePtr()->Setx (x);
    dis_int->KinePtr()->Sety (y);
    double xsec = fDISXSecModel->XSec(dis_int,kPSxyfE); 
    delete dis_int;
    xsec = TMath::Max(0., xsec);
    LOG("gvldtest", pINFO) 
       << "d2xsec_dxdy(CCDIS; " 
       << "E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
       << ", x = " << x << ", y = " << y << ") = "
       << xsec/(1E-38*units::cm2) << " 1E-38 cm^2";
    d2SigDIS_dxdy = xsec;
  }

  //
  // Calculate DIS charm cross-section
  //
  double d2SigDISc_dxdy = 0; 
  if(fDISCharmXSecModel != 0) {
    // Note: Not setting quark ID. See comments above.
    Interaction * dis_int = 
        Interaction::DISCC(tgt_pdg, nucleon_pdg, lepton_pdg, E);
    dis_int->ExclTagPtr()->SetCharm();
    dis_int->KinePtr()->SetW (W);
    dis_int->KinePtr()->SetQ2(Q2);
    dis_int->KinePtr()->Setx (x);
    dis_int->KinePtr()->Sety (y);
    double xsec = fDISCharmXSecModel->XSec(dis_int,kPSxyfE); 
    delete dis_int;
    xsec = TMath::Max(0., xsec);
    LOG("gvldtest", pINFO) 
       << "d2xsec_dxdy(CCDIS;charm; " 
       << "E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
       << ", x = " << x << ", y = " << y << ") = "
       << xsec/(1E-38*units::cm2) << " 1E-38 cm^2";
    d2SigDISc_dxdy = xsec;
  }

  // Return total d2sigma/dxdy
  double d2Sig_dxdy = d2SigDIS_dxdy + d2SigDISc_dxdy + d2SigRES_dxdy;
  return d2Sig_dxdy;
}
//____________________________________________________________________________
