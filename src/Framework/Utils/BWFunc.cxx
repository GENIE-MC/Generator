//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - November 22, 2004

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 01, 2016 - Libo Jiang
   Add W dependence to Delta->N gamma 

*/
//____________________________________________________________________________

#include <cassert>

#include <TMath.h>

#include "Framework/Utils/BWFunc.h"
#include "Framework/Conventions/Constants.h"

using namespace genie;
using namespace genie::constants;

//
double genie::utils::bwfunc::BreitWignerLGamma(
               double W, int L, double mass, double width0, double norm)
{
//Inputs:
// - W:      Invariant mass (GeV)
// - L:      Resonance orbital angular momentum
// - mass:   Resonance mass (GeV)
// - width0: Resonance width
// - norm:   Breit Wigner norm

  //-- sanity checks
  assert(mass   >  0);
  assert(width0 >  0);
  assert(norm   >  0);
  assert(W      >  0);
  assert(L      >= 0);

  //-- auxiliary parameters
  double mN    = kNucleonMass;
  double mPi   = kPi0Mass;
  double m_2   = TMath::Power(mass, 2);
  double mN_2  = TMath::Power(mN,   2);
  double W_2   = TMath::Power(W,    2);
  
  double m=mass;
  //m_aux1 m_aux2
  double m_aux1= TMath::Power(mN+mPi, 2);
  double m_aux2= TMath::Power(mN-mPi, 2);

  double BRPi0    = 0.994; //Npi Branching Ratio
  double BRgamma0 = 0.006; //Ngamma Branching Ratio

  double widPi0   = width0*BRPi0;  
  double widgamma0= width0*BRgamma0;  

  //-- calculate the L-dependent resonance width
  double EgammaW= (W_2-mN_2)/(2*W);
  double Egammam= (m_2-mN_2)/(2*m);


  if(EgammaW<0) {
//  cout<< "Two small W!!! W is lower than one Nucleon Mass!!!!"<<endl;
  return 0;
  }
  //pPiW pion momentum 
  double pPiW     = 0;
  // 
  if(W_2>m_aux1) pPiW   = TMath::Sqrt((W_2-m_aux1)*(W_2-m_aux2))/(2*W);

  double pPim   = TMath::Sqrt((m_2-m_aux1)*(m_2-m_aux2))/(2*m);

  //double TPiW   = pPiW*TMath::Power(pPiW*rDelta, 2*L)/(1+TMath::Power(pPiW*rDelta, 2));
  //double TPim   = pPim*TMath::Power(pPim*rDelta, 2*L)/(1+TMath::Power(pPim*rDelta, 2));

  // Form factors
  //double fgammaW= 1/(TMath::Power(1+EgammaW*EgammaW/0.706, 2)*(1+EgammaW*EgammaW/3.519));
  //double fgammam= 1/(TMath::Power(1+Egammam*Egammam/0.706, 2)*(1+Egammam*Egammam/3.519));
  double fgammaW= 1/(TMath::Power(1+EgammaW*EgammaW/0.706, 2));
  double fgammam= 1/(TMath::Power(1+Egammam*Egammam/0.706, 2));


  double EgammaW_3=TMath::Power(EgammaW, 3);
  double Egammam_3=TMath::Power(Egammam, 3);
  double fgammaW_2=TMath::Power(fgammaW, 2);
  double fgammam_2=TMath::Power(fgammam, 2);

  //double width = widPi0*(TPiW/TPim)+widgamma0*(EgammaW_3*fgammaW_2/(Egammam_3*fgammam_2));
  double width = widPi0*TMath::Power((pPiW/pPim),3)+widgamma0*(EgammaW_3*fgammaW_2/(Egammam_3*fgammam_2));
  //-- calculate the Breit Wigner function for the input W
  double width_2 = TMath::Power( width,  2);
  double W_m_2   = TMath::Power( W-mass, 2);

  double bw = (0.5/kPi) * (width/norm) / (W_m_2 + 0.25*width_2);

  return bw;
}
//______________________________________________________________________
double genie::utils::bwfunc::BreitWignerL(
               double W, int L, double mass, double width0, double norm)
{
//Inputs:
// - W:      Invariant mass (GeV)
// - L:      Resonance orbital angular momentum
// - mass:   Resonance mass (GeV)
// - width0: Resonance width
// - norm:   Breit Wigner norm

  //-- sanity checks
  assert(mass   >  0);
  assert(width0 >  0);
  assert(norm   >  0);
  assert(W      >  0);
  assert(L      >= 0);

  //-- auxiliary parameters
  double mN    = kNucleonMass;
  double mPi   = kPi0Mass;
  double m_2   = TMath::Power(mass, 2);
  double mN_2  = TMath::Power(mN,   2);
  double mPi_2 = TMath::Power(mPi,  2);
  double W_2   = TMath::Power(W,    2);

  //-- calculate the L-dependent resonance width
  double qpW_2 = ( TMath::Power(W_2 - mN_2 - mPi_2, 2) - 4*mN_2*mPi_2 );
  double qpM_2 = ( TMath::Power(m_2 - mN_2 - mPi_2, 2) - 4*mN_2*mPi_2 );
  if(qpW_2 < 0) qpW_2 = 0;
  if(qpM_2 < 0) qpM_2 = 0;
  double qpW   = TMath::Sqrt(qpW_2) / (2*W);
  double qpM   = TMath::Sqrt(qpM_2) / (2*mass);
  double width = width0 * TMath::Power( qpW/qpM, 2*L+1 );

  //-- calculate the Breit Wigner function for the input W
  double width_2 = TMath::Power( width,  2);
  double W_m_2   = TMath::Power( W-mass, 2);

  double bw = (0.5/kPi) * (width/norm) / (W_m_2 + 0.25*width_2);
  return bw;
}
//______________________________________________________________________
double genie::utils::bwfunc::BreitWigner(
                       double W, double mass, double width, double norm)
{
//Inputs:
// - W:     Invariant mass (GeV)
// - mass:  Resonance mass (GeV)
// - width: Resonance width
// - norm:  Breit Wigner norm

  //-- sanity checks
  assert(mass  >  0);
  assert(width >  0);
  assert(norm  >  0);
  assert(W     >  0);

  //-- auxiliary parameters
  double width_2 = TMath::Power( width,  2);
  double W_m_2   = TMath::Power( W-mass, 2);

  //-- calculate the Breit Wigner function for the input W
  double bw = (0.5/kPi) * (width/norm) / (W_m_2 + 0.25*width_2);
  return bw;
}
//______________________________________________________________________

