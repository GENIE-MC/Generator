//____________________________________________________________________________
/*
  Copyright (c) 2003-2020, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Author: Iker de Icaza <i.de-icaza-astiz \at sussex.ac.uk>
  University of Sussex

  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
  University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Physics/DarkNeutrino/XSection/EngelFormFactor.h"
#include "Framework/Messenger/Messenger.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"

#include "Framework/Utils/StringUtils.h"

using namespace genie;


//____________________________________________________________________________
EngelFormFactor::EngelFormFactor() :
Algorithm("genie::EngelFormFactor")
{

}
//____________________________________________________________________________
EngelFormFactor::EngelFormFactor(string config) :
Algorithm("genie::EngelFormFactor", config)
{

}
//____________________________________________________________________________
EngelFormFactor::~EngelFormFactor()
{

}
//____________________________________________________________________________
double EngelFormFactor::FormFactor(const double Q2, const Target & target) const {

  if(!target.IsValidNucleus()) return 0.;

  // const double kappa = TMath::Sqrt(2.*M*Er); // Q^2 \sim -2\kappa^2
  // const double kappa2 = kappa * kappa;

  const double M = target.Mass(); // units: GeV
  const double A = target.A();

  // const double Er = (Q*Q)/(-2.*M);
  // const double kappa2 = -Q * Q;
  // const double kappa = TMath::Abs(kappa2);
  // const double Q2 = Q * Q;
  const double s = 1.*units::fm;
  const double s2 = s * s;
  const double R = 1.2*TMath::Power(A, 1./3.)*units::fm;
  const double r = TMath::Sqrt(R*R - 5.*s*s);
  const double qr = TMath::Sqrt(Q2) * r;

  const double f1 = 3. * TMath::Exp(-.5*Q2*s2) * TMath::Power(qr, -3.);
  const double f2 = TMath::Sin(qr) - qr*TMath::Cos(qr);
  return f1 * f2;
}
//____________________________________________________________________________
void EngelFormFactor::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void EngelFormFactor::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void EngelFormFactor::LoadConfig(void)
{

  // // load R33 parameters
  // this -> GetParamVect( "DV-Coefficient", fFBCs ) ;

  // GetParam( "DV-Radius", fRadius ) ;
  // fRadius *= units::fm ;

  // GetParam( "DV-Nucleus", fPDG ) ;

  // LOG("EngelFormFactor", pINFO) << "Loaded " << fFBCs.size() << " coeffictients for nucleus " << fPDG ; 

}
//____________________________________________________________________________
