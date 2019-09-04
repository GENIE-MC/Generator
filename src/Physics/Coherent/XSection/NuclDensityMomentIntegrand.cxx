//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at cern.ch>
         University of Liverpool & STFC Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/Coherent/XSection/NuclDensityMomentIntegrand.h"

using namespace genie;

//____________________________________________________________________________
utils::gsl::wrap::NuclDensityMomentIntegrand::NuclDensityMomentIntegrand(
 int A, int k):
ROOT::Math::IBaseFunctionOneDim()
{
  fA  = A;
  fK  = k;

  if(fK < 0) {
     LOG("Nuclear", pDEBUG)
        << "Sure you want to calculate an inverse nuclear density moment ("
        << "E[r^{" << fK << "}]) ?";
  }
  if(fA <= 1) {
     LOG("Nuclear", pWARN)
        << "The atomic mass number A should be >1 (input value was: "
        << fA << ")";
  }
}
//____________________________________________________________________________
utils::gsl::wrap::NuclDensityMomentIntegrand::~NuclDensityMomentIntegrand()
{

}
//____________________________________________________________________________
unsigned int utils::gsl::wrap::NuclDensityMomentIntegrand::NDim(void) const
{
  return 1;
}
//____________________________________________________________________________
double utils::gsl::wrap::NuclDensityMomentIntegrand::DoEval(double r) const
{
  if(fA <= 1) return 0;

  double rho  = utils::nuclear::Density(r,fA);
  double rkth = TMath::Power(r,fK);

  double integrand = rho * rkth;
  return integrand;
}
//____________________________________________________________________________
ROOT::Math::IBaseFunctionOneDim *
  utils::gsl::wrap::NuclDensityMomentIntegrand::Clone(void) const
{
  return new utils::gsl::wrap::NuclDensityMomentIntegrand(fA, fK);
}
//____________________________________________________________________________
