//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - June 15, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Physics/Hadronization/CollinsSpillerFragm.h"
#include "Physics/Hadronization/FragmentationFunctions.h"

using namespace genie;

//___________________________________________________________________________
CollinsSpillerFragm::CollinsSpillerFragm() :
FragmentationFunctionI("genie::CollinsSpillerFragm")
{

}
//___________________________________________________________________________
CollinsSpillerFragm::CollinsSpillerFragm(string config) :
FragmentationFunctionI("genie::CollinsSpillerFragm", config)
{

}
//___________________________________________________________________________
CollinsSpillerFragm::~CollinsSpillerFragm()
{
  delete fFunc;
}
//___________________________________________________________________________
double CollinsSpillerFragm::Value(double z) const
{
// Evaluate the fragmentation function

  if(z<0 || z>1) return 0;
  return fFunc->Eval(z);
}
//___________________________________________________________________________
double CollinsSpillerFragm::GenerateZ(void) const
{
// Return a random number using the fragmentation function as PDF

  return fFunc->GetRandom();
}
//___________________________________________________________________________
void CollinsSpillerFragm::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->BuildFunction();
}
//___________________________________________________________________________
void CollinsSpillerFragm::Configure(string config)
{
  Algorithm::Configure(config);
  this->BuildFunction();
}
//___________________________________________________________________________
void CollinsSpillerFragm::BuildFunction(void)
{
  fFunc = new TF1("fFunc",genie::utils::frgmfunc::collins_spiller_func,0,1,2);

  fFunc->SetParNames("Norm","Epsilon");

  double N = -1. ;
  GetParam( "CSFrag-Norm", N, false ) ;

  double e = 0. ;
  GetParam( "CSFrag-Epsilon", e, false ) ;

  // if the normalization parameter was left negative, explicitly normalize
  // the fragmentation function
  if(N<0) {
    N=1;
    fFunc->SetParameters(N,e);
    double I = fFunc->Integral(0,1);
    assert(I>0);
    N = 1./I;
  } 
  fFunc->SetParameters(N,e);
}
//___________________________________________________________________________



