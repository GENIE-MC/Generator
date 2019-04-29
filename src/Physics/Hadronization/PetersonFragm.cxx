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

#include <TROOT.h>

#include "Physics/Hadronization/PetersonFragm.h"
#include "Physics/Hadronization/FragmentationFunctions.h"

using namespace genie;

//___________________________________________________________________________
PetersonFragm::PetersonFragm() :
FragmentationFunctionI("genie::PetersonFragm")
{

}
//___________________________________________________________________________
PetersonFragm::PetersonFragm(string config) :
FragmentationFunctionI("genie::PetersonFragm", config)
{
  this->BuildFunction();
}
//___________________________________________________________________________
PetersonFragm::~PetersonFragm()
{
  delete fFunc;
}
//___________________________________________________________________________  
double PetersonFragm::Value(double z) const
{
// Evaluate the fragmentation function

  if(z<0 || z>1) return 0;
  return fFunc->Eval(z);
}
//___________________________________________________________________________
double PetersonFragm::GenerateZ(void) const
{
// Return a random number using the fragmentation function as PDF

  return fFunc->GetRandom();
}
//___________________________________________________________________________
void PetersonFragm::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->BuildFunction();
}
//___________________________________________________________________________
void PetersonFragm::Configure(string config)
{
  Algorithm::Configure(config);
  this->BuildFunction();
}
//___________________________________________________________________________
void PetersonFragm::BuildFunction(void) 
{
  fFunc = new TF1("fFunc",genie::utils::frgmfunc::peterson_func,0,1,2);

  fFunc->SetParNames("Norm","Epsilon");
  // stop ROOT from deleting this object of its own volition
  gROOT->GetListOfFunctions()->Remove(fFunc);

  double N = -1. ;
  GetParam( "PetFrag-Norm", N, false ) ;

  //The xml files says that Epsilon is not optional, while the
  //orignal behaviour of this function assigned a 0 default value
  //if the xml file was not providing such a number.
  // Original behaviour maintained
  double e = 0. ;
  GetParam( "PetFrag-Epsilon", e, false ) ;

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

