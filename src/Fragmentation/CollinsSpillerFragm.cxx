//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - June 15, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Fragmentation/CollinsSpillerFragm.h"
#include "Fragmentation/fragmentation_functions.h"

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
  //-- get fragmentation function parameters from the config. registry

  double N = (fConfig->Exists("norm"))    ? fConfig->GetDouble("norm")    : 0;
  double e = (fConfig->Exists("epsilon")) ? fConfig->GetDouble("epsilon") : 0;

  //-- evaluate the fragmentation function

  assert( z > 0 && z < 1);

  fFunc->SetParameters(N,e);

  double D = fFunc->Eval(z);

  return D;
}
//___________________________________________________________________________
double CollinsSpillerFragm::GenerateZ(void) const
{
  //-- get fragmentation function parameters from the config. registry

  double N = (fConfig->Exists("norm"))    ? fConfig->GetDouble("norm")    : 0;
  double e = (fConfig->Exists("epsilon")) ? fConfig->GetDouble("epsilon") : 0;

  fFunc->SetParameters(N,e);

  double Rndm = fFunc->GetRandom();

  return Rndm;
}
//___________________________________________________________________________
void CollinsSpillerFragm::BuildFunction(void)
{
  fFunc = new TF1("fFunc",genie::collins_spiller_fragmentation_function,0,1,2);

  fFunc->SetParNames("Norm","Epsilon");
}
//___________________________________________________________________________



