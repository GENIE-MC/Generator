//____________________________________________________________________________
/*!

\class    genie::PetersonFragm

\brief    The Peterson fragmentation function.

          Is a concrete implementation of the FragmentationFunctionI interface.
          
\ref      C.Peterson et al., Phys.Rev.D23, 56 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 15, 2004

*/ 
//____________________________________________________________________________

#include "Fragmentation/PetersonFragm.h"
#include "Fragmentation/fragmentation_functions.h"

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
double PetersonFragm::GenerateZ(void) const
{
  //-- get fragmentation function parameters from the config. registry
  
  double N = (fConfig->Exists("norm"))    ? fConfig->GetDouble("norm")    : 0;
  double e = (fConfig->Exists("epsilon")) ? fConfig->GetDouble("epsilon") : 0;

  fFunc->SetParameters(N,e);

  double Rndm = fFunc->GetRandom();

  return Rndm;
}
//___________________________________________________________________________
void PetersonFragm::BuildFunction(void) 
{
  fFunc = new TF1("fFunc",genie::peterson_fragmentation_function,0,1,2);

  fFunc->SetParNames("Norm","Epsilon");
}
//___________________________________________________________________________

