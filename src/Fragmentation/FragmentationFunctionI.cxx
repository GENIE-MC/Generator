//____________________________________________________________________________
/*!

\class    genie::FragmentationFunctionI

\brief    Pure abstract base class.
          Defines the FragmentationFunctionI interface to be implemented by
          any algorithmic class implementing a fragmentation function.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 15, 2004

*/ 
//____________________________________________________________________________

#include "Fragmentation/FragmentationFunctionI.h"

using namespace genie;

//___________________________________________________________________________
FragmentationFunctionI::FragmentationFunctionI() :
Algorithm()
{

}
//___________________________________________________________________________
FragmentationFunctionI::FragmentationFunctionI(const char * param_set) :
Algorithm(param_set)
{

}
//___________________________________________________________________________
FragmentationFunctionI::~FragmentationFunctionI()
{

}
//___________________________________________________________________________  
