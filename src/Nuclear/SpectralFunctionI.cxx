//____________________________________________________________________________
/*!

\class    genie::SpectralFunctionI

\brief    Pure abstract base class. Defines the SpectralFunctionI interface
          to be implemented by any algorithmic class implementing a spectral
          function.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include "Nuclear/SpectralFunctionI.h"

using namespace genie;

//___________________________________________________________________________
SpectralFunctionI::SpectralFunctionI() :
Algorithm()
{

}
//___________________________________________________________________________
SpectralFunctionI::SpectralFunctionI(const char * param_set) :
Algorithm(param_set)
{

}
//___________________________________________________________________________
SpectralFunctionI::~SpectralFunctionI()
{

}
//___________________________________________________________________________  
