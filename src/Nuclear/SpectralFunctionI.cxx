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
SpectralFunctionI::SpectralFunctionI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
SpectralFunctionI::SpectralFunctionI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
SpectralFunctionI::~SpectralFunctionI()
{

}
//___________________________________________________________________________  
