//____________________________________________________________________________
/*!

\class    genie::SpectralFunctionLDA

\brief    A realistic spectral function computed using the Local Density
          Aproximation.

          Is a concrete implementation of the SpectralFunctionI interface.
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 07, 2004

*/ 
//____________________________________________________________________________

#include "Nuclear/SpectralFunctionLDA.h"

using namespace genie;

//____________________________________________________________________________
SpectralFunctionLDA::SpectralFunctionLDA() :
SpectralFunctionI()
{
  fName = "genie::SpectralFunctionLDA";
}
//____________________________________________________________________________
SpectralFunctionLDA::SpectralFunctionLDA(const char * param_set) :
SpectralFunctionI(param_set)
{
  fName = "genie::SpectralFunctionLDA";
}
//____________________________________________________________________________
SpectralFunctionLDA::~SpectralFunctionLDA()
{

}
//____________________________________________________________________________
double SpectralFunctionLDA::Prob(double P_nucleon, double E_nucleus)
{
  return 0;
}
//____________________________________________________________________________

