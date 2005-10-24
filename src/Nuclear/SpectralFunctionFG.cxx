//____________________________________________________________________________
/*!

\class    genie::SpectralFunctionFG

\brief    A Fermi Gas equivalent spectral function.

          Is a concrete implementation of the SpectralFunctionI interface.
          
\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 07, 2004

*/
//____________________________________________________________________________

#include "Nuclear/SpectralFunctionFG.h"

using namespace genie;

//____________________________________________________________________________
SpectralFunctionFG::SpectralFunctionFG() :
SpectralFunctionI("genie::SpectralFunctionFG")
{

}
//____________________________________________________________________________
SpectralFunctionFG::SpectralFunctionFG(string config) :
SpectralFunctionI("genie::SpectralFunctionFG", config)
{

}
//____________________________________________________________________________
SpectralFunctionFG::~SpectralFunctionFG()
{

}
//____________________________________________________________________________
double SpectralFunctionFG::Prob(double P_nucleon, double E_nucleus)
{
  return 0;
}
//____________________________________________________________________________

