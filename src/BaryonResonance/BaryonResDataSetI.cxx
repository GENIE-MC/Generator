//____________________________________________________________________________
/*!

\class    genie::BaryonResDataSetI

\brief    Pure abstract base class. Defines the BaryonResDataSetI interface
          to be implemented by any algorithmic class computing (or merely
          retrieving) baryon resonance data.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BaryonResDataSetI.h"

using namespace genie;

//____________________________________________________________________________
BaryonResDataSetI::BaryonResDataSetI() :
Algorithm()
{

}
//____________________________________________________________________________
BaryonResDataSetI::BaryonResDataSetI(const char * param_set) :
Algorithm(param_set)
{

}
//____________________________________________________________________________
BaryonResDataSetI::~BaryonResDataSetI()
{

}
//____________________________________________________________________________





