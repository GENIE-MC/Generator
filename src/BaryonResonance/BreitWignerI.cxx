//____________________________________________________________________________
/*!

\class    genie::BreitWignerI

\brief    Pure abstract base class. Defines the BreitWignerI interface to
          be implemented by any algorithmic class modeling a Breit Wigner
          function.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 04, 2004

*/
//____________________________________________________________________________

#include "BaryonResonance/BreitWignerI.h"

using namespace genie;

//______________________________________________________________________
BreitWignerI::BreitWignerI() : Algorithm()
{

}
//______________________________________________________________________
BreitWignerI::BreitWignerI(const char * param_set) :
Algorithm(param_set)
{

}
//______________________________________________________________________
BreitWignerI::~BreitWignerI()
{

}
//______________________________________________________________________

