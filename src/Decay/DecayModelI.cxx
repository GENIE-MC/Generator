//____________________________________________________________________________
/*!

\class    genie::DecayModelI

\brief    Pure abstract base class. Defines the DecayModelI interface to be
          implemented by any algorithmic class decaying a particle.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 20, 2004

*/
//____________________________________________________________________________

#include "Decay/DecayModelI.h"

using namespace genie;

//____________________________________________________________________________
DecayModelI::DecayModelI() :
Algorithm()
{

}
//____________________________________________________________________________
DecayModelI::DecayModelI(const char * param_set) :
Algorithm(param_set)
{

}
//____________________________________________________________________________
DecayModelI::~DecayModelI() 
{ 

}
//____________________________________________________________________________



