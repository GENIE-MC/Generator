/*!___________________________________________________________________________

\class    genie::RSHelicityAmplModel

\brief    Pure abstract base class. Defines the RSHelicityAmplModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 10, 2004

____________________________________________________________________________*/

#include "ReinSeghal/RSHelicityAmplModelI.h"

using namespace genie;

//____________________________________________________________________________
RSHelicityAmplModelI::RSHelicityAmplModelI() :
Algorithm()
{

}
//____________________________________________________________________________
RSHelicityAmplModelI::RSHelicityAmplModelI(const char * param_set) :
Algorithm(param_set)
{

}
//____________________________________________________________________________
RSHelicityAmplModelI::~RSHelicityAmplModelI()
{

}
//____________________________________________________________________________
