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
RSHelicityAmplModelI::RSHelicityAmplModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
RSHelicityAmplModelI::RSHelicityAmplModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
RSHelicityAmplModelI::~RSHelicityAmplModelI()
{

}
//____________________________________________________________________________
