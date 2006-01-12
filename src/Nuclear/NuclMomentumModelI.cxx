//____________________________________________________________________________
/*!

\class    genie::NuclMomentumModelI

\brief    Pure abstract base class.
          Defines the NuclMomentumModelI interface to be implemented by
          any algorithmic class generating momenta for nucleons within nuclei

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#include "Nuclear/NuclMomentumModelI.h"

using namespace genie;

//____________________________________________________________________________
NuclMomentumModelI::NuclMomentumModelI() :
Algorithm()
{

}
//____________________________________________________________________________
NuclMomentumModelI::NuclMomentumModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
NuclMomentumModelI::NuclMomentumModelI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
NuclMomentumModelI::~NuclMomentumModelI()
{

}
//____________________________________________________________________________

