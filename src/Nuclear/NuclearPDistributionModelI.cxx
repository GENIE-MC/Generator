//____________________________________________________________________________
/*!

\class    genie::NuclearPDistributionModelI

\brief    Pure abstract base class.
          Defines the NuclearPDistributionModelI interface to be implemented by
          any algorithmic class generating momenta for nucleons within nuclei

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#include "Nuclear/NuclearPDistributionModelI.h"

using namespace genie;

//____________________________________________________________________________
NuclearPDistributionModelI::NuclearPDistributionModelI() :
Algorithm()
{

}
//____________________________________________________________________________
NuclearPDistributionModelI::NuclearPDistributionModelI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
NuclearPDistributionModelI::NuclearPDistributionModelI(
                                                string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
NuclearPDistributionModelI::~NuclearPDistributionModelI()
{

}
//____________________________________________________________________________

