//____________________________________________________________________________
/*!

\class    genie::MultiplicityProbModelI

\brief    Pure abstract base class.
          Defines the MultiplicityProbModelI interface to be implemented by
          any algorithmic class computing multiplicity probability
          distributions for a hadronization model.

          Used in a Strategy Pattern together with the MultProbDistribution
          class.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

*/
//____________________________________________________________________________

#include "Fragmentation/MultiplicityProbModelI.h"

using namespace genie;

//____________________________________________________________________________
MultiplicityProbModelI::MultiplicityProbModelI() :
Algorithm()
{

}
//____________________________________________________________________________
MultiplicityProbModelI::MultiplicityProbModelI(const char * param_set) :
Algorithm(param_set)
{

}
//____________________________________________________________________________
MultiplicityProbModelI::~MultiplicityProbModelI()
{

}
//____________________________________________________________________________

