//____________________________________________________________________________
/*!

\class    genie::ELFormFactorsModelI

\brief    Pure abstract base class. Defines the ELFormFactorsModelI interface
          to be implemented by any algorithmic class computing Elastic Form
          Factors.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "Base/ELFormFactorsModelI.h"

using namespace genie;

//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI() :
Algorithm()
{

}
//____________________________________________________________________________
ELFormFactorsModelI::ELFormFactorsModelI(const char * param_set) :
Algorithm(param_set)
{

}
//____________________________________________________________________________
ELFormFactorsModelI::~ELFormFactorsModelI()
{

}
//____________________________________________________________________________





