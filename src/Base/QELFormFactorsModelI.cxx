//____________________________________________________________________________
/*!

\class    genie::QELFormFactorsModelI

\brief    Pure abstract base class. Defines the QELFormFactorsModelI interface
          to be implemented by any algorithmic class computing Quasi-Elastic
          Form Factors.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004
 
*/
//____________________________________________________________________________

#include "Base/QELFormFactorsModelI.h"

using namespace genie;

//____________________________________________________________________________
QELFormFactorsModelI::QELFormFactorsModelI() :
Algorithm()
{

}
//____________________________________________________________________________
QELFormFactorsModelI::QELFormFactorsModelI(const char * param_set) :
Algorithm(param_set)
{

}
//____________________________________________________________________________
QELFormFactorsModelI::~QELFormFactorsModelI()
{

}
//____________________________________________________________________________





