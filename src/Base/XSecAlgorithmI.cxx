//____________________________________________________________________________
/*!

\class    genie::XSecAlgorithmI

\brief    Cross Section Calculation Interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#include "Base/XSecAlgorithmI.h"

using namespace genie;

//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI() :
Algorithm()
{

}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(const char * param_set) :
Algorithm(param_set)
{

}
//___________________________________________________________________________
XSecAlgorithmI::~XSecAlgorithmI()
{

}
//___________________________________________________________________________

