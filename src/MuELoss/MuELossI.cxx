//____________________________________________________________________________
/*!

\class    genie::MuELossI

\brief    Cross Section Calculation Interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  December 10, 2003

*/
//____________________________________________________________________________

#include "MuELoss/MuELossI.h"

using namespace genie;
using namespace genie::mueloss;

//___________________________________________________________________________
MuELossI::MuELossI() :
Algorithm()
{

}
//___________________________________________________________________________
MuELossI::MuELossI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
MuELossI::MuELossI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
MuELossI::~MuELossI()
{

}
//___________________________________________________________________________

