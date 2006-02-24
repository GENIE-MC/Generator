//____________________________________________________________________________
/*!

\class    genie::MinimizerI

\brief    Numerical minimization/maximazation algorithm ABC

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#include "Numerical/MinimizerI.h"

using namespace genie;

//____________________________________________________________________________
MinimizerI::MinimizerI() :
Algorithm()
{

}
//____________________________________________________________________________
MinimizerI::MinimizerI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
MinimizerI::MinimizerI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
MinimizerI::~MinimizerI()
{

}
//____________________________________________________________________________
