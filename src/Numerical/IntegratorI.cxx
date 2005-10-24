//____________________________________________________________________________
/*!

\class    genie::IntegratorI

\brief    

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 01, 2004

*/
//____________________________________________________________________________

#include "Numerical/IntegratorI.h"

using namespace genie;

//____________________________________________________________________________
IntegratorI::IntegratorI() :
Algorithm()
{

}
//____________________________________________________________________________
IntegratorI::IntegratorI(string name) :
Algorithm(name)
{

}
//____________________________________________________________________________
IntegratorI::IntegratorI(string name, string config) :
Algorithm(name, config)
{

}
//____________________________________________________________________________
IntegratorI::~IntegratorI()
{

}
//____________________________________________________________________________
