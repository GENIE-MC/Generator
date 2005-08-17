//____________________________________________________________________________
/*!

\class   genie::EventGeneratorI

\brief   Defines the EventGeneratorI interface.

         The concrete implementations of this interface are Event Record
         Visitors (subclasses of the EventRecordVisitorI pABC) that,
         additionally, declare a 'Validity Context'. \n

         The declared validity context is used for selecting the appropriate
         concrete EventGeneratorI to generate the interacion at hand using
         the 'chain-of-responsibility' design pattern.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 22, 2004

*/
//____________________________________________________________________________

#include "EVGCore/EventGeneratorI.h"

using namespace genie;

//___________________________________________________________________________
EventGeneratorI::EventGeneratorI() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
EventGeneratorI::EventGeneratorI(const char * param_set) :
EventRecordVisitorI(param_set)
{

}
//___________________________________________________________________________
EventGeneratorI::~EventGeneratorI()
{

}
//___________________________________________________________________________
