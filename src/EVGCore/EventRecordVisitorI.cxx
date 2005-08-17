//____________________________________________________________________________
/*!

\class   genie::EventRecordVisitorI

\brief   Defines the EventRecordVisitorI interface.
         Concrete implementations of this interface use the 'Visitor' Design
         Pattern to perform an operation on an EventRecord.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 04, 2004

*/
//____________________________________________________________________________

#include "EVGCore/EventRecordVisitorI.h"

using namespace genie;

//___________________________________________________________________________
EventRecordVisitorI::EventRecordVisitorI() :
Algorithm()
{

}
//___________________________________________________________________________
EventRecordVisitorI::EventRecordVisitorI(const char * param_set) :
Algorithm(param_set)
{

}
//___________________________________________________________________________
EventRecordVisitorI::~EventRecordVisitorI()
{

}
//___________________________________________________________________________
