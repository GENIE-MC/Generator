//____________________________________________________________________________
/*!

\class   genie::EventGeneratorList

\brief   An abstract algorithmic class which can load a list of EventGenerator
         objects specified by its XML configuration file. \n
         To be subclassed by concrete algorithmic objects operating on a list
         of EventGenerators, like an Interaction Selector or the EventGenerator
         Chain of Responsibility.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created January 25, 2004

*/
//____________________________________________________________________________

#include "EVGCore/EventGeneratorList.h"
#include "EVGCore/EventGeneratorI.h"
#include "Messenger/Messenger.h"

using std::endl;
using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const EventGeneratorList & evgl)
 {
   evgl.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
EventGeneratorList::EventGeneratorList()
{

}
//___________________________________________________________________________
EventGeneratorList::~EventGeneratorList()
{

}
//___________________________________________________________________________
void EventGeneratorList::Print(ostream & stream) const
{
  EventGeneratorList::const_iterator iter;

  for(iter = this->begin(); iter != this->end(); ++iter) {

    const EventGeneratorI * evg = *iter;

    if(evg) stream << *evg;
    else    stream << "\n********* NULL EVENT GENERATOR *********" << endl;
  }
}
//___________________________________________________________________________
