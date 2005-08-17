//____________________________________________________________________________
/*!

\class   genie::InteractionList

\brief

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 13, 2005

*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"

using std::endl;
using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const InteractionList & intl)
 {
   intl.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
InteractionList::InteractionList()
{

}
//___________________________________________________________________________
InteractionList::~InteractionList()
{
  InteractionList::const_iterator iter;

  for(iter = this->begin(); iter != this->end(); ++iter) {
    Interaction * interaction = *iter;
    delete interaction;
    interaction = 0;
  }
}
//___________________________________________________________________________
void InteractionList::Print(ostream & stream) const
{
  InteractionList::const_iterator iter;

  for(iter = this->begin(); iter != this->end(); ++iter) {

    Interaction * interaction = *iter;

    if(interaction) stream << *interaction;
    else            stream << "\n******** NULL INTERACTION ********" << endl;
  }
}
//___________________________________________________________________________
