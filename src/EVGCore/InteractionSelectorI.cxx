//____________________________________________________________________________
/*!

\class   genie::InteractionSelectorI

\brief   Defines the InteractionSelectorI interface to be implemented by
         algorithms selecting interactions to be generated.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 05, 2004

*/
//____________________________________________________________________________

#include "EVGCore/InteractionSelectorI.h"
#include "Interaction/Interaction.h"

using namespace genie;

//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI() :
Algorithm()
{

}
//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
InteractionSelectorI::InteractionSelectorI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
InteractionSelectorI::~InteractionSelectorI()
{

}
//___________________________________________________________________________
