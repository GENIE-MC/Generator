//____________________________________________________________________________
/*!

\class   genie::InteractionFilter

\brief   Filters interactions out of InteractionLists.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 15, 2005

*/
//____________________________________________________________________________

#include "EVGCore/InteractionFilter.h"
#include "Interaction/Interaction.h"

using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const InteractionFilter & filter)
 {
   filter.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
InteractionFilter::InteractionFilter()
{
  this->Initialize();
}
//___________________________________________________________________________
InteractionFilter::InteractionFilter(const InteractionFilter & filter)
{
  fIncludeCC  = filter.fIncludeCC;
  fIncludeNC  = filter.fIncludeNC;
  fIncludeEL  = filter.fIncludeEL;
  fIncludeQEL = filter.fIncludeQEL;
  fIncludeRES = filter.fIncludeRES;
  fIncludeDIS = filter.fIncludeDIS;
  fIncludeCOH = filter.fIncludeCOH;
}
//___________________________________________________________________________
InteractionFilter::~InteractionFilter()
{

}
//___________________________________________________________________________
bool InteractionFilter::FilterInteraction(
                                       const Interaction & interaction) const
{
  const ProcessInfo & proc_info = interaction.GetProcessInfo();

  if ( proc_info.IsElastic()       && !fIncludeEL  ) return true;
  if ( proc_info.IsQuasiElastic()  && !fIncludeQEL ) return true;
  if ( proc_info.IsDeepInelastic() && !fIncludeDIS ) return true;
  if ( proc_info.IsResonant()      && !fIncludeRES ) return true;
  if ( proc_info.IsCoherent()      && !fIncludeCOH ) return true;
  if ( proc_info.IsWeakCC()        && !fIncludeCC  ) return true;
  if ( proc_info.IsWeakNC()        && !fIncludeNC  ) return true;

  return false;
}
//___________________________________________________________________________
void InteractionFilter::Initialize()
{
  fIncludeCC  = true;
  fIncludeNC  = true;
  fIncludeEL  = true;
  fIncludeQEL = true;
  fIncludeRES = true;
  fIncludeDIS = true;
  fIncludeCOH = true;
}
//___________________________________________________________________________
void InteractionFilter::Print(ostream & stream) const
{
  stream << "\n Interaction Filter:";

  stream << "\n [-] Interaction types:";
  stream << "\n  |";
  stream << "\n  |-----o  Include CC ?.............." << fIncludeCC;
  stream << "\n  |-----o  Include NC ?.............." << fIncludeNC;
  stream << "\n";

  stream << "\n [-] Scatering types:";
  stream << "\n  |";
  stream << "\n  |-----o  Include EL  ?............." << fIncludeEL;
  stream << "\n  |-----o  Include QEL ?............." << fIncludeQEL;
  stream << "\n  |-----o  Include RES ?............." << fIncludeRES;
  stream << "\n  |-----o  Include DIS ?............." << fIncludeDIS;
  stream << "\n  |-----o  Include COH ?............." << fIncludeCOH;
  stream << "\n";
}
//___________________________________________________________________________
