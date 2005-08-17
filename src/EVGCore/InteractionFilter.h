//____________________________________________________________________________
/*!

\class   genie::InteractionFilter

\brief   Filters interactions out of InteractionLists.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 15, 2005

*/
//____________________________________________________________________________

#ifndef _INTERACTION_FILTER_H_
#define _INTERACTION_FILTER_H_

#include <ostream>

#include "Interaction/InitialState.h"

using std::ostream;

namespace genie {

class Interaction;

class InteractionFilter {

public :

  InteractionFilter();
  InteractionFilter(const InteractionFilter & filter);
  ~InteractionFilter();

  void SwitchChargedCurrent (bool on_off) { fIncludeCC  = on_off; }
  void SwitchNeutralCurrent (bool on_off) { fIncludeNC  = on_off; }
  void SwitchElastic        (bool on_off) { fIncludeEL  = on_off; }
  void SwitchQuasiElastic   (bool on_off) { fIncludeQEL = on_off; }
  void SwitchResonance      (bool on_off) { fIncludeRES = on_off; }
  void SwitchDeepInelastic  (bool on_off) { fIncludeDIS = on_off; }
  void SwitchCoherentPion   (bool on_off) { fIncludeCOH = on_off; }

  bool FilterInteraction (const Interaction & interaction) const;

  void Print(ostream & stream) const;

  friend ostream & operator<< (ostream& stream, const InteractionFilter & filter);
    
private:

  void Initialize();
  
  bool fIncludeCC;
  bool fIncludeNC;
  bool fIncludeEL;
  bool fIncludeQEL;
  bool fIncludeRES;
  bool fIncludeDIS;
  bool fIncludeCOH;
};

}      // genie namespace

#endif // _INTERACTION_FILTER_H_
