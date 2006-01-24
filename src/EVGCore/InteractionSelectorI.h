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

#ifndef _INTERACTION_SELECTOR_I_H_
#define _INTERACTION_SELECTOR_I_H_

#include "Algorithm/Algorithm.h"

class TLorentzVector;

namespace genie {

class XSecAlgorithmMap;
class EventRecord;

class InteractionSelectorI : public Algorithm {

public :

  virtual ~InteractionSelectorI();

  //!  Define the InteractionSelectorI interface
  virtual EventRecord * SelectInteraction
     (const XSecAlgorithmMap * xsmp, const TLorentzVector & p4) const = 0;

protected:

  InteractionSelectorI();
  InteractionSelectorI(string name);
  InteractionSelectorI(string name, string config);
};

}      // genie namespace

#endif // _INTERACTION_SELECTOR_I_H_
