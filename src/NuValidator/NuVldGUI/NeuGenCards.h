//_____________________________________________________________________________
/*!

\class    genie::nuvld::NeuGenCards

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <C.V.Andreopoulos@rl.ac.uk>

\created  January 12, 2004
*/
//_____________________________________________________________________________

#ifndef _NEUGEN_CARDS_H_
#define _NEUGEN_CARDS_H_

#include <TObject.h>

#include "Facades/NeuGenConfig.h"
#include "Facades/NeuGenInputs.h"

using namespace genie::nuvld::facades;

namespace genie {
namespace nuvld {

class NeuGenCards : public TObject {

public:

   friend class NeuGenConfigDialog;
   friend class NeuGenInputsDialog;

   static NeuGenCards *  Instance (void);

   NeuGenInputs * CurrInputs (void) const { return _inputs; }
   NeuGenConfig * CurrConfig (void) const { return _config; }

   void SetInputs ( const NeuGenInputs * inputs );
   void SetConfig ( const NeuGenConfig * config );

private:

   NeuGenCards();
   NeuGenCards(const NeuGenCards & cards);
   ~NeuGenCards();

   void   SetDefaultConfig(void);

   static NeuGenCards * _self;

   NeuGenConfig *   _config;
   NeuGenInputs *   _inputs;

   ClassDef(NeuGenCards, 0)
};

} // nuvld namespace
} // genie namespace

#endif

