//____________________________________________________________________________
/*!

\class   genie::RESResonanceSelector

\brief   Generates a baryon resonance for (v+N->Resonance->pi+X) events &
         adds it to the event record.

         Is a concrete implementation of the VtxGeneratorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#ifndef _RES_RESONANCE_SELECTOR_H_
#define _RES_RESONANCE_SELECTOR_H_

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResList.h"
#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class RESResonanceSelector : public EventRecordVisitorI {

public :

  RESResonanceSelector();
  RESResonanceSelector(const char * param_set);
  ~RESResonanceSelector();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  Resonance_t SelectResonance   (GHepRecord * event_rec) const;
  void        AddResonance      (GHepRecord * event_rec) const;
  void        FillResonanceList (BaryonResList& reslist) const;
  int         ResQ              (const Interaction * in) const;
};

}      // genie namespace

#endif // _RES_RESONANCE_SELECTOR_H_
