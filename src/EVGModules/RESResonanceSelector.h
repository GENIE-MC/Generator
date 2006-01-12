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

class XSecAlgorithmI;

class RESResonanceSelector : public EventRecordVisitorI {

public :

  RESResonanceSelector();
  RESResonanceSelector(string config);
  ~RESResonanceSelector();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadSubAlg     (void);
  void LoadConfigData (void);

  Resonance_t SelectResonance   (GHepRecord * event_rec) const;
  void        AddResonance      (GHepRecord * event_rec) const;
  int         ResQ              (const Interaction * in) const;

  BaryonResList          fResList;
  const XSecAlgorithmI * fXSecAlg;
};

}      // genie namespace

#endif // _RES_RESONANCE_SELECTOR_H_
