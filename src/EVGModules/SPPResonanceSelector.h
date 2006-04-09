//____________________________________________________________________________
/*!

\class   genie::SPPResonanceSelector

\brief   Generates an intermediate baryon resonance for exclusive interactions
         proceeding through resonance productions and adds it to the event
         record. The resonance is selected based on its contribution to the
         selected exclusive reaction cross section.
         Is a concrete implementation of the VtxGeneratorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#ifndef _SPP_RESONANCE_SELECTOR_H_
#define _SPP_RESONANCE_SELECTOR_H_

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResList.h"
#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class XSecAlgorithmI;

class SPPResonanceSelector : public HadronicSystemGenerator {

public :
  SPPResonanceSelector();
  SPPResonanceSelector(string config);
  ~SPPResonanceSelector();

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

  BaryonResList          fResList;
  const XSecAlgorithmI * fXSecAlg;
};

}      // genie namespace

#endif // _SPP_RESONANCE_SELECTOR_H_
