//____________________________________________________________________________
/*!

\class    genie::RSPPResonanceSelector

\brief    Generates an intermediate baryon resonance for exclusive interactions
          proceeding through resonance productions and adds it to the event
          record. The resonance is selected based on its contribution to the
          selected exclusive reaction cross section.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  November 18, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _RSPP_RESONANCE_SELECTOR_H_
#define _RSPP_RESONANCE_SELECTOR_H_

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResList.h"
#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class RSPPResonanceSelector : public HadronicSystemGenerator {

public :
  RSPPResonanceSelector();
  RSPPResonanceSelector(string config);
  ~RSPPResonanceSelector();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:
  void LoadConfig (void);

  Resonance_t SelectResonance   (GHepRecord * event_rec) const;
  void        AddResonance      (GHepRecord * event_rec) const;

  BaryonResList fResList; ///< baryon resonances taken into account
};

}      // genie namespace
#endif // _RSPP_RESONANCE_SELECTOR_H_
