//____________________________________________________________________________
/*!

\class    genie::RSPPResonanceSelector

\brief    Generates an intermediate baryon resonance for exclusive interactions
          proceeding through resonance productions and adds it to the event
          record. The resonance is selected based on its contribution to the
          selected exclusive reaction cross section.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  November 18, 2004

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // override the Algorithm::Configure methods to load configuration
  // data to private data members
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
