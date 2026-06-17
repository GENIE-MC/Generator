//____________________________________________________________________________
/*!

\class    genie::DMRESHadronicSystemGenerator

\brief    Generates the 'final state' hadronic system in DM RES interactions. (same as v RES)
          It adds the remnant nucleus (if any), the pre-selected resonance
          and the resonance decay products at the GHEP record.
          Unlike the SPP thread, in the RES thread the resonance is specified
          at the time an interaction is selected but its decay products not
          (semi-inclusive resonance reactions). The off the mass-shell baryon
          resonance is decayed using a phase space generator. All kinematically
          available decay channels are being used (not just 1 pi channels).
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Zachary W. Orr <zachorrt@colostate.edu>
          Colorado State University

          based on code by
          Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory


\created  November 23, 2004

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _DMRES_HADRONIC_SYSTEM_GENERATOR_H_
#define _DMRES_HADRONIC_SYSTEM_GENERATOR_H_

#include "Physics/Common/HadronicSystemGenerator.h"

namespace genie {

class Decayer;

class DMRESHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  DMRESHadronicSystemGenerator();
  DMRESHadronicSystemGenerator(string config);
 ~DMRESHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig                (void);
  int  GetResonancePdgCode       (GHepRecord * evrec) const;
  void AddResonance              (GHepRecord * evrec, int pdgc) const;
  // void AddResonanceDecayProducts (GHepRecord * evrec, int pdgc) const;

  const EventRecordVisitorI * fResonanceDecayer;
};

}      // genie namespace

#endif // _DMRES_HADRONIC_SYSTEM_GENERATOR_H_
