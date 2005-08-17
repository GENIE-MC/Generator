//____________________________________________________________________________
/*!

\class   genie::RESHadronicSystemGenerator

\brief   Generates the 'final state' hadronic system in v RES interactions.

         It creates the GHepParticle entries for the target nucleus (if any)
         and the res. decay products and they are added to the GHEP record. \n

         The resonance decay products should be known at the time this visitor
         acts on th event record (this visitor should run on event generation
         threads initiated by selecting an exlusive RES channel).

         It does not handle the propagation of generated hadrons out of the
         nuclear medium and it does not handle decays of unstable particles
         (these would be handled by other event record visitors later on the
         event generation threads). So the 'final state' might not be final
         after all.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 23, 2004

*/
//____________________________________________________________________________

#ifndef _RES_HADRONIC_SYSTEM_GENERATOR_H_
#define _RES_HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class RESHadronicSystemGenerator : public HadronicSystemGenerator {

public :

  RESHadronicSystemGenerator();
  RESHadronicSystemGenerator(const char * param_set);
  ~RESHadronicSystemGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;

private:

  void AddResonanceDecayProducts (GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _RES_HADRONIC_SYSTEM_GENERATOR_H_
