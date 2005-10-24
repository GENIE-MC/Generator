//____________________________________________________________________________
/*!

\class   genie::IMDTargetRemnantGenerator

\brief   Generates all the non - primarly lepton final state particles in v
         IMD interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 17, 2005

*/
//____________________________________________________________________________

#ifndef _IMD_TARGET_REMNANT_GENERATOR_H_
#define _IMD_TARGET_REMNANT_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class IMDTargetRemnantGenerator : public EventRecordVisitorI {

public :

  IMDTargetRemnantGenerator();
  IMDTargetRemnantGenerator(string config);
  ~IMDTargetRemnantGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord (GHepRecord * evrec) const;

private:
  void AddElectronNeutrino     (GHepRecord * evrec) const;
  void AddTargetNucleusRemnant (GHepRecord * evrec) const;
};

}      // genie namespace

#endif // _IMD_TARGET_REMNANT_GENERATOR_H_
