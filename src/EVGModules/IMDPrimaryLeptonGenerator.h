//____________________________________________________________________________
/*!

\class   genie::IMDPrimaryLeptonGenerator

\brief   Generates the final state primary lepton in inverse muon decay (IMD)
         neutrino interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 13, 2005

*/
//____________________________________________________________________________

#ifndef _IMD_PRIMARY_LEPTON_GENERATOR_H_
#define _IMD_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class IMDPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :

  IMDPrimaryLeptonGenerator();
  IMDPrimaryLeptonGenerator(const char * param_set);
  ~IMDPrimaryLeptonGenerator();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _IMD_PRIMARY_LEPTON_GENERATOR_H_
