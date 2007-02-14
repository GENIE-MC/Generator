//____________________________________________________________________________
/*!

\class    genie::PauliBlocker

\brief    Examines whether the generated event should be Pauli blocked.
          Is a concerete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 08, 2004

\cpright  Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PAULI_BLOCKER_H_
#define _PAULI_BLOCKER_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class FermiMomentumTable;

class PauliBlocker : public EventRecordVisitorI {

public :
  PauliBlocker();
  PauliBlocker(string config);
 ~PauliBlocker();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

   void LoadKFTable(void);

   const FermiMomentumTable * fKFTable;
   string fKFTableName;
};

}      // genie namespace

#endif // _PAULI_BLOCKER_H_
