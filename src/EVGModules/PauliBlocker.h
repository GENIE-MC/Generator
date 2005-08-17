//____________________________________________________________________________
/*!

\class   genie::PauliBlocker

\brief   Examines whether the generated event should be Pauli blocked.

         Is a concerete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 08, 2004

*/
//____________________________________________________________________________

#ifndef _PAULI_BLOCKER_H_
#define _PAULI_BLOCKER_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class PauliBlocker : public EventRecordVisitorI {

public :

  PauliBlocker();
  PauliBlocker(const char * param_set);
  ~PauliBlocker();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _PAULI_BLOCKER_H_
