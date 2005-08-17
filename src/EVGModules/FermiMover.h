//____________________________________________________________________________
/*!

\class   genie::FermiMover

\brief   It visits the event record & computes a Fermi motion momentum for
         initial state nucleons bound in nuclei.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 08, 2004

*/
//____________________________________________________________________________

#ifndef _FERMI_MOVER_H_
#define _FERMI_MOVER_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class FermiMover : public EventRecordVisitorI {

public :

  FermiMover();
  FermiMover(const char * param_set);
  ~FermiMover();

  //-- implement the EventRecordVisitorI interface

  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace

#endif // _FERMI_MOVER_H_
