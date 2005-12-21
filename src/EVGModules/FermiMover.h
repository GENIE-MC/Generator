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

class NuclearPDistributionModelI;

class FermiMover : public EventRecordVisitorI {

public :

  FermiMover();
  FermiMover(string config);
  ~FermiMover();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  const NuclearPDistributionModelI *  fNuclPModel;
};

}      // genie namespace

#endif // _FERMI_MOVER_H_
