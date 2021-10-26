//____________________________________________________________________________
/*!
*/
//____________________________________________________________________________

#ifndef _PHOTON_COH_GENERATOR_H_
#define _PHOTON_COH_GENERATOR_H_

#include <TPythia6.h>
#include <TComplex.h>

#include "Framework/EventGen/EventRecordVisitorI.h"

namespace genie {

class PhotonCOHGenerator : public EventRecordVisitorI {

public :
  PhotonCOHGenerator();
  PhotonCOHGenerator(string config);
 ~PhotonCOHGenerator();

  // implement the EventRecordVisitorI interface
  void Initialize         (void)               const;
  void ProcessEventRecord (GHepRecord * evrec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);

  void Configure(string config);

private:

  void LoadConfig(void);

  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class

};

}      // genie namespace
#endif // _PHOTON_COH_GENERATOR_H_