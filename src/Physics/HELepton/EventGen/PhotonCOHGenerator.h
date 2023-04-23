//____________________________________________________________________________
/*!

\class    genie::PhotonCOHGenerator

\brief    Generator for W boson production.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC & Harvard University

\created  Dec 8, 2021

\cpright  Copyright (c) 2003-2023, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
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