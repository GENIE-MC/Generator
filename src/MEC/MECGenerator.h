//____________________________________________________________________________
/*!

\class    genie::MECGenerator

\brief    

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Sep. 22, 2008

\cpright  Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_GENERATOR_H_
#define _MEC_GENERATOR_H_

#include "EVGCore/EventRecordVisitorI.h"

namespace genie {

class XSecAlgorithmI;

class MECGenerator : public EventRecordVisitorI {

public :
  MECGenerator();
  MECGenerator(string config);
 ~MECGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig          (void);
  void AddFinalStateLepton (GHepRecord * event_rec) const;
  void SelectKinematics    (GHepRecord * event_rec) const;
  void AddNucleonCluster   (GHepRecord * event_rec) const;
  void AddTargetRemnant    (GHepRecord * event_rec) const;
  void DecayNucleonCluster (GHepRecord * event_rec) const;

  mutable const XSecAlgorithmI * fXSecModel;
};

}      // genie namespace
#endif // _MEC_GENERATOR_H_
