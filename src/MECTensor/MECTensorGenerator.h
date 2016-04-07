//____________________________________________________________________________
/*!

\class    genie::MECTensorGenerator

\brief    Simulate the primary MEC interaction from Nieves' MEC model

\author   Jackie Schwehr
           Heavily based on Steve Dytman and Costas Andreopoulous' MEC Generator

\created  2014/09/15

\notes    2015/01/14 all neutrino flavors, carbon and oxygen ~jschwehr


*/
//____________________________________________________________________________

#ifndef _MEC_TENSOR_GENERATOR_H_
#define _MEC_TENSOR_GENERATOR_H_

#include <TGenPhaseSpace.h>

#include "EVGCore/EventRecordVisitorI.h"
#include "PDG/PDGCodeList.h"


namespace genie {

class XSecAlgorithmI;
class NuclearModelI;

class MECTensorGenerator : public EventRecordVisitorI {

public :
  MECTensorGenerator();
  MECTensorGenerator(string config);
 ~MECTensorGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void        LoadConfig                  (void);
  void        AddNucleonCluster           (GHepRecord * event) const;
  void        AddTargetRemnant            (GHepRecord * event) const;
  void        SelectLeptonKinematics      (GHepRecord * event) const;
  void        GenerateInitialHadrons      (GHepRecord * event) const;
  //void        RecoilNucleonCluster        (GHepRecord * event) const;
  void        DecayNucleonCluster         (GHepRecord * event) const;
  PDGCodeList NucleonClusterConstituents  (int pdgc)           const;
  
  mutable const XSecAlgorithmI * fXSecModel;
  mutable TGenPhaseSpace         fPhaseSpaceGenerator;
  const NuclearModelI *          fNuclModel;
};

}      // genie namespace
#endif // _MEC_GENERATOR_H_
