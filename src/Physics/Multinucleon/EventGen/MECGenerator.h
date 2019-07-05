//____________________________________________________________________________
/*!

\class    genie::MECGenerator

\brief    Simulate the primary MEC interaction

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Steve Dytman <dytman+ \at pitt.edu>
          Pittsburgh University

\created  Sep. 22, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _MEC_GENERATOR_H_
#define _MEC_GENERATOR_H_

#include <TGenPhaseSpace.h>

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/ParticleData/PDGCodeList.h"

namespace genie {

class XSecAlgorithmI;
class NuclearModelI;

class MECGenerator : public EventRecordVisitorI {

public :
  MECGenerator();
  MECGenerator(string config);
 ~MECGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void    LoadConfig                        (void);
  void    AddNucleonCluster                 (GHepRecord * event) const;
  void    AddTargetRemnant                  (GHepRecord * event) const;
  void    GenerateFermiMomentum             (GHepRecord * event) const;
  void    SelectEmpiricalKinematics         (GHepRecord * event) const;
  void    AddFinalStateLepton               (GHepRecord * event) const;
  void    RecoilNucleonCluster              (GHepRecord * event) const;
  void    DecayNucleonCluster               (GHepRecord * event) const;
  void    SelectNSVLeptonKinematics         (GHepRecord * event) const;
  void    GenerateNSVInitialHadrons         (GHepRecord * event) const;
  PDGCodeList NucleonClusterConstituents    (int pdgc)           const;
  
  mutable const XSecAlgorithmI * fXSecModel;
  mutable TGenPhaseSpace         fPhaseSpaceGenerator;
  const NuclearModelI *          fNuclModel;

  double fQ3Max;
};

}      // genie namespace
#endif // _MEC_GENERATOR_H_
