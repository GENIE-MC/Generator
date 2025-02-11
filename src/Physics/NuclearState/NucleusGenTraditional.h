//____________________________________________________________________________
/*!

\class    genie::NucleusGenTraditional

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  October 17, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________


#ifndef _NUCLEUS_GEN_TRADITIONAL_H_
#define _NUCLEUS_GEN_TRADITIONAL_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"

namespace genie {


class NucleusGenTraditional : public EventRecordVisitorI {

public :
  NucleusGenTraditional();
  NucleusGenTraditional(string config);
 ~NucleusGenTraditional();

  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfig (void);

  const EventRecordVisitorI *fFermiMover;
  const EventRecordVisitorI *fVertexGenerator;

};

}      // genie namespace

#endif // _NUCLEUS_GEN_TRADITIONAL_H_
