//____________________________________________________________________________
/*!

  \class    genie::NucleusGenerator

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

#ifndef _NUCLEUS_GENERATOR_H_
#define _NUCLEUS_GENERATOR_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Interaction/Target.h"
#include "Physics/NuclearState/SRCNuclearRecoil.h"
#include "Physics/NuclearState/SecondNucleonEmissionI.h"

namespace genie {


  class NucleusGenerator : public EventRecordVisitorI {

    public :
      NucleusGenerator();
      NucleusGenerator(string config);
      ~NucleusGenerator();

      //-- implement the EventRecordVisitorI interface
      void ProcessEventRecord(GHepRecord * event_rec) const;

      //-- overload the Algorithm::Configure() methods to load private data
      //   members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:

      void LoadConfig (void);

      const EventRecordVisitorI *fNucleusGen;
  };
}      // genie namespace
#endif // _NUCLEUS_GENERATOR_H_
