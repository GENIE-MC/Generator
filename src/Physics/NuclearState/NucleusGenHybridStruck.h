//____________________________________________________________________________
/*!

  \class    genie::NucleusGenHybridStruck

  \brief    It visits the event record & combines the FermoMover and VertexGenerator
  computes a Fermi motion momentum and position for initial state nucleons 
  bound in nuclei. It is configured in NucleusGenerator. This new interface
  is compatible with previous implments.
  Is a concrete implementation of the EventRecordVisitorI interface.

  \author   Liang Liu <liangliu \at fnal.gov>
  Fermi National Accelerator Laboratory

  \created  Jan 17, 2025

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________


#ifndef _NUCLEUS_GEN_TRADITIONAL_H_
#define _NUCLEUS_GEN_TRADITIONAL_H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Physics/NuclearState/NucleusGenI.h"
#include "Framework/GHEP/GHepParticle.h"

namespace genie {


  class NucleusGenHybridStruck : public NucleusGenI {

    public :
      NucleusGenHybridStruck();
      NucleusGenHybridStruck(string config);
      ~NucleusGenHybridStruck();

      //-- implement the EventRecordVisitorI interface
      void ProcessEventRecord(GHepRecord * event_rec) const;

      void GenerateCluster(GHepRecord *event_rec) const;
      void setInitialStateVertex   (GHepRecord * evrec) const;
      void setInitialStateMomentum   (GHepRecord * evrec) const;
      void BindHitNucleon(Interaction& interaction, double& Eb, QELEvGen_BindingMode_t hitNucleonBindingMode) const;
      void GenerateNucleon(Interaction* interaction, ResamplingHitNucleon_t resampling_mode) const;
      bool isRPValid(double r, double p, const Target & tgt) const;
      void SetHitNucleonOnShellMom(TVector3 p3) const;

      //-- overload the Algorithm::Configure() methods to load private data
      //   members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

    private:

      void LoadConfig (void);

      TVector3 GetVertex(Interaction* interaction) const;
      void setClusterVertex(GHepRecord * evrec) const;
      void setINCLVertex(GHepRecord * evrec) const;

      // Hybrid nuclear model combines the following 
      // 1. INCL vertex  + GENIE nuclear model + INCL FSI
      // 2. GENIE vertex + GENIE nuclear model + INCL FSI
      // 3. GENIE vertex + GENIE nuclear model + GENIE FSI
      // The third option should be identical with prior GENIE
      // nuclear model in genie
      const EventRecordVisitorI *fFermiMover;
      // vertex model in genie
      const EventRecordVisitorI *fVertexGenerator;
      // nuclear model of INCL
      const NucleusGenI *fNucleusGen; 

      bool fINCLVertex;
      bool fINCLFSI;
  };

}      // genie namespace

#endif // _NUCLEUS_GEN_TRADITIONAL_H_
