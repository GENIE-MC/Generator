//____________________________________________________________________________
/*!

  \class    genie::INCLNucleus

  \brief    INCLXX nuclear model. Implements the NuclearModelI 
  interface.

  \ref      

  \author   Liang Liu, (liangliu@fnal.gov)

  \created  Oct. 2024

  \cpright  Copyright (c) 2003-2024, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________


#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

#ifndef _INCL_NUCLEUS_H_
#define _INCL_NUCLEUS_H_

// ROOT
#include "TVector3.h"


// INCL++
// For configuration
#include "G4INCLConfig.hh"
#include "G4INCLNucleus.hh"
// GENIE
#include "Framework/GHEP/GHepRecord.h"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPropagationModel.hh"
#include "G4INCLStandardPropagationModel.hh"
#include "G4INCLCascadeAction.hh"
#include "G4INCLEventInfo.hh"
#include "G4INCLGlobalInfo.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"
#include "G4INCLRootFinder.hh"



namespace genie {

  class INCLNucleus {

    public: 
      static INCLNucleus * Instance (void);
      void initialize(const GHepRecord * evrec);
      void reset(const GHepRecord * evrec);

      TVector3 getHitNucleonPosition();
      TVector3 getHitNucleonMomentum();
      double   getHitNucleonEnergy();
      double   getHitNucleonMass();
      double   getMass();
      double   getRemovalEnergy();
      G4INCL::Nucleus * getNuclues();
      G4INCL::StandardPropagationModel * getPropagationModel();
      G4INCL::Particle *getHitParticle();
      G4INCL::Config *getConfig(){return theConfig_;}

    private:
      INCLNucleus();
      ~INCLNucleus();
      void init();

      static INCLNucleus *fInstance;

      TVector3 v3_; // position of initial nucleon 
      TVector3 p3_; // fermi momentum of initial nucleon
      double  energy_; // off-shell energy of initial nucleon
      G4INCL::Config *theConfig_;
      G4INCL::Nucleus *nucleus_;
      G4INCL::StandardPropagationModel *propagationModel_;
      G4INCL::CascadeAction *cascadeAction_;
      G4INCL::Particle *hitNucleon_;

      int nucleon_index_;

      double maxUniverseRadius_;
      double maxInteractionDistance_;
      double minRemnantSize_;


  };

}         // genie namespace
#endif    // _INCL_NUCLEUS_H_
#endif // __GENIE_INCL_ENABLED__
