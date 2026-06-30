//____________________________________________________________________________
/*!

  \class    genie::INCLNucleus

  \brief    INCLXX nuclear model. Implements the NuclearModelI 
  interface.

  \ref      

  \author   Liang Liu, (liangliu@fnal.gov)

  \created  Oct. 2024

  \cpright  Copyright (c) 2003-2026, The GENIE Collaboration
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
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Interaction/Target.h"
#include "Physics/NuclearState/NuclearModel.h"
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
#include "G4INCLParticleTable.hh"


namespace genie {

  class INCLNucleus {

    public: 
      static INCLNucleus * Instance (void);
      void reset(const Target * tgt);
      void configure();

      // TODO: these function is related to single nucleon, might need to refactor it
      TVector3 getHitNucleonPosition();
      TVector3 getHitNucleonMomentum();
      double   getHitNucleonEnergy();
      double   getHitNucleonMass();
      double   getMass();
      double   getRemovalEnergy();


      //============================================================
      G4INCL::Config *getConfig(){return theConfig_;}
      G4INCL::StandardPropagationModel * getPropagationModel();
      G4INCL::Nucleus * getNuclues();
      G4INCL::Particle *getHitParticle();
      std::shared_ptr<G4INCL::Cluster>  getHitNNCluster();
      double getMaxUniverseRadius() {return maxUniverseRadius_;}


      // void setHitParticle(const int pdg, TVector3 &posi);
      void setHitNNCluster(const int pdg1, const int pdg2, TVector3 &posi);

      // set INCL configurations
      void setINCLXXDataFilePath(std::string str){ INCLXXDataFilePath_ = str; }
      void setABLAXXDataFilePath(std::string str){ ablaxxDataFilePath_ = str; }
      void setABLA07DataFilePath(std::string str){ abla07DataFilePath_ = str; }
      void setGEMINIXXDataFilePath(std::string str){ geminixxDataFilePath_ = str; }
      void setDeExcitationType(G4INCL::DeExcitationType deExType){ deExcitationType_ = deExType; }
      void setPotentialType(G4INCL::PotentialType type) { potentialType_ = type; }
      void setPauliType(G4INCL::PauliType type) { pauliType_ = type; }
      void setPauliString(std::string str) {pauliString_ = str;}

      void setLocalEnergyBBType(G4INCL::LocalEnergyType type) {localEnergyTypeBB_ = type;}
      void setLocalEnergyPiType(G4INCL::LocalEnergyType type) {localEnergyTypePi_ = type;}
      void setHadronizationTime(const double t) { hadronizationTime_=t; }
      void setClusterAlgorithmType(const G4INCL::ClusterAlgorithmType c){
        clusterAlgorithmType_ = c;
      }
      void setClusterAlgorithmString(const std::string str){
        clusterAlgorithmString_ = str;
      }


      bool isRPValid(double r, double p);

      void ResamplingHitNucleon();
      TVector3 ResamplingVertex(const int pdg);

      void setHybridModel(NuclearModel_t model){
        model_type_ = model;
      }

      NuclearModel_t getHybridModel(){
        return model_type_;
      }

    private:
      INCLNucleus();
      ~INCLNucleus();
      void initialize(const Target * tgt);
      void initUniverseRadius(const int A, const int Z);
      G4INCL::Particle* getNucleon(const int pdg);
      std::shared_ptr<G4INCL::Cluster> getNNCluster(const int pdg1, const int pdg2);
      static INCLNucleus *fInstance;

      G4INCL::Config *theConfig_;
      G4INCL::Nucleus *nucleus_;
      G4INCL::Particle *hitNucleon_;
      // NN cluster for MEC channel
      std::shared_ptr<G4INCL::Cluster>  clusterNN_;

      // index of nucleon inside nucleus
      // 1p1h, 2p2h, might need NpNh
      int nucleon_index_;
      int cluster_index1_;
      int cluster_index2_;


      G4INCL::StandardPropagationModel *propagationModel_;
      // TODO: Using official GENIE action
      // G4INCL::CascadeAction *cascadeAction_;
      const G4INCL::NuclearDensity *theDensity;
      const G4INCL::NuclearPotential::INuclearPotential *thePotential;

      double maxUniverseRadius_;
      double minRemnantSize_;
      double hadronizationTime_;

      std::string INCLXXDataFilePath_;
      std::string abla07DataFilePath_;
      std::string ablaxxDataFilePath_;
      std::string geminixxDataFilePath_;

      G4INCL::DeExcitationType deExcitationType_;
      G4INCL::PotentialType potentialType_;
      G4INCL::PauliType pauliType_;
      std::string pauliString_;

      G4INCL::LocalEnergyType localEnergyTypeBB_;
      G4INCL::LocalEnergyType localEnergyTypePi_;
      std::string clusterAlgorithmString_;
      G4INCL::ClusterAlgorithmType clusterAlgorithmType_;

      NuclearModel_t model_type_;
  };

}         // genie namespace
#endif    // _INCL_NUCLEUS_H_
#endif // __GENIE_INCL_ENABLED__
