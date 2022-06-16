//----------------------------------------------------------------------------
/*!

  This is a module for GENIE to read in (gnumi) flux ntuples and construct NHL fluxes
  on the fly. 
  
  Core loop: + Open dk2nu entry
             + Get ancestry information
	     + Assume decay to NHL (discard SM nu info, is unneeded)
	     + Calculate NHL production mode based on parameter space read from config
	     + Calculate kinematics of NHL
	     + Return NHL as SimpleNHL object.

\namespace  genie::NHL::NHLFluxCreator

\brief      Calculates NHL production kinematics & vertex

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    April 25th, 2022

\cpright    ??? - TBD

 */
//----------------------------------------------------------------------------
/*
  TODO: Make BBox from geometry file! (unit BBox if no geom-file?)
        Fix fLT (see above)
        Fix zm in case !parentHistCentre && IPdev.Mag() < fLT/2 (in this case zm = 0.0)
	  ==> parentHitsCentre should become parentInDetector
 */
//----------------------------------------------------------------------------

#ifndef _NHL_FLUXCREATOR_H_
#define _NHL_FLUXCREATOR_H_

// -- C++ includes
#include <array>
#include <cassert>
#include <map>
#include <sstream>

// -- ROOT includes
#include "TChain.h"
#include "TDecayChannel.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

// -- GENIE includes
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

#include "Tools/Flux/GNuMIFlux.h"

#include "Physics/NeutralHeavyLepton/NHLBRFunctions.h"
#include "Physics/NeutralHeavyLepton/NHLDecayVolume.h"
#include "Physics/NeutralHeavyLepton/NHLDecayUtils.h"
#include "Physics/NeutralHeavyLepton/NHLEnums.h"
#include "Physics/NeutralHeavyLepton/NHLKinUtils.h"
#include "Physics/NeutralHeavyLepton/SimpleNHL.h"

#endif // #ifndef _NHL_FLUXCREATOR_H_

const double kRDET = 1.0; // calculate fluxes per m^2

namespace genie{

  namespace NHL{
    
    class SimpleNHL;
    
    namespace NHLFluxCreator{

      // workhorse method
      void MakeTupleFluxEntry( int iEntry, genie::flux::GNuMIFluxPassThroughInfo * gnmf, std::string finpath, int run );
      void FillNonsense( int iEntry, genie::flux::GNuMIFluxPassThroughInfo * gnmf, int run );

      // init
      void OpenFluxInput( std::string finpath );
      void InitialiseTree();
      void InitialiseMeta();
      //RETHERE Should produce histograms for validation!
      void InitialiseHists();

      // Custom NHL kinematics, POT scaling, production probabilities
      double ScalePOT( double sm_pot );
      // returns NHL 4-momentum from random decay in same frame as p4par
      TLorentzVector NHLEnergy( genie::NHL::NHLProd_t nhldm, TLorentzVector p4par );
      // gets random point in BBox and returns separation to it in BEAM FRAME
      TVector3 PointToRandomPointInBBox( TVector3 detO_beam );
      TVector3 GetBoostBetaVec( TLorentzVector parp4 );

      void ReadBRs();
      std::map< genie::NHL::NHLProd_t, double > GetProductionProbs( int parPDG );
      
      // Obtain detector dimensions + position
      // RETHERE: BBox isn't good enough! But roll with it for now
      void MakeBBox();
      TVector3 ApplyUserRotation( TVector3 vec, bool doBackwards = false );
      
      // calculate detector acceptance (== solid angle of proj of det onto unit-radius sphere / (4pi))
      // NOTE THIS IS A LAB FRAME (==GEOMETRICAL) ACCEPTANCE!!!!
      // detO == detector BBox centre wrt NHL prod vertex, L{x,y,z} BBox length on each axis. Both [m]
      double CalculateDetectorAcceptanceSAA( TVector3 detO );
      // collimation effect calc, returns NHL_acc / geom_acc
      double CalculateAcceptanceCorrection( TLorentzVector p4par, TLorentzVector p4NHL, double SMECM, double zm, double zp );
      double labangle( double * x, double * par ); // function formula for correction
      // get minimum and maximum deviation from parent momentum to hit detector, [deg]
      double GetAngDeviation( TLorentzVector p4par, TVector3 detO, bool seekingMax );
      // returns 1.0 / (area of flux calc)
      double CalculateAreaNormalisation();

      /*
      // Obtain translation from beam to detector frame
      std::vector< double > GetBeam2UserTranslation();
      // Obtain rotation from beam to detector frame (Euler angles, extrinsic x-z-x)
      std::vector< double > GetBeam2UserRotation();
      // package the above two in a single call
      void GetBeam2UserTransformation();
      */

      // get delay on top of a SM neutrino.
      double GetTimeDelay(); //args?? 

      double GetPOTFromMeta( TChain * cmeta );
      void   LoopEntries( TChain * cflux, TChain * cmeta );

      // current path to keep track of what is loaded
      extern std::string fCurrPath;
      
      // maps to keep P( production )
      extern std::map< genie::NHL::NHLProd_t, double > dynamicScores; // map in use
      extern std::map< genie::NHL::NHLProd_t, double > dynamicScores_pion;
      extern std::map< genie::NHL::NHLProd_t, double > dynamicScores_kaon;
      extern std::map< genie::NHL::NHLProd_t, double > dynamicScores_muon;
      extern std::map< genie::NHL::NHLProd_t, double > dynamicScores_neuk;
      
      extern double BR_pi2mu, BR_pi2e, BR_K2mu, BR_K2e, BR_K3mu, BR_K3e, BR_K03mu, BR_K03e;

      extern bool doProduceHists;
      extern bool isParentOnAxis;

      extern TFile * fin;
      extern TChain * ctree, * cmeta;
      extern bool isTreeInit, isMetaInit, isBoxInit;

      extern double fLx, fLy, fLz;   //BBox length [m]
      extern double fCx, fCy, fCz;   //BBox centre wrt NHL prod [m]
      extern double fAx1, fAz, fAx2; //Euler angles, extrinsic x-z-x. Ax2 then Az then Ax1 
      extern double fDx, fDy, fDz;   //NHL production point [m]

      //std::vector< double > trVec, roVec;

      extern double parentMass, parentMomentum, parentEnergy; // GeV

      extern TGenPhaseSpace fPhaseSpaceGenerator;

      // functions for acceptance correction
      extern TF1 * fNHL, * fSMv;

      // tree variables. Add as per necessary.
      extern double potnum;                             ///< N POT for this SM-v
      extern int    decay_ptype;                        ///< PDG code of parent
      extern double decay_vx, decay_vy, decay_vz;       ///< coordinates of prod vtx [cm]
      extern double decay_pdpx, decay_pdpy, decay_pdpz; ///< final parent momentum [GeV]
      extern double decay_necm;                         ///< SM v CM energy [GeV]
      extern double decay_nimpwt;                       ///< Importance weight from beamsim

      // meta variables. Add as necessary
      extern int    job;                                ///< beamsim MC job number
      extern double pots;                               ///< how many pot in this job?

    } // namespace NHLFluxCreator
      
  } // namespace NHL
} // namespace genie
