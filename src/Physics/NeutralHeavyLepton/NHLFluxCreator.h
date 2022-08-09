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

  From 23/Jul/2022: Merge NHLFluxReader into this class.
  Reads in flux histograms and returns appropriate stuff.

\class      genie::NHL::NHLFluxCreator

\brief      Calculates NHL production kinematics & vertex.
            Is a concrete implementation of the EventRecordVisitorI interface

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
#include <iomanip> // for momentum balance stream
#include <map>
#include <sstream>

// -- ROOT includes
#include "TChain.h"
#include "TDecayChannel.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TH1.h"
#include "TH3.h"
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
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
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
    
    class NHLFluxCreator : public EventRecordVisitorI {

    public:

      NHLFluxCreator();
      NHLFluxCreator(string config);
      ~NHLFluxCreator();

      //-- implement the EventRecordVisitorI interface
      void ProcessEventRecord(GHepRecord * event_rec) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);

      // set input path
      void SetInputPath( std::string finpath ) const;
      // get N(flux input entries)
      int GetNEntries() const;

      void SetFirstEntry( int iFirst ) const;
      void SetGeomFile( string geomfile ) const;
      void ImportBoundingBox( TGeoBBox * box ) const;

      // FluxReader-inherited functions
      // only if not using dk2nu!
      int SelectMass( const double mN ) const;
      std::string SelectFile( std::string fin, const double mN ) const; // find but don't open the file

      /// get the histogram and energy from it
      TH1F * GetFluxHist1F( std::string fin, int masspoint, bool isParticle ) const;
      /// sample production vertex from this histogram
      TH3D * GetFluxHist3D( std::string fin, std::string dirName, std::string hName ) const;
      std::vector< double > * GenerateVtx3X( TH3D * prodVtxHist ) const;

      flux::GNuMIFluxPassThroughInfo * RetrieveGNuMIFluxPassThroughInfo() const;

    private:

      void LoadConfig(void);

      // workhorse methods
      void MakeTupleFluxEntry( int iEntry, genie::flux::GNuMIFluxPassThroughInfo * gnmf, std::string finpath ) const;
      void FillNonsense( int iEntry, genie::flux::GNuMIFluxPassThroughInfo * gnmf ) const;

      // init
      void OpenFluxInput( std::string finpath ) const;
      void InitialiseTree() const;
      void InitialiseMeta() const;

      // returns NHL 4-momentum from random decay in same frame as p4par
      TLorentzVector NHLEnergy( genie::NHL::NHLProd_t nhldm, TLorentzVector p4par ) const;
      // gets random point in BBox and returns separation to it in BEAM FRAME
      TVector3 PointToRandomPointInBBox( TVector3 detO_beam ) const;

      void ReadBRs() const;
      std::map< genie::NHL::NHLProd_t, double > GetProductionProbs( int parPDG ) const;
      
      // Obtain detector dimensions + position
      // RETHERE: BBox isn't good enough! But roll with it for now
      void MakeBBox() const;
      TVector3 ApplyUserRotation( TVector3 vec, bool doBackwards = false ) const;
      TVector3 ApplyUserRotation( TVector3 vec, TVector3 oriVec, std::vector<double> rotVec, bool doBackwards = false ) const;
      
      // calculate detector acceptance (== solid angle of proj of det onto unit-radius sphere / (4pi))
      // NOTE THIS IS A LAB FRAME (==GEOMETRICAL) ACCEPTANCE!!!!
      // detO == detector BBox centre wrt NHL prod vertex, L{x,y,z} BBox length on each axis. Both [m]
      double CalculateDetectorAcceptanceSAA( TVector3 detO ) const;
      // collimation effect calc, returns NHL_acc / geom_acc
      double CalculateAcceptanceCorrection( TLorentzVector p4par, TLorentzVector p4NHL, double SMECM, double zm, double zp ) const;
      static double labangle( double * x, double * par ); // function formula for correction
      // get minimum and maximum deviation from parent momentum to hit detector, [deg]
      double GetAngDeviation( TLorentzVector p4par, TVector3 detO, bool seekingMax ) const;
      void GetAngDeviation( TLorentzVector p4par, TVector3 detO, TGeoManager * gm, double &zm, double &zp ) const;
      // returns 1.0 / (area of flux calc)
      double CalculateAreaNormalisation();

      // current path to keep track of what is loaded
      mutable std::string fCurrPath = "";
      // and which entry we're on
      mutable int iCurrEntry = 0;
      // which one was first?
      mutable int fFirstEntry = 0;
      // out of how many?
      mutable int fNEntries = 0;
      
      // maps to keep P( production )
      mutable std::map< genie::NHL::NHLProd_t, double > dynamicScores; // map in use
      mutable std::map< genie::NHL::NHLProd_t, double > dynamicScores_pion;
      mutable std::map< genie::NHL::NHLProd_t, double > dynamicScores_kaon;
      mutable std::map< genie::NHL::NHLProd_t, double > dynamicScores_muon;
      mutable std::map< genie::NHL::NHLProd_t, double > dynamicScores_neuk;
      
      mutable double BR_pi2mu, BR_pi2e, BR_K2mu, BR_K2e, BR_K3mu, BR_K3e, BR_K03mu, BR_K03e;

      mutable bool isParentOnAxis = true;
      mutable bool fUseBeamMomentum = false; // use this if your detector hall is parallel to tgt hall
      mutable TGeoManager * fGeoManager = 0;
      mutable TGeoVolume * fTopVol = 0;
      mutable string fGeomFile;

      mutable TChain * ctree = 0, * cmeta = 0;

      mutable double fMass; // NHL mass, GeV
      mutable std::vector< double > fU4l2s; // couplings
      
      mutable double fLx, fLy, fLz;
      mutable double fLxR, fLyR, fLzR; // BBox side [m]

      mutable std::vector< double > fB2UTranslation, fB2URotation;
      mutable std::vector< double > fDetRotation; // rotation of detector wrt tgt hall
      mutable std::vector< double > fDetOffset; // offset of det centre wrt geom file origin
      mutable double fCx, fCy, fCz;   // BBox centre wrt NHL prod [m]
      mutable double fAx1, fAz, fAx2; // Euler angles, extrinsic x-z-x. Tgt-hall to beam
      mutable double fBx1, fBz, fBx2; // Tgt-hall to detector frame

      mutable double fDx, fDy, fDz; //NHL production point [m]

      //std::vector< double > trVec, roVec;

      mutable double parentMass, parentMomentum, parentEnergy; // GeV

      // tree variables. Add as per necessary.
      mutable double potnum;                             ///< N POT for this SM-v
      mutable int    decay_ptype;                        ///< PDG code of parent
      mutable double decay_vx, decay_vy, decay_vz;       ///< coordinates of prod vtx [cm]
      mutable double decay_pdpx, decay_pdpy, decay_pdpz; ///< final parent momentum [GeV]
      mutable double decay_necm;                         ///< SM v CM energy [GeV]
      mutable double decay_nimpwt;                       ///< Importance weight from beamsim

      // meta variables. Add as necessary
      mutable int    job;                                ///< beamsim MC job number
      mutable double pots;                               ///< how many pot in this job?

      mutable bool fIsConfigLoaded = false;

      mutable bool fUsingDk2nu = true;

    }; // class NHLFluxCreator
      
  } // namespace NHL
} // namespace genie
