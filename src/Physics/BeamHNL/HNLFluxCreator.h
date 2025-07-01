//----------------------------------------------------------------------------
/*!

  This is a module for GENIE to read in hadron flux ntuples and construct HNL fluxes
  on the fly. 
  
  Core loop: + Open flux entry
             + Get ancestry information
	     + Assume decay to HNL (discard SM nu info, is unneeded)
	     + Calculate HNL production mode based on parameter space read from config
	     + Calculate kinematics of HNL
	     + Return HNL as SimpleHNL object.

\class      genie::hnl::FluxCreator

\brief      Calculates HNL production kinematics & production vertex.
            Is a concrete implementation of the FluxRecordVisitorI interface

\author     John Plows <komninos-john.plows@physics.ox.ac.uk>

\created    April 25th, 2022

\cpright    Copyright (c) 2003-2025, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org

 */
//----------------------------------------------------------------------------

#ifndef _HNL_FLUXCREATOR_H_
#define _HNL_FLUXCREATOR_H_

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

#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>

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

#include "Physics/BeamHNL/HNLFluxRecordVisitorI.h"

#include "Physics/BeamHNL/HNLBRCalculator.h"
#include "Physics/BeamHNL/HNLDecayUtils.h"
#include "Physics/BeamHNL/HNLEnums.h"
#include "Physics/BeamHNL/HNLFluxContainer.h"
#include "Physics/BeamHNL/HNLKinUtils.h"
#include "Physics/BeamHNL/SimpleHNL.h"

const double kRDET = 1.0; // calculate fluxes per m^2

namespace genie{

  namespace hnl{
    
    class SimpleHNL;
    class FluxContainer;
    
    class FluxCreator : public FluxRecordVisitorI {

    public:

      FluxCreator();
      FluxCreator(string name);
      FluxCreator(string name, string config);
      ~FluxCreator();

      //-- implement the FluxRecordVisitorI interface
      void ProcessEventRecord(GHepRecord * event_rec) const;

      // overload the Algorithm::Configure() methods to load private data
      // members from configuration options
      void Configure(const Registry & config);
      void Configure(string config);
      
      // set the input path for a flux
      void SetInputFluxPath( std::string finpath ) const;
      // set path to geometry file
      void SetGeomFile( string geomfile ) const;
      // get N(flux input entries)
      int GetNFluxEntries() const;
      // set first entry for read-in from chain
      void SetFirstFluxEntry( int iFirst ) const;

      // get flux info
      FluxContainer RetrieveFluxInfo() const;

      // return information about frames
      std::vector< double > GetB2UTranslation() const { return fB2UTranslation; }
      std::vector< double > GetB2URotation() const { return fB2URotation; }
      std::vector< double > GetDetOffset() const { return fDetOffset; }
      std::vector< double > GetDetRotation() const { return fDetRotation; }

    private:

      void LoadConfig(void);

      // if using root geom, let this module know
      void SetUsingRootGeom( bool IsUsingRootGeom ) const;
      void ImportBoundingBox( TGeoBBox * box ) const;

      void SetCurrentEntry( int iCurr ) const;

      // workhorse methods
      FluxContainer MakeTupleFluxEntry( int iEntry, std::string finpath ) const;
      void FillNonsense( int iEntry, genie::hnl::FluxContainer & gnmf ) const;

      // init
      void OpenFluxInput( std::string finpath ) const;
      void InitialiseTree() const;
      void InitialiseMeta() const;

      // returns HNL 4-momentum from random decay in same frame as p4par
      TLorentzVector HNLEnergy( genie::hnl::HNLProd_t hnldm, TLorentzVector p4par ) const;
      // gets random point in BBox and returns separation to it in BEAM FRAME
      TVector3 PointToRandomPointInBBox( ) const;

      void ReadBRs() const;
      std::map< genie::hnl::HNLProd_t, double > GetProductionProbs( int parPDG ) const;
      
      // Obtain detector dimensions + position
      // RETHERE: BBox isn't good enough! But roll with it for now
      void MakeBBox() const;
      TVector3 ApplyUserRotation( TVector3 vec, bool doBackwards = false ) const;
      TVector3 ApplyUserRotation( TVector3 vec, TVector3 oriVec, std::vector<double> rotVec, bool doBackwards = false ) const;
      
      // calculate detector acceptance (== solid angle of proj of det onto unit-radius sphere / (4pi))
      // NOTE THIS IS A LAB FRAME (==GEOMETRICAL) ACCEPTANCE!!!!
      // detO == detector BBox centre wrt HNL prod vertex, L{x,y,z} BBox length on each axis. Both [m]
      double CalculateDetectorAcceptanceSAA( TVector3 detO ) const;
      // collimation effect calc, returns HNL_acc / geom_acc
      double CalculateAcceptanceCorrection( TLorentzVector p4par, TLorentzVector p4HNL, double SMECM, double zm, double zp ) const;
      static double labangle( double * x, double * par ); // function formula for correction
      // get minimum and maximum deviation from parent momentum to hit detector, [deg]
      double GetAngDeviation( TLorentzVector p4par, TVector3 detO, bool seekingMax ) const;
      void GetAngDeviation( TLorentzVector p4par, TVector3 detO, double &zm, double &zp ) const;
      // returns 1.0 / (area of flux calc)
      double CalculateAreaNormalisation();

      // utility function -- is a copy of TGeoChecker::CheckPoint() but doesn't output to cout
      std::string CheckGeomPoint( Double_t x, Double_t y, Double_t z ) const;

      // current path to keep track of what is loaded
      mutable std::string fCurrPath = ""; mutable bool fPathLoaded = false;
      // and which entry we're on
      mutable int iCurrEntry = 0;
      // which one was first?
      mutable int fFirstEntry = 0;
      // out of how many?
      mutable int fNEntries = 0;
      
      // maps to keep P( production )
      mutable std::map< genie::hnl::HNLProd_t, double > dynamicScores; // map in use
      mutable std::map< genie::hnl::HNLProd_t, double > dynamicScores_pion;
      mutable std::map< genie::hnl::HNLProd_t, double > dynamicScores_kaon;
      mutable std::map< genie::hnl::HNLProd_t, double > dynamicScores_muon;
      mutable std::map< genie::hnl::HNLProd_t, double > dynamicScores_neuk;
      
      mutable double BR_pi2mu, BR_pi2e, BR_K2mu, BR_K2e, BR_K3mu, BR_K3e, BR_K03mu, BR_K03e;

      mutable bool isParentOnAxis = true;
      mutable TGeoVolume * fTopVol = 0;
      mutable string fGeomFile = "";
      mutable bool fIsUsingRootGeom = false;

      mutable TChain * ctree = 0, * cmeta = 0;

      mutable double fMass; // HNL mass, GeV
      mutable std::vector< double > fU4l2s; // couplings
      mutable bool fIsMajorana;
      //mutable int fType; // for hist fluxes. 0 ==> only particle, 1 ==> only anti, 2 ==> mix of both
      //mutable double fMinE = 0.0, fMaxE = 100.0, fAngDev = 0.0; // for hist fluxes
      
      mutable double fLx, fLy, fLz;
      mutable double fLxR, fLyR, fLzR; // BBox side [m]

      mutable std::vector< double > fB2UTranslation, fB2URotation;
      mutable std::vector< double > fDetRotation; // rotation of detector wrt tgt hall
      mutable std::vector< double > fDetOffset; // offset of det centre wrt geom file origin
      mutable double fCx, fCy, fCz;   // BBox centre wrt HNL prod [m]
      mutable double fAx1, fAz, fAx2; // Euler angles, extrinsic x-z-x. Tgt-hall to beam
      mutable double fBx1, fBz, fBx2; // Tgt-hall to detector frame

      mutable double fDx, fDy, fDz; //HNL production point [m]

      //std::vector< double > trVec, roVec;

      mutable bool doPol = true, fixPol = false;
      mutable double fLPx, fLPy, fLPz; // direction of co-produced lepton == polarisation vector
      mutable double fLPE;
      mutable std::vector< double > fFixedPolarisation;

      mutable int fLepPdg; // pdg code of co-produced lepton
      mutable int fNuPdg; // pdg code of SM neutrino from same decay type
      mutable double parentMass, parentMomentum, parentEnergy; // GeV
      mutable double fECM, fSMECM; // GeV
      mutable double fZm, fZp; // deg

      mutable int fProdChan, fNuProdChan;

      mutable TVector3 fTargetPoint, fTargetPointUser;

      static const int maxArray = 30, maxC = 100;

      // tree variables. Add as per necessary.
      mutable double potnum;                             ///< N POT for this SM-v
      mutable int    decay_ptype;                        ///< PDG code of parent
      mutable double decay_vx, decay_vy, decay_vz;       ///< coordinates of prod vtx [cm]
      mutable double decay_pdpx, decay_pdpy, decay_pdpz; ///< final parent momentum [GeV]
      mutable double decay_necm;                         ///< SM v CM energy [GeV]
      mutable double decay_nimpwt;                       ///< Importance weight from beamsim

      mutable int arSize, anArSize, trArSize;
      mutable int djob;
      mutable double ppvx, ppvy, ppvz;
      mutable int decay_norig, decay_ndecay, decay_ntype;
      mutable double decay_ppdxdz, decay_ppdydz, decay_pppz, decay_ppenergy;
      mutable int decay_ppmedium;
      mutable double decay_muparpx, decay_muparpy, decay_muparpz, decay_mupare;

      mutable double nuray_px[maxArray], nuray_py[maxArray], nuray_pz[maxArray], nuray_E[maxArray], nuray_wgt[maxArray];
      
      mutable int ancestor_pdg[maxArray];
      mutable double ancestor_startx[maxArray], ancestor_starty[maxArray], ancestor_startz[maxArray], ancestor_startt[maxArray];
      mutable double ancestor_startpx[maxArray], ancestor_startpy[maxArray], ancestor_startpz[maxArray];
      mutable double ancestor_stoppx[maxArray], ancestor_stoppy[maxArray], ancestor_stoppz[maxArray];
      mutable double ancestor_polx[maxArray], ancestor_poly[maxArray], ancestor_polz[maxArray];
      mutable double ancestor_pprodpx[maxArray], ancestor_pprodpy[maxArray], ancestor_pprodpz[maxArray];
      mutable int ancestor_nucleus[maxArray];
      mutable char ancestor_proc[maxArray*maxC], ancestor_ivol[maxArray*maxC], ancestor_imat[maxArray*maxC];

      mutable double tgtexit_tvx, tgtexit_tvy, tgtexit_tvz;
      mutable double tgtexit_tpx, tgtexit_tpy, tgtexit_tpz;
      mutable int tgtexit_tptype, tgtexit_tgen;

      mutable double traj_trkx[maxArray], traj_trky[maxArray], traj_trkz[maxArray];
      mutable double traj_trkpx[maxArray], traj_trkpy[maxArray], traj_trkpz[maxArray];
      
      // meta variables. Add as necessary
      mutable int    job;                                ///< beamsim MC job number
      mutable double pots;                               ///< how many pot in this job?

      mutable int mArSize;
      mutable char beamsim[maxC], physics[maxC], physcuts[maxC];
      mutable char tgtcfg[maxC], horncfg[maxC], dkvolcfg[maxC];
      mutable double beam0x, beam0y, beam0z;
      mutable double beamhwidth, beamvwidth;
      mutable double beamdxdz, beamdydz;
      mutable double location_x[maxArray], location_y[maxArray], location_z[maxArray];
      mutable char location_name[maxArray*maxC];

      mutable genie::hnl::FluxContainer fGnmf;

      mutable double POTScaleWeight;
      mutable std::vector<double> fScales;
      
      mutable bool fDoingOldFluxCalc = false;
      mutable bool fRerollPoints = false;
      mutable double fRadius; // m
      mutable bool fSupplyingBEAM = false;
      mutable bool fIsConfigLoaded = false;

      mutable bool fUsingDk2nu = true;
      mutable string fFinPath, fProdHist;
      mutable TH1D * fSpectrum = 0, * fIntegrals = 0;

    }; // class FluxCreator
      
  } // namespace hnl
} // namespace genie

#endif // #ifndef _HNL_FLUXCREATOR_H_
