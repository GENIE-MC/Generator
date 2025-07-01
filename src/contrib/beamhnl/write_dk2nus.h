//____________________________________________________________________________/*
/*!
  ! This is a script to generate flat root trees from dk2nu tuples.
  ! It copies the dk2nu structure but does away with members of bsim
  ! so that GENIE HNL simulation can avoid having dk2nu as a compile-time package

\author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
         University of Oxford

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
 */
//____________________________________________________________________________

#ifndef write_dk2nus_h
#define write_dk2nus_h

#include <string>
#include <iostream>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <algorithm>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TChain.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TLorentzVector.h"

#include "FluxDrivers/GNuMIFlux.h"
#include "dk2nu/genie/GDk2NuFlux.h"

#include "tree/dk2nu.h"
#include "tree/calcLocationWeights.h"
#include "tree/dkmeta.h"

TRandom3 generator(0);
const Int_t setID = 0;

//==============================================================================
// Control
//==============================================================================
const int maxArray = 30;
const int maxC     = 100;

//==============================================================================
// Function Declaration
//==============================================================================
void   InitialiseMetaBranches(TTree * meta);
void   InitialiseTreeBranches(TTree * tree);
void   FillMetaBranches(bsim::DkMeta * dkmeta);
void   FillTreeBranches(bsim::Dk2Nu * dk2nu);
void   LoopEntries(TChain* cflux, TChain* dflux, bool grid, bool debug);
void   RootifyChar( std::string rfch, char fdch[maxC] );
void   LoadDetectorPosition( bool grid, genie::GFluxI * gfluxdriver );

//==============================================================================
// Directories, locations
//==============================================================================
const std::string USER            = std::getenv("USER") != NULL ? string(std::getenv("USER")) : "user";
const std::string OUTDIR          = std::getenv("OUTDIR") != NULL ? string(std::getenv("OUTDIR")) : string(std::getenv("PWD"));
const std::string SAMPLEDK2NU     = std::getenv("INDIR") != NULL ? string(std::getenv("INDIR"))+"/sample_dk2nu.root" : string(std::getenv("PWD"))+"/sample_dk2nu.root";

const std::string INDIR_GRID      = "";
const std::string OUTDIR_GRID     = "";

const std::string INDIR_DEBUG     = "./DEBUG";

const std::string GDK2NU_PSET = "MINERVA-v10r8";

//==============================================================================
// Input files: playlist location + length
//==============================================================================
const int m_maxNFiles             = 100;

//==============================================================================
// Branch names
// For more documentation about the meanings of these variables, visit
// http://cdcvs.fnal.gov/redmine/projects/lbne-beamsim/wiki/Dk2nu_ntuples
//==============================================================================

// --- meta
TTree * dkMeta = 0;

int         mArSize     = 0;               // Size of location arrays
int         mJob        = 0;               // Simulation job number
double      mPots       = 0.0;             // Number of POT simulated
char        mBeamsim[maxC];                // Describes the simulation that generated the file
char        mPhysics[maxC];                // Describes the geant version and physics list
char        mPhyscuts[maxC];               // Describes geant4 tracking cuts
char        mTgtcfg[maxC];                 // Describes target configuation
char        mHorncfg[maxC];                // Describes horn configuration
char        mDkvolcfg[maxC];               // Describes decay pipe configuration
double      mBeam0x     = 0.0;             // Initial x position of primary p in cm
double      mBeam0y     = 0.0;             // Initial y position of primary p in cm
double      mBeam0z     = 0.0;             // Initial z position of primary p in cm
double      mBeamhwidth = 0.0;             // Beam horizontal radius (rms) in cm
double      mBeamvwidth = 0.0;             // Beam vertical radius (rms) in cm
double      mBeamdxdz   = 0.0;             // Initial p dx/dz
double      mBeamdydz   = 0.0;             // Initial p dy/dz
double      mLocationDotX[maxArray];       // x position of locations in cm
double      mLocationDotY[maxArray];       // y position of locations in cm
double      mLocationDotZ[maxArray];       // z position of locations in cm
char      mLocationDotName[maxC*maxArray]; // Location names
//std::vector<std::string> mLocationDotName; // works but is SLOW compared to char hack
//std::vector<char>  mVintnames;             // Names associated with Vint
//std::vector<char>  mVdblnames;             // Names associated with Vdbl

// --- dk2nu
TTree * dkTree = 0;

int         dArSize           = 0;         // Size of location arrays
int         dAnArSize         = 0;         // Size of ancestor arrays
int         dTrArSize         = 0;         // Size of traj     arrays
// - - -
int         dJob              = 0;         // Simulation job number
double      dPotnum           = 0.0;       // N POT for this v (0 job-beginning, total POT job-end)
double      dPpvx             = 0.0;       // x component of v parent production vertex in cm
double      dPpvy             = 0.0;       // y component of v parent production vertex in cm
double      dPpvz             = 0.0;       // z component of v parent production vertex in cm
// - - -
int         dDecayDotNorig    = 0;         // v origin (1,2,3 = tgt/baffle, muon decay, other)
int         dDecayDotNdecay   = 0;         // Decay code of decay that produced v
int         dDecayDotNtype    = 0;         // GEANT particle code of v
double      dDecayDotVx       = 0.0;       // x component of v vertex position in cm
double      dDecayDotVy       = 0.0;       // y component of v vertex position in cm
double      dDecayDotVz       = 0.0;       // z component of v vertex position in cm
double      dDecayDotPdpx     = 0.0;       // x component of final parent momentum in GeV
double      dDecayDotPdpy     = 0.0;       // y component of final parent momentum in GeV
double      dDecayDotPdpz     = 0.0;       // z component of final parent momentum in GeV
double      dDecayDotPpdxdz   = 0.0;       // px/pz of parent at parent production point
double      dDecayDotPpdydz   = 0.0;       // py/pz of parent at parent production point
double      dDecayDotPppz     = 0.0;       //    pz of parent at parent production point in GeV
double      dDecayDotPpenergy = 0.0;       //     E of parent at parent production point in GeV
int         dDecayDotPpmedium = 0;         // empty branch
int         dDecayDotPtype    = 0;         // GEANT particle code of parent
double      dDecayDotMuparpx  = 0.0;       // (parent == mu) ? grandparent px in GeV : -99999
double      dDecayDotMuparpy  = 0.0;       // (parent == mu) ? grandparent py in GeV : -99999
double      dDecayDotMuparpz  = 0.0;       // (parent == mu) ? grandparent pz in GeV : -99999
double      dDecayDotMupare   = 0.0;       // (parent == mu) ? grandparent  E in GeV : -99999
double      dDecayDotNecm     = 0.0;       // v E in parent rest frame in GeV
double      dDecayDotNimpwt   = 0.0;       // Importance weight
// - - -
double      dNurayDotPx[maxArray];         // v px in GeV for each location in meta
double      dNurayDotPy[maxArray];         // v py in GeV for each location in meta
double      dNurayDotPz[maxArray];         // v pz in GeV for each location in meta
double      dNurayDotE[maxArray];          // v  E in GeV for each location in meta
double      dNurayDotWgt[maxArray];        // weights to make v flux spectra for each location in meta
// - - -
int         dAncestorDotPdg[maxArray];     // PDG code of ancestor
double      dAncestorDotStartx[maxArray];  // x component of ancestor start position in cm
double      dAncestorDotStarty[maxArray];  // y component of ancestor start position in cm
double      dAncestorDotStartz[maxArray];  // z component of ancestor start position in cm
double      dAncestorDotStartt[maxArray];  // t component of ancestor start position in (ns ?)
double      dAncestorDotStartpx[maxArray]; // ancestor initial px in GeV
double      dAncestorDotStartpy[maxArray]; // ancestor initial py in GeV
double      dAncestorDotStartpz[maxArray]; // ancestor initial pz in GeV
double      dAncestorDotStoppx[maxArray];  // ancestor final   px in GeV 
double      dAncestorDotStoppy[maxArray];  // ancestor final   py in GeV
double      dAncestorDotStoppz[maxArray];  // ancestor final   pz in GeV
double      dAncestorDotPolx[maxArray];    // empty branch
double      dAncestorDotPoly[maxArray];    // empty branch
double      dAncestorDotPolz[maxArray];    // empty branch
double      dAncestorDotPprodpx[maxArray]; // parent px prior to secondaries production (meaning?)
double      dAncestorDotPprodpy[maxArray]; // parent py prior to secondaries production (meaning?)
double      dAncestorDotPprodpz[maxArray]; // parent pz prior to secondaries production (meaning?)
int         dAncestorDotNucleus[maxArray]; // PDG code of nucleus where the interaction happened
char   dAncestorDotProc[maxArray*maxC];    // Describes processes that created each ancestor
char   dAncestorDotIvol[maxArray*maxC];    // Describes volume   where each ancestor was created
char   dAncestorDotImat[maxArray*maxC];    // Describes material where each ancestor was created
// - - -
double      dTgtexitDotTvx    = 0.0;       // x position of parent target exit in cm
double      dTgtexitDotTvy    = 0.0;       // y position of parent target exit in cm
double      dTgtexitDotTvz    = 0.0;       // z position of parent target exit in cm
double      dTgtexitDotTpx    = 0.0;       // parent px at target exit in GeV
double      dTgtexitDotTpy    = 0.0;       // parent py at target exit in GeV
double      dTgtexitDotTpz    = 0.0;       // parent pz at target exit in GeV
int         dTgtexitDotTptype = 0;         // GEANT particle code of ancestor that exited target
int         dTgtexitDotTgen   = 0;         // Generation number of v
// - - -
double      dTrajDotTrkx[maxArray];        // ?
double      dTrajDotTrky[maxArray];        // ?
double      dTrajDotTrkz[maxArray];        // ?
double      dTrajDotTrkpx[maxArray];       // ?
double      dTrajDotTrkpy[maxArray];       // ?
double      dTrajDotTrkpz[maxArray];       // ?
// - - -
//int         dFlagbits         = 0;         // extra user container
//std::vector<int>        dVint;             // ppfx specific
//std::vector<double>     dVdbl;             // ppfx specific

void LoadDetectorPosition(bool grid, genie::GFluxI* gfluxdriver)
{
  genie::flux::GFluxFileConfigI* gffconfig = 
    dynamic_cast<genie::flux::GFluxFileConfigI*>(gfluxdriver); 

  // change FLUXFILE_GRID to your favourite grid location
  std::string FLUXFILE_GRID = "";
  std::string sample_rootfile = grid ? FLUXFILE_GRID : SAMPLEDK2NU;
  sample_rootfile = gSystem->ExpandPathName(sample_rootfile.c_str());

  // this is generalized out in GFluxFileConfigI ...
  gffconfig->LoadBeamSimData(sample_rootfile, GDK2NU_PSET);

  std::cout << "\nDetector position loaded." << std::endl;
}

#endif //ifndef write_dk2nus_h
