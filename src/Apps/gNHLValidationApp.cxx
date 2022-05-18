//________________________________________________________________________________________
/*!

\program gevald_nhl

\brief   Validation app for NHL generation.

         *** Synopsis :

         gevald_nhl [-h]
                    [-r run#]
                     -n n_of_events
		     -f path/to/flux/files
                    [-E nhl_energy]
		    [-g geometry (ROOT file)]
		    [-L geometry_length_units]
		    [-A geometry_angle_units]
		    [-T geometry_time_units]
		    [-o output_event_file_prefix]
		    [--seed random_number_seed]

         *** Options :

           [] Denotes an optional argument

           -h
              Prints out the gevald_nhl syntax and exits.
           -r
              Specifies the MC run number [default: 1000].
           -n
              Specifies how many events to generate.
           -f
              Input NHL flux.
	   -E
	      NHL energy for monoenergetic NHL
           -g
              Input detector geometry.
              If a geometry is specified, NHL decay vertices will be distributed
              in the desired detector volume.
              Using this argument, you can pass a ROOT file containing your
              detector geometry description.
	   -L
	      Input geometry length units, eg 'm', 'cm', 'mm', ...
              [default: 'mm']
	   -A
	      Input geometry angle units: 'deg', 'rad', or 'mrad'
	      [default: 'rad']
	   -T
	      Input geometry time units, eg 's', 'us', 'ns', ...
	      [default: 'ns']
	   -o
              Sets the prefix of the output event file.
              The output filename is built as:
              [prefix].[run_number].[event_tree_format].[file_format]
              The default output filename is:
              gntp.[run_number].ghep.root
              This cmd line arguments lets you override 'gntp'
           --seed
              Random number seed.

\author  Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
         University of Liverpool & STFC Rutherford Appleton Laboratory

	 John Plows <komninos-john.plows \at cern.ch>
	 University of Oxford

\created May 17, 2022

\cpright Copyright (c) 2003-2022, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>

#include <TSystem.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpWriter.h"

#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"
#include "Physics/NeutralHeavyLepton/NHLDecayUtils.h"
#include "Physics/NeutralHeavyLepton/NHLDecayVolume.h"
#include "Physics/NeutralHeavyLepton/NHLFluxCreator.h"
#include "Physics/NeutralHeavyLepton/NHLFluxReader.h"
#include "Physics/NeutralHeavyLepton/NHLPrimaryVtxGenerator.h"
#include "Physics/NeutralHeavyLepton/NHLProductionMode.h"
#include "Physics/NeutralHeavyLepton/SimpleNHL.h"

#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/CmdLnArgParser.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;
using namespace genie::NHL;

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX__
#include "Tools/Flux/GCylindTH1Flux.h"
#include "Tools/Flux/GNuMIFlux.h"
#include <TH1.h>
#endif // #ifdef __GENIE_FLUX_DRIVERS_ENABLED__

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#define __CAN_USE_ROOT_GEOM__
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

typedef enum t_NHLValidation {
  
  kValNone            = 0,
  kValFluxFromDk2nu   = 1,
  kValFluxFromHists   = 2,
  kValProdEvtRates    = 3,
  kValProdKinematics  = 4,
  kValDecayEvtRates   = 5,
  kValDecayKinematics = 6,
  kValEntryExit       = 7,
  kValVtxPlacement    = 8,
  kValParentGun       = 9,
  kValFullSim         = 10
  
} NHLValidation_t;

// function prototypes
void  GetCommandLineArgs (int argc, char ** argv);
void  ReadInConfig       (void);
void  PrintSyntax        (void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
int   TestFluxFromDk2nu  (void);
int   InitialiseTupleFlux(string finpath);
void  MakeNHLFromTuple   (int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf, string finpath, int run);
#endif

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomAUnits   = "rad";   // default geometry angle units
string          kDefOptGeomTUnits   = "ns";    // default geometry time units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units

NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format
string          kDefOptEvFilePrefix = "gntp";
string          kDefOptFluxFilePath = "./input-flux.root";

string          kDefOptSName   = "genie::EventGenerator";
string          kDefOptSConfig = "NeutralHeavyLepton";

//
Long_t           gOptRunNu        = 1000;                // run number
int              gOptNev          = 10;                  // number of events to generate

NHLValidation_t  gOptValidationMode = kValNone;          // which test to run

// ---- this section set in config

double           gCfgMassNHL      = -1;                  // NHL mass
double           gCfgECoupling    = -1;                  // |U_e4|^2
double           gCfgMCoupling    = -1;                  // |U_m4|^2
double           gCfgTCoupling    = -1;                  // |U_t4|^2

NHLProd_t        gCfgProdMode     = kNHLProdNull;        // NHL production mode
NHLDecayMode_t   gCfgDecayMode    = kNHLDcyNull;         // NHL decay mode
std::vector< NHLDecayMode_t > gCfgIntChannels;           // decays to un-inhibit

double           gCfgUserOx       = -1;                  // user origin x in beam coords
double           gCfgUserOy       = -1;                  // user origin y in beam coords
double           gCfgUserOz       = -1;                  // user origin z in beam coords
double           gCfgUserAx1      = -1;                  // left-most extrinsic Euler angle (x)
double           gCfgUserAz       = -1;                  // middle extrinsic Euler angle (z)
double           gCfgUserAx2      = -1;                  // right-most extrinsic Euler angle (x)

double           gCfgBoxLx        = -1;                  // detector length on user x
double           gCfgBoxLy        = -1;                  // detector length on user y
double           gCfgBoxLz        = -1;                  // detector length on user z

double           gCfgParentEnergy = -1;                  // parent energy in GeV
double           gCfgParentTheta  = -1;                  // wrt beam axis
double           gCfgParentPhi    = -1;                  // 0 == x, pi/2 == y
double           gCfgParentOx     = -1;                  // user x of parent origin
double           gCfgParentOy     = -1;                  // user y of parent origin
double           gCfgParentOz     = -1;                  // user z of parent origin

// ---- end config section

double           gOptEnergyNHL    = -1;                  // NHL energy in GeV

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
string           gOptFluxFilePath = kDefOptFluxFilePath; // where flux files live
bool             gOptIsUsingDk2nu = false;               // using flat dk2nu files?
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
bool             gOptIsMonoEnFlux = true;                // do we have monoenergetic flux?

string           gOptEvFilePrefix = kDefOptEvFilePrefix; // event file prefix
bool             gOptUsingRootGeom = false;              // using root geom?
string           gOptRootGeom;                           // input ROOT file with realistic detector geometry

string           lunits, aunits, tunits;
double           gOptGeomLUnits = 0;                     // input geometry length units
double           gOptGeomAUnits = 0;                     // input geometry angle units
double           gOptGeomTUnits = 0;                     // input geometry time units
#ifdef __CAN_USE_ROOT_GEOM__
TGeoManager *    gOptRootGeoManager = 0;                 // the workhorse geometry manager
TGeoVolume  *    gOptRootGeoVolume  = 0;
#endif // #ifdef __CAN_USE_ROOT_GEOM__
string           gOptRootGeomTopVol = "";                // input geometry top event generation volume

long int         gOptRanSeed = -1;                       // random number seed

NHLPrimaryVtxGenerator * nhlgen = 0;

//_________________________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  // Get the validation configuration
  ReadInConfig();

  // Init messenger and random number seed
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);

  LOG( "gevald_nhl", pFATAL )
    << "\nI am a work in progress. Hello world!"
    << "\n."
    << "\n. ."
    << "\n. . .";

  switch( gOptValidationMode ){
    
  case kValFluxFromDk2nu: return TestFluxFromDk2nu(); break;
  default: LOG( "gevald_nhl", pFATAL ) << "I didn't recognise this mode. Goodbye world!"; break;

  }

  return 0;
}
//_________________________________________________________________________________________
//............................................................................
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
int TestFluxFromDk2nu()
{
  assert( gOptIsUsingDk2nu );

  string foutName("test_flux_dk2nu.root");
  
  LOG( "gevald_nhl", pINFO )
    << "\n\nTesting flux prediction from dk2nu input files."
    << "\nWill produce 1 ROOT file (" << foutName << ") with:"
    << "\n--> Energy spectra for the NHL (total and broken down by parent)"
    << "\n--> Production vertex locations"
    << "\n--> Counters for each production mode"
    << "\n--> Spectrum of acceptance correction as function of parent boost factor"
    << "\n--> Boost factor spectrum of parents broken down by type";

  int maxFluxEntries = InitialiseTupleFlux( gOptFluxFilePath );
  if( gOptNev > maxFluxEntries ){
    LOG( "gevald_nhl", pWARN )
      << "You have asked for " << gOptNev << " events, but only provided "
      << maxFluxEntries << " flux entries. Truncating events to " << maxFluxEntries << ".";
    gOptNev = maxFluxEntries;
  }

  TFile * fout = TFile::Open( foutName.c_str(), "RECREATE" );
  TH1D hEAll, hEPion, hEKaon, hEMuon, hENeuk;
  TH3D hProdVtxPos;
  TH1D hCounters;
  TH2D hAccCorrVsBoostBeta;
  TH1D hBAll, hBPion, hBKaon, hBMuon, hBNeuk;

  hEAll  = TH1D( "hEAll",  "NHL energy - all parents", 100, 0., 100. );
  hEPion = TH1D( "hEPion", "NHL energy - pion parent", 100, 0., 100. );
  hEKaon = TH1D( "hEKaon", "NHL energy - kaon parent", 100, 0., 100. );
  hEMuon = TH1D( "hEMuon", "NHL energy - muon parent", 100, 0., 100. );
  hENeuk = TH1D( "hENeuk", "NHL energy - neuk parent", 100, 0., 100. );
      
  hProdVtxPos = TH3D( "hProdVtxPos", "NHL production vertex (user coordinates, cm)",
		      200, -100., 100., 200, -100., 100., 1100, -110000., 0.);
  
  hCounters = TH1D( "hCounters", "NHL production channel counters", 11, 0, 10 );
  
  hAccCorrVsBoostBeta = TH2D( "hAccCorrVsBoostBeta", "Acceptance correction vs boost beta",
			      100, 0., 1., 500, 0., 50.);
  
  hBAll  = TH1D( "hBAll",  "Boost beta - all parents", 100, 0., 1. );
  hBPion = TH1D( "hBPion", "Boost beta - pion parent", 100, 0., 1. );
  hBKaon = TH1D( "hBKaon", "Boost beta - kaon parent", 100, 0., 1. );
  hBMuon = TH1D( "hBMuon", "Boost beta - muon parent", 100, 0., 1. );
  hBNeuk = TH1D( "hBNeuk", "Boost beta - neuk parent", 100, 0., 1. );

  int parPDG;
  TLorentzVector p4NHL;
  TLorentzVector x4NHL;
  int nPion2Muon = 0, nPion2Electron = 0, nKaon2Muon = 0,
    nKaon2Electron = 0, nKaon3Muon = 0, nKaon3Electron = 0,
    nNeuk3Muon = 0, nNeuk3Electron = 0, nMuon3Numu = 0,
    nMuon3Nue = 0, nMuon3Nutau = 0;
  double betaMag;
  double accCorr;
  double nimpwt; // hadroproduction importance weight

  flux::GNuMIFluxPassThroughInfo * gnmf = new flux::GNuMIFluxPassThroughInfo();

  int ievent = 0;
  while(true)
    {
      if( gOptNev >= 10000 ){
	if( ievent % (gOptNev / 1000) == 0 ){
	  int irat = ievent / (gOptNev / 1000);
	  std::cerr << Form("%2.2f", 0.1 * irat) << " % ( " << ievent << " / "
		    << gOptNev << " ) \r" << std::flush;
	}
      }
      
      if( ievent == gOptNev ) break;

      MakeNHLFromTuple( ievent, gnmf, gOptFluxFilePath, gOptRunNu );
      
      // reject nonsense
      if( gnmf->necm < 0 ){
	
	LOG( "gevald_nhl", pDEBUG )
	  << "Skipping nonsense for event " << ievent << " (was this parent too light?)";

      } else {
	// now to make stuff from this... i.e. fill histos
	
	// first get the channel
	int iChannel = gnmf->ndecay - 30; int typeMod = 1;
	if( iChannel % 2 == 0 ){ iChannel--; typeMod = -1; }
	NHLGNuMIProd_t gChannel = ( NHLGNuMIProd_t ) iChannel;
	NHLProd_t pChannel;

	switch( gChannel ){
	case kNHLGProdNeuk3Electron:
	  parPDG = kPdgK0L; nNeuk3Electron++; pChannel = kNHLProdNeuk3Electron; break;
	case kNHLGProdNeuk3Muon:
	  parPDG = kPdgK0L; nNeuk3Muon++; pChannel = kNHLProdNeuk3Muon; break;
	case kNHLGProdKaon2Electron:
	  parPDG = kPdgKP; nKaon2Electron++; pChannel = kNHLProdKaon2Electron; break;
	case kNHLGProdKaon2Muon:
	  parPDG = kPdgKP; nKaon2Muon++; pChannel = kNHLProdKaon2Muon; break;
	case kNHLGProdKaon3Electron:
	  parPDG = kPdgKP; nKaon3Electron++; pChannel = kNHLProdKaon3Electron; break;
	case kNHLGProdKaon3Muon:
	  parPDG = kPdgKP; nKaon3Muon++; pChannel = kNHLProdKaon3Muon; break;
	case kNHLGProdMuon3Nue:
	  parPDG = kPdgMuon; nMuon3Nue++; pChannel = kNHLProdMuon3Nue; break;
	case kNHLGProdMuon3Numu:
	  parPDG = kPdgMuon; nMuon3Numu++; pChannel = kNHLProdMuon3Numu; break;
	case kNHLGProdMuon3Nutau:
	  parPDG = kPdgMuon; nMuon3Nutau++; pChannel = kNHLProdMuon3Nutau; break;
	case kNHLGProdPion2Electron:
	  parPDG = kPdgPiP; nPion2Electron++; pChannel = kNHLProdPion2Electron; break;
	case kNHLGProdPion2Muon:
	  parPDG = kPdgPiP; nPion2Muon++; pChannel = kNHLProdPion2Muon; break;
	default:
	  LOG( "gevald_nhl", pERROR )
	    << "Unknown decay mode " << iChannel << " at entry " << gnmf->evtno;
	  parPDG = 0; break;
	}
	
	if( parPDG == 0 ){ ievent++; continue; }
	
	double MPar = PDGLibrary::Instance()->Find( parPDG )->Mass();
	TVector3 p3par( gnmf->xpoint, gnmf->ypoint, gnmf->zpoint );
	double EPar = std::sqrt( p3par.Mag2() + MPar*MPar );
	TLorentzVector p4par( p3par.Px(), p3par.Py(), p3par.Pz(), EPar );
	
	TVector3 boost_beta = p4par.BoostVector();
	betaMag = boost_beta.Mag();
	
	p4NHL = gnmf->fgP4User;
	x4NHL = gnmf->fgX4User; // in cm

	double acceptance = gnmf->fgXYWgt; // full acceptance
	accCorr    = gnmf->nwtnear;   // just the correction
	nimpwt = gnmf->nimpwt;
	
	// fill the histos!
	hEAll.Fill( p4NHL.E(), acceptance * nimpwt );
	hProdVtxPos.Fill( x4NHL.X(), x4NHL.Y(), x4NHL.Z(), nimpwt );
	hAccCorrVsBoostBeta.Fill( betaMag, accCorr, nimpwt );
	hBAll.Fill( betaMag, nimpwt );
	
	switch( parPDG ){
	case kPdgPiP:
	  hEPion.Fill( p4NHL.E(), acceptance * nimpwt ); 
	  hBPion.Fill( betaMag, nimpwt );
	  break;
	case kPdgKP:
	  hEKaon.Fill( p4NHL.E(), acceptance * nimpwt ); 
	  hBKaon.Fill( betaMag, nimpwt );
	  break;
	case kPdgK0L:
	  hENeuk.Fill( p4NHL.E(), acceptance * nimpwt ); 
	  hBNeuk.Fill( betaMag, nimpwt );
	  break;
	case kPdgMuon:
	  hEMuon.Fill( p4NHL.E(), acceptance * nimpwt ); 
	  hBMuon.Fill( betaMag, nimpwt );
	  break;
	}
	
	LOG( "gevald_nhl", pDEBUG )
	  << " *** Output for event no " << ievent << "... ***"
	  << "\nparPDG = " << parPDG
	  << "\np4par = " << utils::print::P4AsString( &p4par ) << " [GeV] "
	  << "\nbetaVec = " << utils::print::Vec3AsString( &boost_beta )
	  << " : mag = " << betaMag
	  << "\np4NHL = " << utils::print::P4AsString( &p4NHL ) << " [GeV] "
	  << "\nx4NHL = " << utils::print::X4AsString( &x4NHL ) << " [cm]"
	  << "\nacceptance = " << acceptance << " : accCorr = " << accCorr
	  << "\nnimpwt = " << nimpwt
	  << "\nProduction channel = " << utils::nhl::ProdAsString( pChannel );

      } // if not nonsense

      ievent++;
    }

  // fill hCounters; each bin corresponds to the NHLProd_t with index iBin-1
  hCounters.SetBinContent( 1,  nPion2Muon );
  hCounters.SetBinContent( 2,  nPion2Electron );
  hCounters.SetBinContent( 3,  nKaon2Muon );
  hCounters.SetBinContent( 4,  nKaon2Electron );
  hCounters.SetBinContent( 5,  nKaon3Muon );
  hCounters.SetBinContent( 6,  nKaon3Electron );
  hCounters.SetBinContent( 7,  nNeuk3Muon );
  hCounters.SetBinContent( 8,  nNeuk3Electron );
  hCounters.SetBinContent( 9,  nMuon3Numu );
  hCounters.SetBinContent( 10, nMuon3Nue );
  hCounters.SetBinContent( 11, nMuon3Nutau );

  LOG( "gevald_nhl", pDEBUG )
    << "Showing production channel stats:"
    << "\n" << utils::nhl::ProdAsString( kNHLProdPion2Muon ) << ": " << nPion2Muon
    << "\n" << utils::nhl::ProdAsString( kNHLProdPion2Electron ) << ": " << nPion2Electron
    << "\n" << utils::nhl::ProdAsString( kNHLProdKaon2Muon ) << ": " << nKaon2Muon
    << "\n" << utils::nhl::ProdAsString( kNHLProdKaon2Electron ) << ": " << nKaon2Electron
    << "\n" << utils::nhl::ProdAsString( kNHLProdKaon3Muon ) << ": " << nKaon3Muon
    << "\n" << utils::nhl::ProdAsString( kNHLProdKaon3Electron ) << ": " << nKaon3Electron
    << "\n" << utils::nhl::ProdAsString( kNHLProdNeuk3Muon ) << ": " << nNeuk3Muon
    << "\n" << utils::nhl::ProdAsString( kNHLProdNeuk3Electron ) << ": " << nNeuk3Electron
    << "\n" << utils::nhl::ProdAsString( kNHLProdMuon3Numu ) << ": " << nMuon3Numu
    << "\n" << utils::nhl::ProdAsString( kNHLProdMuon3Nue ) << ": " << nMuon3Nue
    << "\n" << utils::nhl::ProdAsString( kNHLProdMuon3Nutau ) << ": " << nMuon3Nutau;

  fout->Write();
  fout->Close();

  return 0;
}
//_________________________________________________________________________________________
int InitialiseTupleFlux( std::string finpath )
{
  LOG( "gevald_nhl", pDEBUG )
    << "Opening input flux now from path " << finpath.c_str();

  NHLFluxCreator::OpenFluxInput( gOptFluxFilePath );
  assert( NHLFluxCreator::tree && NHLFluxCreator::meta && NHLFluxCreator::tree->GetEntries() > 0 );
  return NHLFluxCreator::tree->GetEntries();
}
//_________________________________________________________________________________________
void MakeNHLFromTuple( int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf, std::string finpath, int run )
{
  // This genereates a full NHL from the flux tuples
  // by interfacing with NHLFluxCreator

  NHLFluxCreator::MakeTupleFluxEntry( iEntry, gnmf, finpath, run );
  
  LOG( "gevald_nhl", pDEBUG ) << "MakeNHLFromTuple complete.";
}
//............................................................................
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_
//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevald_nhl", pINFO) << "Parsing command line arguments";

  // Common run options.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
    PrintSyntax();
    exit(0);
  }

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevald_nhl", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevald_nhl", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevald_nhl", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevald_nhl", pFATAL)
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  if( parser.OptionExists('M') ) {
    LOG("gevald_nhl", pDEBUG)
      << "Detecting mode. . .";
    gOptValidationMode = (NHLValidation_t) parser.ArgAsInt('M');
  } else {
    LOG("gevald_nhl", pFATAL)
      << "You must specify a validation mode.";
    PrintSyntax();
    exit(0);
  } // -M

  bool isMonoEnergeticFlux = true;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  if( parser.OptionExists('f') ) {
    LOG("gevald_nhl", pDEBUG)
      << "A flux has been offered. Searching this path: " << parser.ArgAsString('f');
    isMonoEnergeticFlux = false;
    gOptFluxFilePath = parser.ArgAsString('f');
    
    // check if this is dk2nu
    if( gOptFluxFilePath.find( "dk2nu" ) != string::npos ){
      gOptIsUsingDk2nu = true;
      LOG("gevald_nhl", pDEBUG)
	<< "dk2nu flux files detected. Will create flux spectrum dynamically.";
    }
  } else {
    // we need the 'E' option! Log it and pass below
    LOG("gevald_nhl", pINFO)
      << "No flux file offered. Assuming monoenergetic flux.";
  } //-f
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

  if( parser.OptionExists('L') ) {
    lunits = parser.ArgAsString('L');
    LOG("gevald_nhl", pDEBUG) << "Setting length units to " << lunits.c_str();
  } else {
    LOG("gevald_nhl", pDEBUG) << "Using default geometry length units";
    lunits = kDefOptGeomLUnits;
  } // -L
  gOptGeomLUnits = utils::units::UnitFromString(lunits);

  if( parser.OptionExists('A') ) {
    aunits = parser.ArgAsString('A');
    LOG("gevald_nhl", pDEBUG) << "Setting angle units to " << aunits.c_str();
  } else {
    LOG("gevald_nhl", pDEBUG) << "Using default angle length units";
    aunits = kDefOptGeomAUnits;
  } // -A
  gOptGeomAUnits = utils::units::UnitFromString(aunits);

  if( parser.OptionExists('T') ) {
    tunits = parser.ArgAsString('T');
    LOG("gevald_nhl", pDEBUG) << "Setting time units to " << tunits.c_str();
  } else {
    LOG("gevald_nhl", pDEBUG) << "Using default geometry time units";
    tunits = kDefOptGeomTUnits;
  } // -T
  gOptGeomTUnits = utils::units::UnitFromString(tunits);

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevald_nhl", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevald_nhl", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevald_nhl", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevald_nhl", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  //
  // >>> print the command line options
  //

  ostringstream gminfo;
  if (gOptUsingRootGeom) {
    gminfo << "Using ROOT geometry - file: " << gOptRootGeom
           << ", top volume: "
           << ((gOptRootGeomTopVol.size()==0) ? "<master volume>" : gOptRootGeomTopVol)
           << ", length  units: " << lunits;
           // << ", density units: " << dunits;
  }

  LOG("gevald_nhl", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevald_nhl job configuration");

  LOG("gevald_nhl", pNOTICE)
     << "\n @@ Run number    : " << gOptRunNu
     << "\n @@ Random seed   : " << gOptRanSeed
     << "\n @@ NHL mass      : " << gCfgMassNHL << " GeV"
     << "\n @@ Decay channel : " << utils::nhl::AsString(gCfgDecayMode)
     << "\n @@ Flux path     : " << gOptFluxFilePath
     << "\n @@ Geometry      : " << gminfo.str()
     << "\n @@ Statistics    : " << gOptNev << " events";
}
//_________________________________________________________________________________________
void ReadInConfig(void)
{
  LOG("gevald_nhl", pFATAL)
    << "Reading in validation configuration. . .";

  gCfgMassNHL   = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Mass" );
  gCfgECoupling = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Ue42" );
  gCfgMCoupling = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Um42" );
  gCfgTCoupling = utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Ut42" );

  gCfgProdMode  = (NHLProd_t) utils::nhl::GetCfgInt( "NHL", "Validation", "NHL-ProdMode"  );
  gCfgDecayMode = (NHLDecayMode_t) utils::nhl::GetCfgInt( "NHL", "Validation", "NHL-DecayMode" );
  
  gCfgIntChannels = {};
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-2B_mu_pi" ) )
    gCfgIntChannels.emplace_back( kNHLDcyPiMu );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-2B_e_pi" ) )
    gCfgIntChannels.emplace_back( kNHLDcyPiE );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-2B_nu_pi0" ) )
    gCfgIntChannels.emplace_back( kNHLDcyPi0Nu );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-3B_nu_nu_nu" ) )
    gCfgIntChannels.emplace_back( kNHLDcyNuNuNu );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-3B_nu_mu_mu" ) )
    gCfgIntChannels.emplace_back( kNHLDcyNuMuMu );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-3B_nu_e_e" ) )
    gCfgIntChannels.emplace_back( kNHLDcyNuEE );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-3B_nu_mu_e" ) )
    gCfgIntChannels.emplace_back( kNHLDcyNuMuE );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-3B_e_pi_pi0" ) )
    gCfgIntChannels.emplace_back( kNHLDcyPiPi0E );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-3B_mu_pi_pi0" ) )
    gCfgIntChannels.emplace_back( kNHLDcyPiPi0Mu );
  if( utils::nhl::GetCfgBool( "NHL", "InterestingChannels", "NHL-3B_nu_pi0_pi0" ) )
    gCfgIntChannels.emplace_back( kNHLDcyPi0Pi0Nu );

  gCfgUserOx    = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "DetCentreXInBeam" );
  gCfgUserOy    = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "DetCentreYInBeam" );
  gCfgUserOz    = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "DetCentreZInBeam" );

  gCfgUserAx1   = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "EulerExtrinsicX1" );
  gCfgUserAz    = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "EulerExtrinsicZ"  );
  gCfgUserAx2   = utils::nhl::GetCfgDouble( "NHL", "CoordinateXForm", "EulerExtrinsicX2" );

  gCfgBoxLx     = utils::nhl::GetCfgDouble( "NHL", "Validation", "BoxLx" );
  gCfgBoxLy     = utils::nhl::GetCfgDouble( "NHL", "Validation", "BoxLy" );
  gCfgBoxLz     = utils::nhl::GetCfgDouble( "NHL", "Validation", "BoxLz" );

  gCfgParentEnergy = utils::nhl::GetCfgDouble( "NHL", "Validation", "Parent-E"     );
  gCfgParentTheta  = utils::nhl::GetCfgDouble( "NHL", "Validation", "Parent-Theta" );
  gCfgParentPhi    = utils::nhl::GetCfgDouble( "NHL", "Validation", "Parent-Phi"   );

  gCfgParentOx  = utils::nhl::GetCfgDouble( "NHL", "Validation", "Parent-Ox" );
  gCfgParentOy  = utils::nhl::GetCfgDouble( "NHL", "Validation", "Parent-Oy" );
  gCfgParentOz  = utils::nhl::GetCfgDouble( "NHL", "Validation", "Parent-Oz" );

  // now transform the lengths and angles to the correct units
  gCfgUserOx   *= units::m / gOptGeomLUnits;
  gCfgUserOy   *= units::m / gOptGeomLUnits;
  gCfgUserOz   *= units::m / gOptGeomLUnits;

  gCfgUserAx1  *= units::rad / gOptGeomAUnits;
  gCfgUserAz   *= units::rad / gOptGeomAUnits;
  gCfgUserAx2  *= units::rad / gOptGeomAUnits;

  gCfgBoxLx    *= units::m / gOptGeomLUnits;
  gCfgBoxLy    *= units::m / gOptGeomLUnits;
  gCfgBoxLz    *= units::m / gOptGeomLUnits;

  gCfgParentTheta *= units::rad / gOptGeomAUnits;
  gCfgParentPhi   *= units::rad / gOptGeomAUnits;

  gCfgParentOx *= units::m / gOptGeomLUnits;
  gCfgParentOy *= units::m / gOptGeomLUnits;
  gCfgParentOz *= units::m / gOptGeomLUnits;

  ostringstream csts; 
  csts << "Read out the following config:"
       << "\n"
       << "\nNHL mass = " << gCfgMassNHL << " [GeV]"
       << "\n|U_e4|^2 = " << gCfgECoupling
       << "\n|U_m4|^2 = " << gCfgMCoupling
       << "\n|U_t4|^2 = " << gCfgTCoupling
       << "\n"
       << "\nProduction mode = " << utils::nhl::ProdAsString( gCfgProdMode )
       << "\nDecay mode      = " << utils::nhl::AsString( gCfgDecayMode )
       << "\nInteresting decay channels:";
  for( std::vector< NHLDecayMode_t >::iterator chit = gCfgIntChannels.begin();
       chit != gCfgIntChannels.end(); ++chit ){ csts << "\n\t" << utils::nhl::AsString(*chit); }
  csts << "\n"
       << "\nUser origin in beam coordinates = ( " << gCfgUserOx
       << ", " << gCfgUserOy << ", " << gCfgUserOz << " ) [" << lunits.c_str() << "]"
       << "\nEuler extrinsic x-z-x rotation = ( " << gCfgUserAx1
       << ", " << gCfgUserAz << ", " << gCfgUserAx2 << " ) [" << aunits.c_str() << "]"
       << "\nBox dimensions in user x-y-z: " << gCfgBoxLx << " x " << gCfgBoxLy
       << " x " << gCfgBoxLz << " [" << lunits.c_str() << "^3]"
       << "\n"
       << "\nParent energy = " << gCfgParentEnergy << " [GeV]"
       << "\nParent theta = " << gCfgParentTheta << " [" << aunits.c_str() << "]"
       << "\nParent phi   = " << gCfgParentPhi << " [" << aunits.c_str() << "]"
       << "\nParent origin = ( " << gCfgParentOx << ", " << gCfgParentOy << ", "
       << gCfgParentOz << " ) [" << lunits.c_str() << "]";

  LOG("gevald_nhl", pDEBUG) << csts.str();
  
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevald_nhl", pFATAL)
   << "\n **Syntax**"
   << "\n gevald_nhl [-h] "
   << "\n            [-r run#]"
   << "\n             -n n_of_events"
   << "\n            [-f path_to_flux_files]"
   << "\n            [-g geometry_file]"
   << "\n             -M mode:"
   << "\n                1: Flux prediction from dk2nu files. Needs -f option"
   << "\n                2: Flux prediction from histograms.  Needs -f option"
   << "\n                3: NHL production event rates from dk2nu. Needs -f option"
   << "\n                4: NHL production kinematics. Specify selected mode and"
   << "\n                   parent energy in config"
   << "\n                5: NHL decay event rates. Specify modes to uninhibit"
   << "\n                   in config"
   << "\n                6: NHL decay kinematics. Specify selected mode in config"
   << "\n                7: Custom geometry file entry/exit validation.  Needs -g option"
   << "\n                   Define trajectory in config"
   << "\n                8: Particle-gun for decay vertex placement + survival"
   << "\n                   probability (define energy, trajectory, and detector"
   << "\n                   box in config)"
   << "\n                9: Parent particle gun. Define energy and trajectory in config"
   << "\n               10: Full simulation (like gevgen_nhl but with lots of debug!)"
   << "\n"
   << "\n The configuration file lives at $GENIE/config/CommonNHL.xml - see"
   << " <param_set name=\"Validation\">"
   << "\n"
   << "\n Please also read the detailed documentation at http://www.genie-mc.org"
   << "\n or look at the source code: $GENIE/src/Apps/gNHLValidationApp.cxx"
   << "\n";
}
//_________________________________________________________________________________________
