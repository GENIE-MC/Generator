//________________________________________________________________________________________
/*!

\program gevald_hnl

\brief   Validation app for HNL generation.

         *** Synopsis :

         gevald_hnl [-h]
                    [-r run#]
                     -n n_of_events
		     -f path/to/flux/files
                    [-E hnl_energy]
		    [-g geometry (ROOT file)]
		    [-L geometry_length_units]
		    [-A geometry_angle_units]
		    [-T geometry_time_units]
		    [-o output_event_file_prefix]
		    [--seed random_number_seed]

         *** Options :

           [] Denotes an optional argument

           -h
              Prints out the gevald_hnl syntax and exits.
           -r
              Specifies the MC run number [default: 1000].
           -n
              Specifies how many events to generate.
           -f
              Input HNL flux.
	   -E
	      HNL energy for monoenergetic HNL
           -g
              Input detector geometry.
              If a geometry is specified, HNL decay vertices will be distributed
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

#include "Physics/BeamHNL/HNLDecayMode.h"
#include "Physics/BeamHNL/HNLDecayUtils.h"
#include "Physics/BeamHNL/HNLDecayVolume.h"
#include "Physics/BeamHNL/HNLFluxCreator.h"
//#include "Physics/BeamHNL/HNLFluxReader.h"
#include "Physics/BeamHNL/HNLDecayer.h"
#include "Physics/BeamHNL/HNLProductionMode.h"
#include "Physics/BeamHNL/SimpleHNL.h"

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
using namespace genie::HNL;

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

typedef enum t_HNLValidation {
  
  kValNone            = 0,
  kValFluxFromDk2nu   = 1,
  kValFluxFromHists   = 2,
  kValDecay           = 3,
  kValGeom            = 4,
  kValFullSim         = 5
  
} HNLValidation_t;

// function prototypes
void  GetCommandLineArgs (int argc, char ** argv);
void  ReadInConfig       (void);
void  PrintSyntax        (void);
const EventRecordVisitorI * HNLGenerator(void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
int      TestFluxFromDk2nu  (void);
int      TestFluxFromHists  (void);
GFluxI * TH1FluxDriver      (void);
#endif

int      TestDecay          (void);

#ifdef __CAN_USE_ROOT_GEOM__
int      TestGeom           (void);
void     InitBoundingBox    (void);
#endif // #ifdef __CAN_USE_ROOT_GEOM__

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomAUnits   = "rad";   // default geometry angle units
string          kDefOptGeomTUnits   = "ns";    // default geometry time units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units

NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format
string          kDefOptEvFilePrefix = "gntp";
string          kDefOptFluxFilePath = "./input-flux.root";

string          kDefOptSName   = "genie::EventGenerator";
string          kDefOptSConfig = "BeamHNL";

//
Long_t           gOptRunNu        = 1000;                // run number
int              gOptNev          = 10;                  // number of events to generate

HNLValidation_t  gOptValidationMode = kValNone;          // which test to run

// ---- this section set in config

double           gCfgMassHNL      = -1;                  // HNL mass
double           gCfgECoupling    = -1;                  // |U_e4|^2
double           gCfgMCoupling    = -1;                  // |U_m4|^2
double           gCfgTCoupling    = -1;                  // |U_t4|^2
bool             gCfgIsMajorana   = false;
int              gCfgHNLKind      = 2;


HNLProd_t        gCfgProdMode     = kHNLProdNull;        // HNL production mode
HNLDecayMode_t   gCfgDecayMode    = kHNLDcyNull;         // HNL decay mode
std::vector< HNLDecayMode_t > gCfgIntChannels;           // decays to un-inhibit

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

double           gCfgHNLCx        = -1;                  // directional cosine for x for HNL traj
double           gCfgHNLCy        = -1;                  // directional cosine for y for HNL traj
double           gCfgHNLCz        = -1;                  // directional cosine for z for HNL traj

// ---- end config section

double           gOptEnergyHNL    = -1;                  // HNL energy in GeV

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
string           gOptFluxFilePath = kDefOptFluxFilePath; // where flux files live
bool             gOptIsUsingDk2nu = false;               // using flat dk2nu files?
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
bool             gOptIsMonoEnFlux = true;                // do we have monoenergetic flux?

string           gOptEvFilePrefix = kDefOptEvFilePrefix; // event file prefix
bool             gOptUsingRootGeom = false;              // using root geom?
string           gOptRootGeom;                           // input ROOT file with realistic detector geometry

string           geom;
string           lunits, aunits, tunits;
double           gOptGeomLUnits = 0;                     // input geometry length units
double           gOptGeomAUnits = 0;                     // input geometry angle units
double           gOptGeomTUnits = 0;                     // input geometry time units
#ifdef __CAN_USE_ROOT_GEOM__
TGeoManager *    gOptRootGeoManager = 0;                 // the workhorse geometry manager
TGeoVolume  *    gOptRootGeoVolume  = 0;
#endif // #ifdef __CAN_USE_ROOT_GEOM__
string           gOptRootGeomTopVol = "";                // input geometry top event generation volume

// Geometry bounding box and origin - read from the input geometry file (if any)
double fdx = 0; // half-length - x
double fdy = 0; // half-length - y
double fdz = 0; // half-length - z
double fox = 0; // origin - x
double foy = 0; // origin - y
double foz = 0; // origin - z

long int         gOptRanSeed = -1;                       // random number seed

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

  LOG( "gevald_hnl", pFATAL )
    << "\nI am a work in progress. Hello world!"
    << "\n."
    << "\n. ."
    << "\n. . .";

  switch( gOptValidationMode ){
    
  case kValFluxFromDk2nu: return TestFluxFromDk2nu(); break;
  case kValFluxFromHists: return TestFluxFromHists(); break;
  case kValDecay:         return TestDecay();         break;
  case kValGeom:          return TestGeom();          break;
  default: LOG( "gevald_hnl", pFATAL ) << "I didn't recognise this mode. Goodbye world!"; break;
  }

  return 0;
}
//_________________________________________________________________________________________
//............................................................................
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
int TestFluxFromDk2nu()
{
  assert( !gOptIsMonoEnFlux && gOptIsUsingDk2nu );

  string foutName("test_flux_dk2nu.root");
  
  LOG( "gevald_hnl", pINFO )
    << "\n\nTesting flux prediction from dk2nu input files."
    << "\nWill produce 1 ROOT file (" << foutName << ") with:"
    << "\n--> Energy spectra for the HNL (total and broken down by parent)"
    << "\n--> Production vertex locations"
    << "\n--> Counters for each production mode"
    << "\n--> Spectrum of acceptance correction as function of parent boost factor"
    << "\n--> Boost factor spectrum of parents broken down by type";
  
  const Algorithm * algFluxCreator = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLFluxCreator", "Default");

  const HNLFluxCreator * fluxCreator = dynamic_cast< const HNLFluxCreator * >( algFluxCreator );

  fluxCreator->SetInputPath( gOptFluxFilePath );
  fluxCreator->SetGeomFile( gOptRootGeom );
  int maxFluxEntries = fluxCreator->GetNEntries();
  if( gOptNev > maxFluxEntries ){
    LOG( "gevald_hnl", pWARN )
      << "You have asked for " << gOptNev << " events, but only provided "
      << maxFluxEntries << " flux entries. Truncating events to " << maxFluxEntries << ".";
    gOptNev = maxFluxEntries;
  }

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 

  // RETHERE implement top volume option from cmd line
  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );
  TGeoShape * ts  = top_volume->GetShape();
  TGeoBBox *  box = (TGeoBBox *)ts;
  fluxCreator->ImportBoundingBox( box );

  TFile * fout = TFile::Open( foutName.c_str(), "RECREATE" );
  TH1D hEAll, hEPion, hEKaon, hEMuon, hENeuk;
  TH1D hPop, hImpwt;
  TH1D hAcceptanceCorr, hAcceptance, hAcceptNoBCorr;
  TH3D hProdVtxPos;
  TH1D hCounters;
  TH1D hBAll, hBPion, hBKaon, hBMuon, hBNeuk;
  TH1D hParamSpace; // to store mass + couplings

  hEAll  = TH1D( "hEAll",  "HNL energy - all parents", 1000, 0., 100. );
  hEPion = TH1D( "hEPion", "HNL energy - pion parent", 1000, 0., 100. );
  hEKaon = TH1D( "hEKaon", "HNL energy - kaon parent", 1000, 0., 100. );
  hEMuon = TH1D( "hEMuon", "HNL energy - muon parent", 1000, 0., 100. );
  hENeuk = TH1D( "hENeuk", "HNL energy - neuk parent", 1000, 0., 100. );

  hPop   = TH1D( "hPop",   "HNL populations in energy bins", 1000, 0., 100. );
  hImpwt = TH1D( "hImpwt", "HNL importance weights", 1000, 0., 100. );

  static const Double_t accbins[] = { 0.0, 2.5e-7, 5.0e-7, 7.5e-7, 1.0e-6, 2.5e-6, 5.0e-6, 7.5e-6, 1.0e-5, 2.5e-5, 5.0e-5, 7.5e-5, 1.0e-4, 2.5e-4, 5.0e-4, 7.5e-4, 1.0e-3, 2.5e-3, 5.0e-3, 7.5e-3, 1.0e-2, 2.5e-2, 5.0e-2, 7.5e-2, 1.0e-1, 2.5e-1, 5.0e-1, 7.5e-1, 1.0e+0, 2.5e+0, 5.0e+0, 7.5e+0, 1.0e+1 };
  const Int_t naccbins = sizeof(accbins)/sizeof(accbins[0]) - 1;
  hAcceptanceCorr = TH1D( "hAcceptanceCorr", "HNL acceptance correction", naccbins, accbins );
  hAcceptance     = TH1D( "hAcceptance",     "HNL acceptance"           , 200, 0.0, 100.0 );
  hAcceptNoBCorr  = TH1D( "hAcceptNoBCorr",  "HNL SAA * accCorr"        , 200, 0.0, 100.0 );
      
  hProdVtxPos = TH3D( "hProdVtxPos", "HNL production vertex (user coordinates, cm)",
		      200, -100., 100., 200, -100., 100., 1100, -110000., 0.);
  
  hCounters = TH1D( "hCounters", "HNL production channel counters", 11, 0, 11 );
  
  hBAll  = TH1D( "hBAll",  "Boost beta - all parents", 100, 0., 1. );
  hBPion = TH1D( "hBPion", "Boost beta - pion parent", 100, 0., 1. );
  hBKaon = TH1D( "hBKaon", "Boost beta - kaon parent", 100, 0., 1. );
  hBMuon = TH1D( "hBMuon", "Boost beta - muon parent", 100, 0., 1. );
  hBNeuk = TH1D( "hBNeuk", "Boost beta - neuk parent", 100, 0., 1. );

  hParamSpace = TH1D( "hParamSpace", "Parameter space", 5, 0., 5. );

  int parPDG;
  TLorentzVector p4HNL;
  TLorentzVector x4HNL;
  int nPion2Muon = 0, nPion2Electron = 0, nKaon2Muon = 0,
    nKaon2Electron = 0, nKaon3Muon = 0, nKaon3Electron = 0,
    nNeuk3Muon = 0, nNeuk3Electron = 0, nMuon3Numu = 0,
    nMuon3Nue = 0, nMuon3Nutau = 0;
  double nPOT;
  double betaMag;
  double accCorr;
  double nimpwt; // hadroproduction importance weight

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
      
      fluxCreator->SetCurrentEntry( ievent );

      EventRecord * event = new EventRecord;

      fluxCreator->ProcessEventRecord(event);
      flux::GNuMIFluxPassThroughInfo * gnmf = fluxCreator->RetrieveGNuMIFluxPassThroughInfo();
      
      // reject nonsense
      if( gnmf->necm < 0 ){
	
	LOG( "gevald_hnl", pDEBUG )
	  << "Skipping nonsense for event " << ievent << " (was this parent too light?)";

      } else {
	// now to make stuff from this... i.e. fill histos
	
	// first get the channel
	int iChannel = gnmf->ndecay - 30; int typeMod = 1;
	if( iChannel % 2 == 0 ){ iChannel--; typeMod = -1; }
	HNLGNuMIProd_t gChannel = ( HNLGNuMIProd_t ) iChannel;
	HNLProd_t pChannel;

	switch( gChannel ){
	case kHNLGProdNeuk3Electron:
	  parPDG = kPdgK0L; nNeuk3Electron++; pChannel = kHNLProdNeuk3Electron; break;
	case kHNLGProdNeuk3Muon:
	  parPDG = kPdgK0L; nNeuk3Muon++; pChannel = kHNLProdNeuk3Muon; break;
	case kHNLGProdKaon2Electron:
	  parPDG = kPdgKP; nKaon2Electron++; pChannel = kHNLProdKaon2Electron; break;
	case kHNLGProdKaon2Muon:
	  parPDG = kPdgKP; nKaon2Muon++; pChannel = kHNLProdKaon2Muon; break;
	case kHNLGProdKaon3Electron:
	  parPDG = kPdgKP; nKaon3Electron++; pChannel = kHNLProdKaon3Electron; break;
	case kHNLGProdKaon3Muon:
	  parPDG = kPdgKP; nKaon3Muon++; pChannel = kHNLProdKaon3Muon; break;
	case kHNLGProdMuon3Nue:
	  parPDG = kPdgMuon; nMuon3Nue++; pChannel = kHNLProdMuon3Nue; break;
	case kHNLGProdMuon3Numu:
	  parPDG = kPdgMuon; nMuon3Numu++; pChannel = kHNLProdMuon3Numu; break;
	case kHNLGProdMuon3Nutau:
	  parPDG = kPdgMuon; nMuon3Nutau++; pChannel = kHNLProdMuon3Nutau; break;
	case kHNLGProdPion2Electron:
	  parPDG = kPdgPiP; nPion2Electron++; pChannel = kHNLProdPion2Electron; break;
	case kHNLGProdPion2Muon:
	  parPDG = kPdgPiP; nPion2Muon++; pChannel = kHNLProdPion2Muon; break;
	default:
	  LOG( "gevald_hnl", pERROR )
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
	
	p4HNL = gnmf->fgP4User;
	x4HNL = gnmf->fgX4User; // in cm

	double acceptance = gnmf->fgXYWgt; // full acceptance
	double bCorr = gnmf->nenergyn; // boost correction
	accCorr    = gnmf->nwtnear;   // just the acceptance correction
	nimpwt = gnmf->nimpwt;

	double gam = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );
	double gamStar = gnmf->necm / gCfgMassHNL;
	double betaStar = std::sqrt( 1.0 - 1.0 / ( gamStar * gamStar ) );
	double bigBeta = betaMag * gam / betaStar;

	double Ox = -1.0 * x4HNL.X();
	double Oy = -1.0 * x4HNL.Y();
	double Oz = -1.0 * x4HNL.Z();

	double rootArg = gam*gam*Oz*Oz - 
	  (gam*gam - bigBeta*bigBeta)*(Oz*Oz - bigBeta*bigBeta*(Ox*Ox + Oy*Oy));

	double timelikeBit = ( rootArg >= 0.0 ) ? std::sqrt( rootArg ) / ( betaMag * gam * gam * gam ) : 0.0;
	
	double fullTerm = (1.0 - 1.0 / gam) * Oz / ( betaMag * gam ) - timelikeBit;
	if( std::abs( fullTerm ) > 100.0 ) fullTerm *= 100.0 / std::abs( fullTerm );
	
	nPOT = gnmf->norig;
	
	// fill the histos!
	hEAll.Fill( p4HNL.E(), acceptance * nimpwt );
	hProdVtxPos.Fill( x4HNL.X(), x4HNL.Y(), x4HNL.Z(), nimpwt );
	hBAll.Fill( betaMag, nimpwt );
	  
	hPop.Fill( p4HNL.E(), 1.0 );
	hImpwt.Fill( p4HNL.E(), nimpwt );
	hAcceptanceCorr.Fill( accCorr, nimpwt );
	hAcceptance.Fill( p4HNL.E(), nimpwt * acceptance );
	hAcceptNoBCorr.Fill( p4HNL.E(), nimpwt * acceptance / ( bCorr * bCorr ) );
	
	switch( parPDG ){
	case kPdgPiP:
	  hEPion.Fill( p4HNL.E(), acceptance * nimpwt ); 
	  hBPion.Fill( betaMag, nimpwt );
	  break;
	case kPdgKP:
	  hEKaon.Fill( p4HNL.E(), acceptance * nimpwt ); 
	  hBKaon.Fill( betaMag, nimpwt );
	  break;
	case kPdgK0L:
	  hENeuk.Fill( p4HNL.E(), acceptance * nimpwt ); 
	  hBNeuk.Fill( betaMag, nimpwt );
	  break;
	case kPdgMuon:
	  hEMuon.Fill( p4HNL.E(), acceptance * nimpwt ); 
	  hBMuon.Fill( betaMag, nimpwt );
	  break;
	}
	
	LOG( "gevald_hnl", pDEBUG )
	  << " *** Output for event no " << ievent << "... ***"
	  << "\nparPDG = " << parPDG
	  << "\np4par = " << utils::print::P4AsString( &p4par ) << " [GeV] "
	  << "\nbetaVec = " << utils::print::Vec3AsString( &boost_beta )
	  << " : mag = " << betaMag
	  << "\np4HNL = " << utils::print::P4AsString( &p4HNL ) << " [GeV] "
	  << "\nx4HNL = " << utils::print::X4AsString( &x4HNL ) << " [cm]"
	  << "\nacceptance = " << acceptance << " : accCorr = " << accCorr
	  << "\nnimpwt = " << nimpwt
	  << "\nProduction channel = " << utils::hnl::ProdAsString( pChannel );

      } // if not nonsense

      // clean up
      delete event;

      ievent++;
    }

  // fill hCounters; each bin corresponds to the HNLProd_t with index iBin-1
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

  hParamSpace.SetBinContent( 1, 1000.0 * gCfgMassHNL ); // MeV
  hParamSpace.SetBinContent( 2, gCfgECoupling );
  hParamSpace.SetBinContent( 3, gCfgMCoupling );
  hParamSpace.SetBinContent( 4, gCfgTCoupling );
  hParamSpace.SetBinContent( 5, nPOT );

  fout->Write();
  fout->Close();

  return 0;
}
//_________________________________________________________________________________________
int TestFluxFromHists()
{
  assert( !gOptIsMonoEnFlux && !gOptIsUsingDk2nu );

  string foutName("test_flux_hists.root");

  LOG( "gevald_hnl", pINFO )
    << "\n\nTesting flux prediction from precomputed flux spectra."
    << "\nWill produce 1 ROOT file ( " << foutName << ") with:"
    << "\n--> Energy spectrum for the HNL"
    << "\n--> Momentum spectra on x, y, z (user)"
    << "\n--> Angular deviation from parent spectrum"
    << "\n--> Rates of particle vs antiparticle";

  TFile * fout = TFile::Open( foutName.c_str(), "RECREATE" );

  const EventRecordVisitorI * mcgen = HNLGenerator();
  const Algorithm * algFluxCreator = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLFluxCreator", "Default");
  const HNLFluxCreator * fluxCreator = dynamic_cast< const HNLFluxCreator * >( algFluxCreator );

  assert(gCfgMassHNL >= 0.0);
  
  string finPath = gOptFluxFilePath; finPath.append("./histFluxes.root");
  string prodVtxPath = gOptFluxFilePath; prodVtxPath.append("/HNL_vertex_positions.root");
  __attribute__((unused)) int iset = setenv( "PRODVTXDIR", prodVtxPath.c_str(), 1 );
  LOG("gevald_hnl", pDEBUG)
    << "Looking for fluxes in " << finPath.c_str();
  assert( !gSystem->AccessPathName( finPath.c_str()) );
  
  fluxCreator->SetFinPaths( finPath, prodVtxPath );
  fluxCreator->BuildInputFlux();

  TFile * spectrumFile = TFile::Open("./input-flux.root", "READ");
  TDirectory * baseDir = spectrumFile->GetDirectory("");
  std::string fluxName = std::string( "spectrum" );
  assert( baseDir->GetListOfKeys()->Contains( fluxName.c_str() ) );
  TH1D * spectrum = ( TH1D * ) baseDir->Get( fluxName.c_str() );
  assert( spectrum );

  TH1D hHNLPx, hHNLPy, hHNLPz;
  TH1D hHNLAngDev, hHNLPhi;
  TH1D hHNLParticleRates; int nPart = 0, nAntipart = 0;
  TH1D hParamSpace;

  hHNLPx = TH1D( "hHNLPx", "HNL p_x (user coordinates, GeV)", 100, -0.5, 0.5 );
  hHNLPy = TH1D( "hHNLPy", "HNL p_y (user coordinates, GeV)", 100, -0.5, 0.5 );
  hHNLPz = TH1D( "hHNLPz", "HNL p_z (user coordinates, GeV)", 1050, -0.5, 100 );

  double angdev = utils::hnl::GetCfgDouble( "HNL", "InitialState", "HNL-angular_deviation" );
  hHNLAngDev = TH1D( "hHNLAngDev", "HNL angular deviation [deg]", 100, -5.0 * angdev, 5.0 * angdev );
  hHNLPhi    = TH1D( "hHNLPhi", "HNL #phi [deg]", 100, 0.0, 360.0 );

  hHNLParticleRates = TH1D( "hHNLParticleRates", "N (particles && antiparticles)", 2, 0., 2. );
  hParamSpace = TH1D( "hParamSpace", "Parameter space", 5, 0., 5. );

  // Event loop
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

      if( ievent == gOptNev ){ std::cerr << " \n"; break; }

      gOptEnergyHNL = 100.0;
      
      int hpdg = genie::kPdgHNL;
      
      EventRecord * event = new EventRecord;
      Interaction * interaction = Interaction::HNL( hpdg, gOptEnergyHNL, kHNLDcyTEST );
      event->AttachSummary( interaction );

      fluxCreator->SetUsingDk2nu( false );
      fluxCreator->ProcessEventRecord(event);
      
      // now grab momentum from event
      const TLorentzVector * p4HNL = event->Particle(0)->GetP4();
      gOptEnergyHNL = p4HNL->E();

      interaction = event->Summary(); // it's been updated now
      if( event->Particle(0)->Pdg() > 0 && !gCfgIsMajorana ) nPart++;
      else if( event->Particle(0)->Pdg() < 0 && !gCfgIsMajorana ) nAntipart++;
      else{ nPart++; nAntipart++; }

      LOG( "gevald_hnl", pDEBUG )
	<< "*** Event " << ievent << ":"
	<< "\n!*!*!* p4HNL = " << utils::print::P4AsString( p4HNL );
      
      double px = p4HNL->Px(), py = p4HNL->Py(), pz = p4HNL->Pz();
      double theta = TMath::ACos( pz / p4HNL->P() );
      double phi = TMath::ACos( px / ( p4HNL->P() * TMath::Sin( theta ) ) );
      if( py < 0.0 ) phi = 2.0 * constants::kPi - phi;
      if( constants::kPi <= phi && phi < 2.0 * constants::kPi ) theta *= -1;
      
      // fill histos here

      hHNLPx.Fill( px, 1.0 );
      hHNLPy.Fill( py, 1.0 );
      hHNLPz.Fill( pz, 1.0 );
      hHNLAngDev.Fill( theta * 180.0 / constants::kPi, 1.0 );
      hHNLPhi.Fill( phi * 180.0 / constants::kPi, 1.0 );
      
      delete event;
      
      ievent++;
    } // event loop

  hParamSpace.SetBinContent( 1, 1000.0 * gCfgMassHNL ); // MeV
  hParamSpace.SetBinContent( 2, gCfgECoupling );
  hParamSpace.SetBinContent( 3, gCfgMCoupling );
  hParamSpace.SetBinContent( 4, gCfgTCoupling );

  hHNLParticleRates.SetBinContent( 1, nPart );
  hHNLParticleRates.SetBinContent( 2, nAntipart );

  LOG( "gevald_hnl", pDEBUG )
    << "\nnPart, nAntipart = " << nPart << ", " << nAntipart;

  fout->cd();
  hHNLPx.Write();
  hHNLPy.Write();
  hHNLPz.Write();
  hHNLAngDev.Write();
  hHNLPhi.Write();
  hHNLParticleRates.Write();
  hParamSpace.Write();
  fout->Write();
  fout->Close();
  
  return 0;
}
//_________________________________________________________________________________________
GFluxI * TH1FluxDriver(void)
{
  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algFluxCreator = algf->GetAlgorithm("genie::EventGenerator", "BeamHNL");

  const HNLFluxCreator * fluxCreator = 
    dynamic_cast< const HNLFluxCreator * >( algFluxCreator );

  //
  //
  flux::GCylindTH1Flux * flux = new flux::GCylindTH1Flux;
  TH1D * spectrum = 0;

  double emin = 0.0; 
  double emax = utils::hnl::GetCfgDouble( "HNL", "InitialState", "HNL-max-energy" ); 

  TVector3 bdir (0.0,0.0,1.0);
  TVector3 bspot(0.0,0.0,1.0);

  flux->SetNuDirection      (bdir);
  flux->SetBeamSpot         (bspot);
  flux->SetTransverseRadius (-1);
  flux->AddEnergySpectrum   (genie::kPdgHNL, spectrum);

  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  LOG("gevald_hnl", pDEBUG)
    << "Returning flux driver and exiting method.";
  return flux_driver;
}
//............................................................................
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_
//_________________________________________________________________________________________
int TestDecay(void)
{
  string foutName("test_decay.root");

  TFile * fout = TFile::Open( foutName.c_str(), "RECREATE" );

  const EventRecordVisitorI * mcgen = HNLGenerator();
  const Algorithm * algHNLGen = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLDecayer", "Default");
  const HNLDecayer * hnlgen = dynamic_cast< const HNLDecayer * >( algHNLGen );

  LOG( "gevald_hnl", pDEBUG ) << "gOptRootGeom = " << gOptRootGeom;
    
  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 
  
  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );
  TGeoShape * ts  = top_volume->GetShape();
  TGeoBBox *  box = (TGeoBBox *)ts;
  
  LOG( "gevald_hnl", pDEBUG ) << "Imported box.";
  
  const Algorithm * algDkVol = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLDecayVolume", "Default");
  
  const HNLDecayVolume * dkVol = dynamic_cast< const HNLDecayVolume * >( algDkVol );
  dkVol->ImportBoundingBox( box );
  
  SimpleHNL sh = SimpleHNL( "HNLInstance", 0, kPdgHNL, kPdgKP,
			    gCfgMassHNL, gCfgECoupling, gCfgMCoupling, gCfgTCoupling, false );
  std::map< HNLDecayMode_t, double > valMap = sh.GetValidChannels();
  const double CoMLifetime = sh.GetCoMLifetime();

  assert( valMap.size() > 0 ); // must be able to decay to something!
  assert( (*valMap.begin()).first == kHNLDcyNuNuNu );

  LOG( "gevald_hnl", pINFO )
    << "\n\nTesting decay modes for the HNL."
    << "\nWill process gOptNev = " << gOptNev << " x ( N_channels = " << valMap.size()
    << " ) = " << valMap.size() * gOptNev << " events, 1 for each valid channel."
    << "\nWill produce 1 ROOT file ( " << foutName << " ) with:"
    << "\n--> Energy spectrum for the decay products for each channel"
    << "\n--> Rates of each decay channel";

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  // first set the 4-momentum of the HNL
  gOptEnergyHNL = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-Energy" );
  double p3HNL = std::sqrt( gOptEnergyHNL * gOptEnergyHNL - gCfgMassHNL * gCfgMassHNL );
  assert( p3HNL >= 0.0 );
  TLorentzVector * p4HNL = new TLorentzVector( p3HNL * gCfgHNLCx, 
					       p3HNL * gCfgHNLCy, 
					       p3HNL * gCfgHNLCz, gOptEnergyHNL );

  LOG( "gevald_hnl", pDEBUG )
    << "\nUsing HNL with 4-momentum " << utils::print::P4AsString( p4HNL );
  sh.SetEnergy( gOptEnergyHNL );
  sh.SetMomentumDirection( gCfgHNLCx, gCfgHNLCy, gCfgHNLCz );

  TLorentzVector * x4HNL = new TLorentzVector( 1.0, 2.0, 3.0, 0.0 ); // dummy

  // now build array with indices of valid decay modes for speedy access
  HNLDecayMode_t validModes[10] = { kHNLDcyNull, kHNLDcyNull, kHNLDcyNull, kHNLDcyNull, kHNLDcyNull, kHNLDcyNull, kHNLDcyNull, kHNLDcyNull, kHNLDcyNull, kHNLDcyNull };
  double validRates[10] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  std::map< HNLDecayMode_t, double >::iterator vmit = valMap.begin(); int modeIdx = 0;
  std::ostringstream msts;
  for( ; vmit != valMap.end(); ++vmit ){
    validModes[ modeIdx ] = (*vmit).first;
    validRates[ modeIdx ] = (*vmit).second;
    msts << "\n" << utils::hnl::AsString( (*vmit).first );
    modeIdx++;
  }

  LOG( "gevald_hnl", pDEBUG ) << "Here are the modes in order : " << msts.str();

  // declare histos
  // hSpectrum[i][j]: i iterates over HNLDecayMode_t, j over FS particle in same order as event record
  TH1D hSpectrum[10][3], hCMSpectrum[10][3], hRates;
  TH1D hParamSpace = TH1D( "hParamSpace", "Parameter space", 5, 0., 5. );

  hParamSpace.SetBinContent( 1, 1000.0 * gCfgMassHNL ); // MeV
  hParamSpace.SetBinContent( 2, gCfgECoupling );
  hParamSpace.SetBinContent( 3, gCfgMCoupling );
  hParamSpace.SetBinContent( 4, gCfgTCoupling );

  hRates = TH1D( "hRates", "Rates of HNL decay channels", 10, 0, 10 );

  std::string shortModes[10] = { "vvv", "vee", "vmue", "pi0v", "pie", "vmumu", "pimu", "pi0pi0v", 
				 "pipi0e", "pipi0mu" };
  std::string part0names[10] = { "v1", "v", "v", "pi0", "pi", "v", "pi", "pi01", "pi", "pi" };
  std::string part1names[10] = { "v2", "e1", "mu", "v", "e", "mu1", "mu", "pi02", "pi0", "pi0" };
  std::string part2names[10] = { "v3", "e2", "e", "None", "None", "mu2", "None", "v", "e", "mu" };
  std::string partNames[3][10] = { part0names, part1names, part2names };
  for( Int_t iChan = 0; iChan < valMap.size(); iChan++ ){

    std::string shortMode = shortModes[iChan];

    for( Int_t iPart = 0; iPart < 3; iPart++ ){
      std::string ParticleName = partNames[iPart][iChan];
      if( strcmp( ParticleName.c_str(), "None" ) != 0 ){
	hSpectrum[iChan][iPart]   = TH1D( Form( "hSpectrum_%s_%s", shortMode.c_str(), ParticleName.c_str() ),
					  Form( "Fractional energy of particle: %s  in decay: %s", 
						ParticleName.c_str(), 
						(utils::hnl::AsString( validModes[iChan] )).c_str() ), 
					  100, 0., 1.0 );
	hCMSpectrum[iChan][iPart] = TH1D( Form( "hCMSpectrum_%s_%s", shortMode.c_str(), ParticleName.c_str() ),
					  Form( "Rest frame energy of particle: %s  in decay: %s", 
						ParticleName.c_str(), 
						(utils::hnl::AsString( validModes[iChan] )).c_str() ), 
					  100, 0., gCfgMassHNL );
      } // only declare histos of particles that exist in decay
        // and are of allowed decays
    }
  }

  // fill hRates now
  double rnununu = ( (*valMap.find( kHNLDcyNuNuNu )) ).second;
  double rnuee = ( valMap.find( kHNLDcyNuEE ) != valMap.end() ) ? (*(valMap.find( kHNLDcyNuEE ))).second : 0.0;
  double rnumue = ( valMap.find( kHNLDcyNuMuE ) != valMap.end() ) ? (*(valMap.find( kHNLDcyNuMuE ))).second : 0.0;
  double rpi0nu = ( valMap.find( kHNLDcyPi0Nu ) != valMap.end() ) ? (*(valMap.find( kHNLDcyPi0Nu ))).second : 0.0;
  double rpie = ( valMap.find( kHNLDcyPiE ) != valMap.end() ) ? (*(valMap.find( kHNLDcyPiE ))).second : 0.0;
  double rnumumu = ( valMap.find( kHNLDcyNuMuMu ) != valMap.end() ) ? (*(valMap.find( kHNLDcyNuMuMu ))).second : 0.0;
  double rpimu = ( valMap.find( kHNLDcyPiMu ) != valMap.end() ) ? (*(valMap.find( kHNLDcyPiMu ))).second : 0.0;
  double rpi0pi0nu = ( valMap.find( kHNLDcyPi0Pi0Nu ) != valMap.end() ) ? (*(valMap.find( kHNLDcyPi0Pi0Nu ))).second : 0.0;
  double rpipi0e = ( valMap.find( kHNLDcyPiPi0E ) != valMap.end() ) ? (*(valMap.find( kHNLDcyPiPi0E ))).second : 0.0;
  double rpipi0mu = ( valMap.find( kHNLDcyPiPi0Mu ) != valMap.end() ) ? (*(valMap.find( kHNLDcyPiPi0Mu ))).second : 0.0;

  hRates.SetBinContent(  1, rpimu );
  hRates.SetBinContent(  2, rpie );
  hRates.SetBinContent(  3, rpi0nu );
  hRates.SetBinContent(  4, rnununu );
  hRates.SetBinContent(  5, rnumumu );
  hRates.SetBinContent(  6, rnuee );
  hRates.SetBinContent(  7, rnumue );
  hRates.SetBinContent(  8, rpipi0e );
  hRates.SetBinContent(  9, rpipi0mu );
  hRates.SetBinContent( 10, rpi0pi0nu );

  int ievent = 0;
  while( true ){
    
    if( gOptNev >= 10000 ){
      if( ievent % (gOptNev / 1000) == 0 ){
	int irat = ievent / (gOptNev / 1000);
	std::cerr << Form("%2.2f", 0.1 * irat) << " % ( " << ievent << " / "
		  << gOptNev << " ) \r" << std::flush;
      }
    }
    
    if( ievent == gOptNev ){ std::cerr << " \n"; break; }
    
    ostringstream asts;
    for( Int_t iMode = 0; iMode < valMap.size(); iMode++ ){
      if( ievent == 0 ){
	asts
	  << "\nDecay mode " << iMode << " is " << utils::hnl::AsString( validModes[ iMode ] );
      }

      // build an event
      EventRecord * event = new EventRecord;
      Interaction * interaction = Interaction::HNL( genie::kPdgHNL, gOptEnergyHNL, validModes[ iMode ] );
      // set p4 and a dummy vertex so HNLDecayer doesn't attempt to regenerate init state
      interaction->InitStatePtr()->SetProbeP4( *p4HNL );
      event->SetVertex( *x4HNL );

      event->AttachSummary( interaction );
      LOG( "gevald_hnl", pDEBUG )
	<< "Simulating decay with mode " 
	<< utils::hnl::AsString( (HNLDecayMode_t) interaction->ExclTag().DecayMode() );

      // simulate the decay
      hnlgen->ProcessEventRecord( event );

      dkVol->SetStartingParameters( event, CoMLifetime, false, true, gOptRootGeom.c_str() );

      // now get a weight.
      // = exp( - T_{box} / \tau_{HNL} ) = exp( - L_{box} / ( \beta_{HNL} \gamma_{HNL} c ) * h / \Gamma_{tot} )
      // placing the HNL at a point configured by the user
      //dkVol->MakeSDV();
      double ox = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginX" );
      double oy = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginY" );
      double oz = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginZ" );
      ox *= units::m / units::mm; oy *= units::m / units::mm; oz *= units::m / units::mm;
      TVector3 startPoint( ox, oy, oz ); TVector3 entryPoint, exitPoint;
      TLorentzVector tmpVtx( ox, oy, oz, 0.0 );
      event->SetVertex( tmpVtx );
      /*
      std::string dummyGeom = "";
      dkVol->SetStartingParameters( event, 1.0e+20, false, false, dummyGeom );
      */

      dkVol->ProcessEventRecord(event);

      LOG( "gevald_hnl", pDEBUG ) << *event;

      // now fill the histos!
      double wgt = event->Weight();
      hSpectrum[iMode][0].Fill( (event->Particle(1))->E() / gOptEnergyHNL, wgt );
      hSpectrum[iMode][1].Fill( (event->Particle(2))->E() / gOptEnergyHNL, wgt );
      if( event->Particle(3) ) hSpectrum[iMode][2].Fill( (event->Particle(3))->E() / gOptEnergyHNL, wgt );

      // let's also fill the CM spectra
      // get particle 4-momenta and boost back to rest frame!
      TLorentzVector * p4p1 = (event->Particle(1))->GetP4();
      TLorentzVector * p4p2 = (event->Particle(2))->GetP4();
      TLorentzVector * p4p3 = 0;
      if( event->Particle(3) ) p4p3 = (event->Particle(3))->GetP4();

      TVector3 boostVec = p4HNL->BoostVector();

      p4p1->Boost( -boostVec );
      p4p2->Boost( -boostVec );
      if( p4p3 ) p4p3->Boost( -boostVec );

      hCMSpectrum[iMode][0].Fill( p4p1->E(), wgt );
      hCMSpectrum[iMode][1].Fill( p4p2->E(), wgt );
      if( p4p3 ) hCMSpectrum[iMode][2].Fill( p4p3->E(), wgt );

      // clean-up
      delete event;

    } // loop over valid decay channels
    
    if( ievent == 0 ) LOG( "gevald_hnl", pDEBUG ) << asts.str();
    LOG( "gevald_hnl", pDEBUG ) << "Finished event: " << ievent;

    ievent++;
  
  } // event loop

  fout->cd();
  hParamSpace.Write();
  hRates.Write();
  for( Int_t i = 0; i < valMap.size(); i++ ){
    for( Int_t j = 0; j < 3; j++ ){
      std::string ParticleName = partNames[j][i];
      if( strcmp( ParticleName.c_str(), "None" ) != 0 ){
	hSpectrum[i][j].Write();
	hCMSpectrum[i][j].Write();
      }
    }
  }
  fout->Write();
  fout->Close();
  
  delete p4HNL;
  delete x4HNL;
  return 0;
}
//............................................................................
#ifdef __CAN_USE_ROOT_GEOM__
//_________________________________________________________________________________________
int TestGeom(void)
{
  string foutName("test_geom.root");

  TFile * fout = TFile::Open( foutName.c_str(), "RECREATE" );

  LOG( "gevald_hnl", pINFO )
    << "\n\nTesting ROOT geometry for ROOT file " << gOptRootGeom
    << "\nWill produce 1 ROOT file ( " << foutName << " ) with:"
    << "\n--> Specified geometry"
    << "\n--> TTree containing the following branch structure:"
    << "\n    |"
    << "\n    |---- start x,y,z [mm]"
    << "\n    |"
    << "\n    |---- HNL 4-momentum [GeV]"
    << "\n    |"
    << "\n    |---- did intersect detector?"
    << "\n    |"
    << "\n    |---- entry x,y,z [mm]"
    << "\n    |"
    << "\n    |---- exit  x,y,z [mm]"
    << "\n    |"
    << "\n    |---- weight"
    << "\n    |"
    << "\n    |---- HNL lifetime (rest frame) [ns]"
    << "\n    |"
    << "\n    |---- HNL lifetime (lab  frame) [ns]"
    << "\n    |"
    << "\n    |---- decay x,y,z [mm]"
    << "\n--> Weight histogram"
    << "\n--> Travel-length histogram"
    << "\n"
    << "\n(where same coordinate system as user's is used and weight == P(HNL decays in detector) )";

  double dev_start[3]    = { -9999.9, -9999.9, -9999.9 };
  double dev_sphere[2]   = { -9999.9, -9999.9 };
  double use_start[3]    = { -9999.9, -9999.9, -9999.9 };
  double use_entry[3]    = { -9999.9, -9999.9, -9999.9 };
  double use_exit[3]     = { -9999.9, -9999.9, -9999.9 };
  double use_decay[3]    = { -9999.9, -9999.9, -9999.9 };
  double use_momentum[4] = { -9999.9, -9999.9, -9999.9, -9999.9 };
  double use_wgt = -9999.9;
  double use_lifetime = -9999.9, use_CMlifetime = -9999.9;
  bool   didIntersectDet = false;
  
  TTree * outTree = new TTree( "outTree", "Trajectory information tree" );
  outTree->Branch( "startPoint",   use_start,        "startPoint[3]/D"   );
  outTree->Branch( "fourMomentum", use_momentum,     "fourMomentum[4]/D" );
  outTree->Branch( "startDeviate", dev_start,        "startDeviate[3]/D" );
  outTree->Branch( "spherDeviate", dev_sphere,       "spherDeviate[2]/D" );
  outTree->Branch( "didIntersect", &didIntersectDet, "didIntersect/O"    );
  outTree->Branch( "entryPoint",   use_entry,        "entryPoint[3]/D"   );
  outTree->Branch( "exitPoint",    use_exit,         "exitPoint[3]/D"    );
  outTree->Branch( "decayPoint",   use_decay,        "decayPoint[3]/D"   );
  outTree->Branch( "weight",       &use_wgt,         "weight/D"          );
  outTree->Branch( "lifetime_LAB", &use_lifetime,    "lifetime_LAB/D"    );
  outTree->Branch( "lifetime_CM",  &use_CMlifetime,  "lifetime_CM/D"     );

  TH1D hWeight( "hWeight", "Log_{10} ( P( decay in detector ) )", 
		80, -7.0, 1.0  );
  TH1D hLength( "hLength", "Length travelled in detector / max possible length in detector",
		100, 0., 1.0 );

  gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str());
  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );

  // Read geometry bounding box - for vertex position generation
  InitBoundingBox();

  const Algorithm * algDkVol = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLDecayVolume", "Default");
  const HNLDecayVolume * dkVol = dynamic_cast< const HNLDecayVolume * >( algDkVol );

  // get SimpleHNL for lifetime
  SimpleHNL sh = SimpleHNL( "HNLInstance", 0, kPdgHNL, kPdgKP,
			    gCfgMassHNL, gCfgECoupling, gCfgMCoupling, gCfgTCoupling, false );

  // first set the 4-momentum of the HNL
  gOptEnergyHNL = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-Energy" );
  double p3HNL = std::sqrt( gOptEnergyHNL * gOptEnergyHNL - gCfgMassHNL * gCfgMassHNL );
  assert( p3HNL >= 0.0 );
  TLorentzVector * p4HNL = new TLorentzVector( p3HNL * gCfgHNLCx, 
					       p3HNL * gCfgHNLCy, 
					       p3HNL * gCfgHNLCz, gOptEnergyHNL );
  
  LOG( "gevald_hnl", pDEBUG )
    << "\nUsing HNL with 4-momentum " << utils::print::P4AsString( p4HNL );
  sh.SetEnergy( gOptEnergyHNL );
  sh.SetMomentumDirection( gCfgHNLCx, gCfgHNLCy, gCfgHNLCz );

  double betaMag = p4HNL->P() / p4HNL->E();
  double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );

  use_CMlifetime = sh.GetCoMLifetime() / ( units::ns * units::GeV );
  use_lifetime   = sh.GetLifetime() / ( units::ns * units::GeV ); // ns

  const double PGox = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginX" );
  const double PGoy = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginY" );
  const double PGoz = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginZ" ); // m

  const double PGdx = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginDX" );
  const double PGdy = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginDY" );
  const double PGdz = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-OriginDZ" ); // m
  assert( PGdx > 0.0 && PGdy > 0.0 && PGdz > 0.0 );

  double c2 = std::sqrt( std::pow( gCfgHNLCx, 2.0 ) + std::pow( gCfgHNLCy, 2.0 ) + std::pow( gCfgHNLCz, 2.0 ) );
  const double PGcx = gCfgHNLCx / c2;
  const double PGcy = gCfgHNLCy / c2;
  const double PGcz = gCfgHNLCz / c2; // unit-normalised

  const double PGtheta = std::acos( PGcz );
  const double PGphi = ( std::sin( PGtheta ) < controls::kASmallNum ) ? 0.0 :
    ( PGcy >= 0.0 ) ? std::acos( PGcx / PGcz ) : 2.0 * constants::kPi - std::acos( PGcx / PGcz );

  const double PGdtheta = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-DTheta" ) * constants::kPi / 180.0;
  const double PGdphi   = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-DPhi" ) * constants::kPi / 180.0; // rad

  /*
   * The event loop works out a bit differently here.
   * It reads in the origin point and momentum from file and treats them as CV
   * Then it partitions each axis in R^3 (x-y-z, for origin position) and R^2 (theta-phi, for
   * origin momentum) into either 5 points (x,y,z) or 9 points (theta, phi)
   * There is exactly one trajectory that matches both origin position and momentum.
   * For each of these trajectories, the entry and exit points are calculated.
   * The length to decay is also calculated and a histo is filled with L(to decay) / L(max).
   */

  // first make the points
  const int NCARTESIAN = 5;
  const int NSPHERICAL = 9;
  const int NMAX = NCARTESIAN * NCARTESIAN * NCARTESIAN * NSPHERICAL * NSPHERICAL;

  double arr_ox[ NCARTESIAN ] = { PGox - PGdx, PGox - PGdx/2.0, PGox, PGox + PGdx/2.0, PGox + PGdx };
  double arr_oy[ NCARTESIAN ] = { PGoy - PGdy, PGoy - PGdy/2.0, PGoy, PGoy + PGdy/2.0, PGoy + PGdy };
  double arr_oz[ NCARTESIAN ] = { PGoz - PGdz, PGoz - PGdz/2.0, PGoz, PGoz + PGdz/2.0, PGoz + PGdz };

  double arr_theta[ NSPHERICAL ] = { PGtheta - PGdtheta, PGtheta - 3.0/4.0 * PGdtheta, PGtheta - 1.0/2.0 * PGdtheta, PGtheta - 1.0/4.0 * PGdtheta, PGtheta, PGtheta + 1.0/4.0 * PGdtheta, PGtheta + 1.0/2.0 * PGdtheta, PGtheta + 3.0/4.0 * PGdtheta, PGtheta + PGdtheta };
  double arr_phi[ NSPHERICAL ] = { PGphi - PGdphi, PGphi - 3.0/4.0 * PGdphi, PGphi - 1.0/2.0 * PGdphi, PGphi - 1.0/4.0 * PGdphi, PGphi, PGphi + 1.0/4.0 * PGdphi, PGphi + 1.0/2.0 * PGdphi, PGphi + 3.0/4.0 * PGdphi, PGphi + PGdphi };

  // so now we have NCARTESIAN ^3 x NSPHERICAL ^2 points to iterate over. That's 10125 events for 5 and 9

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 
  
  TGeoShape * ts  = top_volume->GetShape();
  
  TGeoBBox *  box = (TGeoBBox *)ts;

  // pass this box to HNLDecayVolume
  dkVol->ImportBoundingBox( box );
  
  int ievent = 0;
  ostringstream asts;
  while( true ){

    if( ievent == NMAX ) break;

    LOG( "gevald_hnl", pDEBUG )
      << "*** Building event = " << ievent;

    EventRecord * event = new EventRecord;
    Interaction * interaction = Interaction::HNL( genie::kPdgHNL, gOptEnergyHNL, HNL::kHNLDcyTEST );
    event->AttachSummary( interaction );

    /*
     * Iterate over events as follows:
     * Least significant --> most significant (with period in brackets)
     * ox (NCART) --> oy (NCART) --> oz (NCART) --> phi (NSPHE) --> theta (NSPHE)
     */
    
    int ix = ievent % NCARTESIAN;
    int iy = ( ievent / NCARTESIAN ) % NCARTESIAN;
    int iz = ( ievent / NCARTESIAN / NCARTESIAN ) % NCARTESIAN;
    int ip = ( ievent / NCARTESIAN / NCARTESIAN / NCARTESIAN ) % NSPHERICAL;
    int it = ( ievent / NCARTESIAN / NCARTESIAN / NCARTESIAN / NSPHERICAL ) % NSPHERICAL;

    double use_ox = arr_ox[ ix ] * units::m / units::mm;
    double use_oy = arr_oy[ iy ] * units::m / units::mm;
    double use_oz = arr_oz[ iz ] * units::m / units::mm;

    dev_start[0] = use_ox - PGox * units::m / units::mm;
    dev_start[1] = use_oy - PGoy * units::m / units::mm;
    dev_start[2] = use_oz - PGoz * units::m / units::mm;

    double use_theta = arr_theta[ it ];
    double use_phi   = arr_phi[ ip ];

    dev_sphere[0] = (use_theta - PGtheta) * 180.0 / constants::kPi; // deg
    dev_sphere[1] = (use_phi - PGphi) * 180.0 / constants::kPi;

    double use_cx = std::cos( use_phi ) * std::sin( use_theta );
    double use_cy = std::sin( use_phi ) * std::sin( use_theta );
    double use_cz = std::cos( use_theta );

    // now, we set the correct start point and momentum.
    p4HNL->SetPxPyPzE( p3HNL * use_cx, p3HNL * use_cy, p3HNL * use_cz, gOptEnergyHNL );

    TVector3 startPoint, momentum;
    TVector3 entryPoint, exitPoint, decayPoint;

    startPoint.SetXYZ( use_ox, use_oy, use_oz );
    momentum.SetXYZ( p4HNL->Px(), p4HNL->Py(), p4HNL->Pz() );

    use_start[0] = use_ox;
    use_start[1] = use_oy;
    use_start[2] = use_oz;

    use_momentum[0] = p4HNL->Px();
    use_momentum[1] = p4HNL->Py();
    use_momentum[2] = p4HNL->Pz();
    use_momentum[3] = p4HNL->E();

    LOG( "gevald_hnl", pDEBUG )
      << "Set start point for this trajectory = " << utils::print::Vec3AsString( &startPoint )
      << " [mm]";
    LOG( "gevald_hnl", pDEBUG )
      << "Set momentum for this trajectory = " << utils::print::Vec3AsString( &momentum )
      << " [GeV/c]";

    TLorentzVector tmpVtx( use_ox, use_oy, use_oz, 0.0);
    event->SetVertex( tmpVtx );
    TLorentzVector tmpMom = *p4HNL;
    GHepParticle ptHNL( genie::kPdgHNL, kIStInitialState, -1, -1, -1, -1, tmpMom, tmpVtx );
    event->AddParticle( ptHNL );
    LOG( "gevald_hnl", pDEBUG ) 
      << "\nProbe p4 = " << utils::print::P4AsString( event->Particle(0)->P4() );
    setenv( "PRODVTXDIR", "NODIR", 1 ); // needed to prevent hnlgen from crashing
    dkVol->SetStartingParameters( event, 1.0e+20, false, true, gOptRootGeom ); LOG( "gevald_hnl", pDEBUG ) << "Line 4";

    dkVol->ProcessEventRecord(event);

    if( event->Vertex()->T() != -999.9 ){
      dkVol->GetInterestingPoints( entryPoint, exitPoint, decayPoint );
      use_wgt = event->Weight();

      // set the branch entries now
      use_entry[0] = entryPoint.X();
      use_entry[1] = entryPoint.Y();
      use_entry[2] = entryPoint.Z();
      
      use_exit[0]  = exitPoint.X();
      use_exit[1]  = exitPoint.Y();
      use_exit[2]  = exitPoint.Z();
      
      use_decay[0] = decayPoint.X();
      use_decay[1] = decayPoint.Y();
      use_decay[2] = decayPoint.Z();
      
      // also fill the histos
      hWeight.Fill( std::log10( use_wgt ), 1.0 );
      didIntersectDet = true;

    } else { // didn't intersect vertex, use nonsense

      use_entry[0] = -999.9;
      use_entry[1] = -999.9;
      use_entry[2] = -999.9;

      use_exit[0] = -999.9;
      use_exit[1] = -999.9;
      use_exit[2] = -999.9;

      use_decay[0] = -999.9;
      use_decay[1] = -999.9;
      use_decay[2] = -999.9;

      didIntersectDet = false;
    }
    outTree->Fill();
    ievent++;
  }

  // save to file and exit
  
  fout->cd();
  outTree->Write();
  hWeight.Write();
  hLength.Write();
  top_volume->Write();
  fout->Close();
    
  return 0;
}
//_________________________________________________________________________________________
void InitBoundingBox(void)
{
// Initialise geometry bounding box, used for generating HNL vertex positions

  LOG("gevald_hnl", pINFO)
    << "Initialising geometry bounding box.";

  fdx = 0; // half-length - x
  fdy = 0; // half-length - y
  fdz = 0; // half-length - z
  fox = 0; // origin - x
  foy = 0; // origin - y
  foz = 0; // origin - z

  if(!gOptUsingRootGeom) return;

  bool geom_is_accessible = ! (gSystem->AccessPathName(gOptRootGeom.c_str()));
  if (!geom_is_accessible) {
    LOG("gevald_hnl", pFATAL)
      << "The specified ROOT geometry doesn't exist! Initialization failed!";
    exit(1);
  }

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 

  // RETHERE implement top volume option from cmd line
  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );
  TGeoShape * ts  = top_volume->GetShape();

  TGeoBBox *  box = (TGeoBBox *)ts;
  
  // pass this box to HNLDecayVolume
  const Algorithm * algDkVol = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLDecayVolume", "Default");
  
  const HNLDecayVolume * dkVol = dynamic_cast< const HNLDecayVolume * >( algDkVol );
  dkVol->ImportBoundingBox( box );

  //get box origin and dimensions (in the same units as the geometry)
  fdx = box->GetDX();
  fdy = box->GetDY();
  fdz = box->GetDZ();
  fox = (box->GetOrigin())[0];
  foy = (box->GetOrigin())[1];
  foz = (box->GetOrigin())[2];

  LOG("gevald_hnl", pINFO)
    << "Before conversion the bounding box has:"
    << "\nOrigin = ( " << fox << " , " << foy << " , " << foz << " )"
    << "\nDimensions = " << fdx << " x " << fdy << " x " << fdz
    << "\n1cm = 1.0 unit";

  // Convert from local to SI units
  fdx *= gOptGeomLUnits;
  fdy *= gOptGeomLUnits;
  fdz *= gOptGeomLUnits;
  fox *= gOptGeomLUnits;
  foy *= gOptGeomLUnits;
  foz *= gOptGeomLUnits;

  LOG("gevald_hnl", pINFO)
    << "Initialised bounding box successfully.";

}
//_________________________________________________________________________________________
#endif // #ifdef __CAN_USE_ROOT_GEOM__
//............................................................................
//_________________________________________________________________________________________
const EventRecordVisitorI * HNLGenerator(void)
{
  //string sname   = "genie::EventGenerator";
  //string sconfig = "BeamHNL";
  AlgFactory * algf = AlgFactory::Instance();

  LOG("gevald_hnl", pINFO)
    << "Instantiating HNL generator.";

  const Algorithm * algmcgen = algf->GetAlgorithm(kDefOptSName, kDefOptSConfig);
  LOG("gevald_hnl", pDEBUG)
    << "Got algorithm " << kDefOptSName.c_str() << "/" << kDefOptSConfig.c_str();;

  const EventRecordVisitorI * mcgen = 
    dynamic_cast< const EventRecordVisitorI * >( algmcgen );
  if(!mcgen) {
     LOG("gevald_hnl", pFATAL) << "Couldn't instantiate the HNL generator";
     gAbortingInErr = true;
     exit(1);
  }

  LOG("gevald_hnl", pINFO)
    << "HNL generator instantiated successfully.";

  return mcgen;
}

//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevald_hnl", pINFO) << "Parsing command line arguments";

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
    LOG("gevald_hnl", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevald_hnl", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevald_hnl", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevald_hnl", pFATAL)
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  if( parser.OptionExists('M') ) {
    LOG("gevald_hnl", pDEBUG)
      << "Detecting mode. . .";
    gOptValidationMode = (HNLValidation_t) parser.ArgAsInt('M');
  } else {
    LOG("gevald_hnl", pFATAL)
      << "You must specify a validation mode.";
    PrintSyntax();
    exit(0);
  } // -M

  bool isMonoEnergeticFlux = true;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  if( parser.OptionExists('f') ) {
    LOG("gevald_hnl", pDEBUG)
      << "A flux has been offered. Searching this path: " << parser.ArgAsString('f');
    isMonoEnergeticFlux = false;
    gOptFluxFilePath = parser.ArgAsString('f');
    
    // check if this is dk2nu
    if( gOptFluxFilePath.find( "dk2nu" ) != string::npos ){
      gOptIsUsingDk2nu = true;
      LOG("gevald_hnl", pDEBUG)
	<< "dk2nu flux files detected. Will create flux spectrum dynamically.";
    }
  } else {
    // we need the 'E' option! Log it and pass below
    LOG("gevald_hnl", pINFO)
      << "No flux file offered. Assuming monoenergetic flux.";
  } //-f
  gOptIsMonoEnFlux = isMonoEnergeticFlux;
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  
#ifdef __CAN_USE_ROOT_GEOM__
  if( parser.OptionExists('g') ) {
    LOG("gevald_hnl", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');
    
    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file =
      ! (gSystem->AccessPathName(geom.c_str()));
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    } else {
      LOG("gevald_hnl", pFATAL)
	<< "Geometry option is not a ROOT file. This is a work in progress; please use ROOT geom.";
      PrintSyntax();
      exit(1);
    }
  } else if( gOptValidationMode == kValGeom ) {
    
    LOG("gevald_hnl", pFATAL)
      << "No geometry option specified - Exiting";
    PrintSyntax();
    exit(1);
  } //-g

  if( parser.OptionExists('L') ) {
    lunits = parser.ArgAsString('L');
    LOG("gevald_hnl", pDEBUG) << "Setting length units to " << lunits.c_str();
  } else {
    LOG("gevald_hnl", pDEBUG) << "Using default geometry length units";
    lunits = kDefOptGeomLUnits;
  } // -L
  gOptGeomLUnits = utils::units::UnitFromString(lunits);

  if( parser.OptionExists('A') ) {
    aunits = parser.ArgAsString('A');
    LOG("gevald_hnl", pDEBUG) << "Setting angle units to " << aunits.c_str();
  } else {
    LOG("gevald_hnl", pDEBUG) << "Using default angle length units";
    aunits = kDefOptGeomAUnits;
  } // -A
  gOptGeomAUnits = utils::units::UnitFromString(aunits);

  if( parser.OptionExists('T') ) {
    tunits = parser.ArgAsString('T');
    LOG("gevald_hnl", pDEBUG) << "Setting time units to " << tunits.c_str();
  } else {
    LOG("gevald_hnl", pDEBUG) << "Using default geometry time units";
    tunits = kDefOptGeomTUnits;
  } // -T
  gOptGeomTUnits = utils::units::UnitFromString(tunits);

#endif // #ifdef __CAN_USE_ROOT_GEOM__

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevald_hnl", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevald_hnl", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevald_hnl", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevald_hnl", pINFO) << "Unspecified random number seed - Using default";
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

  LOG("gevald_hnl", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevald_hnl job configuration");

  LOG("gevald_hnl", pNOTICE)
     << "\n @@ Run number    : " << gOptRunNu
     << "\n @@ Random seed   : " << gOptRanSeed
     << "\n @@ HNL mass      : " << gCfgMassHNL << " GeV"
     << "\n @@ Decay channel : " << utils::hnl::AsString(gCfgDecayMode)
     << "\n @@ Flux path     : " << gOptFluxFilePath
     << "\n @@ Geometry      : " << gminfo.str()
     << "\n @@ Statistics    : " << gOptNev << " events";
}
//_________________________________________________________________________________________
void ReadInConfig(void)
{
  LOG("gevald_hnl", pFATAL)
    << "Reading in validation configuration. . .";

  const Algorithm * algHNLGen = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLDecayer", "Default");
  const HNLDecayer * hnlgen = dynamic_cast< const HNLDecayer * >( algHNLGen );

  SimpleHNL confsh = hnlgen->GetHNLInstance( "BeamHNL" );
  gCfgMassHNL   = confsh.GetMass();
  const std::vector< double > confCoups = confsh.GetCouplings();
  gCfgECoupling = confCoups.at(0);
  gCfgMCoupling = confCoups.at(1);
  gCfgTCoupling = confCoups.at(2);

  gOptEnergyHNL = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-Energy" );

  gCfgHNLCx     = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-cx" );
  gCfgHNLCy     = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-cy" );
  gCfgHNLCz     = utils::hnl::GetCfgDouble( "HNL", "ParticleGun", "PG-cz" );

  double dircosMag2 = std::pow( gCfgHNLCx, 2.0 ) + 
    std::pow( gCfgHNLCy, 2.0 ) + 
    std::pow( gCfgHNLCz, 2.0 );
  double invdircosmag = 1.0 / std::sqrt( dircosMag2 );
  gCfgHNLCx *= invdircosmag;
  gCfgHNLCy *= invdircosmag;
  gCfgHNLCz *= invdircosmag;
  
  gCfgIntChannels = confsh.GetInterestingChannelsVec();

  const Algorithm * algFluxCreator = AlgFactory::Instance()->GetAlgorithm("genie::HNL::HNLFluxCreator", "Default");
  const HNLFluxCreator * fluxCreator = dynamic_cast< const HNLFluxCreator * >( algFluxCreator );

  std::vector< double > UserT = fluxCreator->GetB2UTranslation();
  gCfgUserOx    = UserT.at(0);
  gCfgUserOy    = UserT.at(1);
  gCfgUserOz    = UserT.at(2);

  std::vector< double > UserR = fluxCreator->GetB2URotation();
  gCfgUserAx1   = UserR.at(0);
  gCfgUserAz    = UserR.at(1);
  gCfgUserAx2   = UserR.at(2);

  // now transform the lengths and angles to the correct units
  gCfgUserOx   *= units::m / gOptGeomLUnits;
  gCfgUserOy   *= units::m / gOptGeomLUnits;
  gCfgUserOz   *= units::m / gOptGeomLUnits;

  gCfgUserAx1  *= units::rad / gOptGeomAUnits;
  gCfgUserAz   *= units::rad / gOptGeomAUnits;
  gCfgUserAx2  *= units::rad / gOptGeomAUnits;

  ostringstream csts; 
  csts << "Read out the following config:"
       << "\n"
       << "\nHNL mass = " << gCfgMassHNL << " [GeV]"
       << "\n|U_e4|^2 = " << gCfgECoupling
       << "\n|U_m4|^2 = " << gCfgMCoupling
       << "\n|U_t4|^2 = " << gCfgTCoupling
       << "\n"
       << "\nInteresting decay channels:";
  for( std::vector< HNLDecayMode_t >::iterator chit = gCfgIntChannels.begin();
       chit != gCfgIntChannels.end(); ++chit ){ csts << "\n\t" << utils::hnl::AsString(*chit); }
  csts << "\n"
       << "\nUser origin in beam coordinates = ( " << gCfgUserOx
       << ", " << gCfgUserOy << ", " << gCfgUserOz << " ) [" << lunits.c_str() << "]"
       << "\nEuler extrinsic x-z-x rotation = ( " << gCfgUserAx1
       << ", " << gCfgUserAz << ", " << gCfgUserAx2 << " ) [" << aunits.c_str() << "]"
       << "\nHNL particle-gun directional cosines: ( " << gCfgHNLCx << ", " << gCfgHNLCy 
       << ", " << gCfgHNLCz << ") [ GeV / GeV ]";

  LOG("gevald_hnl", pDEBUG) << csts.str();
  
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevald_hnl", pFATAL)
   << "\n **Syntax**"
   << "\n gevald_hnl [-h] "
   << "\n            [-r run#]"
   << "\n             -n n_of_events"
   << "\n            [-f path_to_flux_files]"
   << "\n            [-g geometry_file]"
   << "\n             -M mode:"
   << "\n                1: Flux prediction from dk2nu files. Needs -f option"
   << "\n                2: Flux prediction from histograms.  Needs -f option"
   << "\n                3: HNL decay validation. Specify an origin point and 4-momentum"
   << "\n                   in the \"ParticleGun\" section in config. Needs -g option."
   << "\n                4: Custom geometry file validation.  Needs -g option"
   << "\n                   Specify origin, momentum, and wiggle room for both of these in the"
   << "\n                   \"ParticleGun\" section in config"
   << "\n                   Regardless of how many events you ask for, this will evaluate 125x81"
   << "\n                   events: 5^3 from wiggling origin and 9^2 from wiggling momentum direction"
   << "\n               10: Full simulation (like gevgen_hnl but with lots of debug!)"
   << "\n"
   << "\n The configuration file lives at $GENIE/config/CommonHNL.xml - see"
   << " <param_set name=\"Validation\">"
   << "\n"
   << "\n Please also read the detailed documentation at http://www.genie-mc.org"
   << "\n or look at the source code: $GENIE/src/Apps/gBeamHNLValidationApp.cxx"
   << "\n";
}
//_________________________________________________________________________________________
