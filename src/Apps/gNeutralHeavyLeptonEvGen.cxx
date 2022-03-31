//________________________________________________________________________________________
/*!

\program gevgen_nhl

\brief   A GENIE-based neutral heavy lepton event generation application.

         *** Synopsis :

         gevgen_nhl [-h]
                   [-r run#]
                    -n n_of_events
		    -f path/to/flux/files
                   [-E nhl_energy]
                    --mass nhl_mass
                   [-m decay_mode]
		   [-g geometry (ROOT file)]
                   [-L geometry_length_units]
                   [-t geometry_top_volume_name]
                   [-o output_event_file_prefix]
                   [--seed random_number_seed]
                   [--message-thresholds xml_file]
                   [--event-record-print-level level]
                   [--mc-job-status-refresh-rate  rate]

         *** Options :

           [] Denotes an optional argument

           -h
              Prints out the gevgen_ndcy syntax and exits.
           -r
              Specifies the MC run number [default: 1000].
           -n
              Specifies how many events to generate.
           --mass
              Specifies the NHL mass (in GeV)
           -m
              NHL decay mode ID:
           -f
              Input NHL flux.
           -g
              Input detector geometry.
              If a geometry is specified, NHL decay vertices will be distributed
              in the desired detector volume.
              Using this argument, you can pass a ROOT file containing your
              detector geometry description.
           -L
              Input geometry length units, eg 'm', 'cm', 'mm', ...
              [default: 'mm']
           -t
              Input 'top volume' for event generation.
              The option be used to force event generation in given sub-detector.
              [default: the 'master volume' of the input geometry]
              You can also use the -t option to switch generation on/off at
              multiple volumes as, for example, in:
              `-t +Vol1-Vol2+Vol3-Vol4',
              `-t "+Vol1 -Vol2 +Vol3 -Vol4"',
              `-t -Vol2-Vol4+Vol1+Vol3',
              `-t "-Vol2 -Vol4 +Vol1 +Vol3"'m
              where:
              "+Vol1" and "+Vol3" tells GENIE to `switch on'  Vol1 and Vol3, while
              "-Vol2" and "-Vol4" tells GENIE to `switch off' Vol2 and Vol4.
              If the very first character is a '+', GENIE will neglect all volumes
              except the ones explicitly turned on. Vice versa, if the very first
              character is a `-', GENIE will keep all volumes except the ones
              explicitly turned off (feature contributed by J.Holeczek).
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

\created February 11, 2020

\cpright Copyright (c) 2003-2022, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//_________________________________________________________________________________________
// TODO: Make a NuMI alternative to input fluxes!
//       Implement geometry read-in for InitBoundingBox()
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
#include "Physics/NeutralHeavyLepton/NHLFluxReader.h"
#include "Physics/NeutralHeavyLepton/NHLPrimaryVtxGenerator.h"
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
using namespace genie::NHL::NHLFluxReader;
using namespace genie::NHL::NHLenums;

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX__
#include "Tools/Flux/GCylindTH1Flux.h"
#include <TH1.h>
#endif // #ifdef __GENIE_FLUX_DRIVERS_ENABLED__

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#define __CAN_USE_ROOT_GEOM__
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>
#endif // #ifdef __GENIE_GEOM_DRIVERS_ENABLED__

// function prototypes
void  GetCommandLineArgs (int argc, char ** argv);
void  PrintSyntax        (void);

int   SelectDecayMode    (std::vector<NHLDecayMode_t> *intChannels, SimpleNHL sh);
const EventRecordVisitorI * NHLGenerator(void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
void     GenerateEventsUsingFlux (void);
GFluxI * TH1FluxDriver           (void);
int      DecideType              (TFile * spectrumFile);
#endif // #ifdef __GENIE_FLUX_DRIVERS_ENABLED__

TLorentzVector GeneratePosition( GHepRecord * event );
#ifdef __CAN_USE_ROOT_GEOM__
void  InitBoundingBox    (void);
#endif // #ifdef __CAN_USE_ROOT_GEOM__

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format
string          kDefOptEvFilePrefix = "gntp";
string          kDefOptFluxFilePath = "./input-flux.root";

//
Long_t           gOptRunNu        = 1000;                // run number
int              gOptNev          = 10;                  // number of events to generate

double           gOptEnergyNHL    = -1;                  // NHL energy
double           gOptMassNHL      = -1;                  // NHL mass
double           gOptECoupling    = -1;                  // |U_e4|^2
double           gOptMCoupling    = -1;                  // |U_m4|^2
double           gOptTCoupling    = -1;                  // |U_t4|^2

bool             gOptIsMajorana   = false;               // is this Majorana?
int              gOptNHLKind      = -1;                  // 0 = nu, 1 = nubar, 2 = mix

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
string           gOptFluxFilePath = kDefOptFluxFilePath; // where flux files live
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
bool             gOptIsMonoEnFlux = true;                // do we have monoenergetic flux?

NHLDecayMode_t   gOptDecayMode    = kNHLDcyNull;         // NHL decay mode

string           gOptEvFilePrefix = kDefOptEvFilePrefix; // event file prefix
bool             gOptUsingRootGeom = false;              // using root geom or target mix?
string           gOptRootGeom;                           // input ROOT file with realistic detector geometry

string           gOptRootGeomTopVol = "";                // input geometry top event generation volume
double           gOptGeomLUnits = 0;                     // input geometry length units
long int         gOptRanSeed = -1;                       // random number seed

// Geometry bounding box and origin - read from the input geometry file (if any)
double fdx = 0; // half-length - x
double fdy = 0; // half-length - y
double fdz = 0; // half-length - z
double fox = 0; // origin - x
double foy = 0; // origin - y
double foz = 0; // origin - z

// HNL lifetime in rest frame
double CoMLifetime = -1.0; // GeV^{-1}
// an array to keep production vertex
double evProdVtx[3] = {0.0, 0.0, 0.0}; // mm
// event weight
double evWeight = 1.0;

//_________________________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  // Init messenger and random number seed
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);

  RandomGen * rnd = RandomGen::Instance();

  // RETHERE do some initial configuration re. couplings in job
  // same for kind and Majorana!
  gOptECoupling = 1.0e-4;
  gOptMCoupling = 1.0e-4;
  gOptTCoupling = 0.0;
  gOptNHLKind = 2; // for mixing
  gOptIsMajorana = false;

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu, gOptRanSeed);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  LOG("gevgen_nhl", pNOTICE)
    << "Initialised Ntuple Writer";

  // add another few branches to tree.
  ntpw.EventTree()->Branch("nhl_mass", &gOptMassNHL, "gOptMassNHL/D");
  ntpw.EventTree()->Branch("nhl_coup_e", &gOptECoupling, "gOptECoupling/D");
  ntpw.EventTree()->Branch("nhl_coup_m", &gOptMCoupling, "gOptMCoupling/D");
  ntpw.EventTree()->Branch("nhl_coup_t", &gOptTCoupling, "gOptTCoupling/D");
  ntpw.EventTree()->Branch("nhl_ismaj", &gOptIsMajorana, "gOptIsMajorana/I");
  ntpw.EventTree()->Branch("nhl_type", &gOptNHLKind, "gOptNHLKind/I");

  // Create a MC job monitor for a periodically updated status file
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  LOG("gevgen_nhl", pNOTICE)
    << "Initialised MC job monitor";

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

#ifdef __CAN_USE_ROOT_GEOM__
  // Read geometry bounding box - for vertex position generation
  if( gOptUsingRootGeom ){
    InitBoundingBox();
  }
#endif // #ifdef __CAN_USE_ROOT_GEOM__

  // Get the nucleon decay generator
  const EventRecordVisitorI * mcgen = NHLGenerator();

  // RETHERE either seek out input flux or loop over some flux tuples
  // WIP
  GFluxI * ff = 0; // only use this if the flux is not monoenergetic!
  TH1D * spectrum = 0;
  TFile * spectrumFile = 0;
  if( !gOptIsMonoEnFlux ){
    ff = TH1FluxDriver();

    // ask the TH1FluxDriver which flux it used
    spectrumFile = TFile::Open("./input-flux.root", "READ");
    TDirectory * baseDir = spectrumFile->GetDirectory("");
    std::string fluxName = std::string( "spectrum" );
    assert( baseDir->GetListOfKeys()->Contains( fluxName.c_str() ) );
    spectrum = ( TH1D * ) baseDir->Get( fluxName.c_str() );
    assert( spectrum && spectrum != NULL );
  }

  // Event loop
  int ievent = 0;
  while (1)
  {
    if( gOptNev >= 10000 ){
      if( ievent % (gOptNev / 1000 ) == 0 ){
	int irat = ievent / ( gOptNev / 1000 );
	std::cerr << 0.1 * irat << " % " << " ( " << ievent
		  << " / " << gOptNev << " ) \r" << std::flush;
      }
    }

    if(ievent == gOptNev) break;
      
     LOG("gevgen_nhl", pNOTICE)
          << " *** Generating event............ " << ievent;

     if( !gOptIsMonoEnFlux ){
       LOG( "gevgen_nhl", pDEBUG )
	 << "Getting energy from flux...";
       gOptEnergyNHL = spectrum->GetRandom();
       unsigned int ien = 0;
       while( gOptEnergyNHL <= gOptMassNHL && ien < controls::kRjMaxIterations ){
	 gOptEnergyNHL = spectrum->GetRandom(); // to prevent binning throwing E <= M
	 ien++;
       }
     }
     assert( gOptEnergyNHL > gOptMassNHL );

     int hpdg = genie::kPdgNHL;
     int typeMod = 1;
     
     // if not Majorana, check if we should be doing mixing
     // if yes, get from fluxes
     // if not, enforce nu vs nubar
     if( !gOptIsMajorana ){
       switch( gOptNHLKind ){
       case 0: // always nu
	 break;
       case 1: // always nubar
	 typeMod = -1;
	 break;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
       case 2:
	 typeMod = DecideType( spectrumFile );
	 break;
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
       default:
	 typeMod = 1; // for now;
       }
       hpdg *= typeMod;
       LOG("gevgen_nhl", pDEBUG) << "typeMod = " << typeMod;
     }

     EventRecord * event = new EventRecord;
     // int target = SelectInitState();
     int decay  = (int) gOptDecayMode;

     if( gOptDecayMode == kNHLDcyNull ){ // select from available modes
       LOG("gevgen_nhl", pNOTICE)
	 << "No decay mode specified - sampling from all available modes.";

       LOG("gevgen_nhl", pDEBUG)
	 << "Couplings are: "
	 << "\n|U_e4|^2 = " << gOptECoupling
	 << "\n|U_m4|^2 = " << gOptMCoupling
	 << "\n|U_t4|^2 = " << gOptTCoupling;

       assert( gOptECoupling >= 0.0 && gOptMCoupling >= 0.0 && gOptTCoupling >= 0.0 );

       // RETHERE assuming all these NHL have K+- parent. This is wrong 
       // (but not very wrong for interesting masses)
       LOG("gevgen_nhl", pDEBUG)
	 << " Building SimpleNHL object ";
       SimpleNHL sh( "NHL", ievent, hpdg, genie::kPdgKP, 
				 gOptMassNHL, gOptECoupling, gOptMCoupling, gOptTCoupling, false );
     
       LOG("gevgen_nhl", pDEBUG)
	 << " Creating interesting channels vector ";
       std::vector< NHLDecayMode_t > * intChannels = new std::vector< NHLDecayMode_t >();
       intChannels->emplace_back( kNHLDcyPiE );
       intChannels->emplace_back( kNHLDcyPiMu );
       //intChannels->emplace_back( kNHLDcyPi0Nu );
       //intChannels->emplace_back( kNHLDcyNuNuNu );
       //intChannels->emplace_back( kNHLDcyNuMuMu );
       //intChannels->emplace_back( kNHLDcyNuEE );
       //intChannels->emplace_back( kNHLDcyNuMuE );
       //intChannels->emplace_back( kNHLDcyPiPi0E );
       //intChannels->emplace_back( kNHLDcyPiPi0Mu );
       //intChannels->emplace_back( kNHLDcyPi0Pi0Nu );

       decay = SelectDecayMode( intChannels, sh );
     }

     Interaction * interaction = Interaction::NHL(typeMod * genie::kPdgNHL, gOptEnergyNHL, decay);
     event->AttachSummary(interaction);

     LOG("gevgen_nhl", pDEBUG)
       << "Note decay mode is " << utils::nhl::AsString(gOptDecayMode);

     // Simulate decay
     mcgen->ProcessEventRecord(event);

     // Generate a position for the decay vertex
     // also currently handles the event weight
     TLorentzVector x4 = GeneratePosition( event );
     event->SetVertex(x4);
     event->SetWeight( evWeight );
     LOG("gevgen_nhl", pDEBUG) << "Weight = " << evWeight;

     LOG("gevgen_nhl", pINFO)
         << "Generated event: " << *event;

     // Add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     delete event;

     ievent++;
  } // event loop

  // Save the generated event tree & close the output file
  ntpw.Save();

  LOG("gevgen_nhl", pNOTICE) << "Done!";

  return 0;
}
//_________________________________________________________________________________________
//............................................................................
#ifdef __CAN_USE_ROOT_GEOM__
void InitBoundingBox(void)
{
// Initialise geometry bounding box, used for generating NHL vertex positions

  LOG("gevgen_nhl", pINFO)
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
    LOG("gevgen_nhl", pFATAL)
      << "The specified ROOT geometry doesn't exist! Initialization failed!";
    exit(1);
  }

  TGeoManager * gm = TGeoManager::Import(gOptRootGeom.c_str());
  TGeoVolume * top_volume = gm->GetTopVolume();
  

  // RETHERE I only make this for some plots. Must remove
  TGeoVolume * tracker = gm->FindVolumeFast("DetectorlvTracker");
  assert(tracker);
  gm->SetTopVolume( tracker );
  TGeoShape * ts  = tracker->GetShape();
  
  //TGeoShape * ts  = top_volume->GetShape();
  TGeoBBox *  box = (TGeoBBox *)ts;
  
  // pass this box to NHLDecayVolume
  NHLDecayVolume::ImportBoundingBox( box );

  //get box origin and dimensions (in the same units as the geometry)
  fdx = box->GetDX(); // half-length
  fdy = box->GetDY(); // half-length
  fdz = box->GetDZ(); // half-length
  fox = (box->GetOrigin())[0];
  foy = (box->GetOrigin())[1];
  foz = (box->GetOrigin())[2];

  LOG("gevgen_nhl", pINFO)
    << "Before conversion the bounding box has:"
    << "\nOrigin = ( " << fox << " , " << foy << " , " << foz << " )"
    << "\nDimensions = " << fdx << " x " << fdy << " x " << fdz;

  // Convert from local to SI units
  fdx *= gOptGeomLUnits;
  fdy *= gOptGeomLUnits;
  fdz *= gOptGeomLUnits;
  fox *= gOptGeomLUnits;
  foy *= gOptGeomLUnits;
  foz *= gOptGeomLUnits;

  LOG("gevgen_nhl", pINFO)
    << "Initialised bounding box successfully.";
}
#endif // #ifdef __CAN_USE_ROOT_GEOM__
//............................................................................
//_________________________________________________________________________________________
//............................................................................
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
GFluxI * TH1FluxDriver(void)
{
  //
  //
  flux::GCylindTH1Flux * flux = new flux::GCylindTH1Flux;
  TH1D * spectrum = 0;

  double emin = 0.0; // RETHERE need to make this configurable.
  double emax = 0.0+100.0; 

  // read in mass of NHL and decide which fluxes to use
  
  assert(gOptMassNHL > 0.0);

  // select mass point

  int closest_masspoint = selectMass( gOptMassNHL );

  LOG("gevgen_nhl", pDEBUG)
    << "Mass inserted: " << gOptMassNHL << " GeV ==> mass point " << closest_masspoint;
  LOG("gevgen_nhl", pDEBUG)
    << "Using fluxes in base path " << gOptFluxFilePath.c_str();
  
  selectFile( gOptFluxFilePath, gOptMassNHL );
  string finPath = NHLFluxReader::fPath; // is it good practice to keep this explicit?
  LOG("gevgen_nhl", pDEBUG)
    << "Looking for fluxes in " << finPath.c_str();
  assert( !gSystem->AccessPathName( finPath.c_str()) );

  // extract specified flux histogram from input root file

  string hFluxName = string( "hHNLFluxCenterAcc" );
  hFluxName.append( Form( "_%d", closest_masspoint ) );

  TH1F *hfluxAllMu    = getFluxHist1F( finPath, hFluxName, kNumu );
  TH1F *hfluxAllMubar = getFluxHist1F( finPath, hFluxName, kNumubar );
  TH1F *hfluxAllE     = getFluxHist1F( finPath, hFluxName, kNue );
  TH1F *hfluxAllEbar  = getFluxHist1F( finPath, hFluxName, kNuebar );

  assert(hfluxAllMu);
  assert(hfluxAllMubar);
  assert(hfluxAllE);
  assert(hfluxAllEbar);

  LOG("gevgen_nhl", pDEBUG)
    << "The histos have entries and max: "
    << "\nNumu:    " << hfluxAllMu->GetEntries() << " entries with max = " << hfluxAllMu->GetMaximum()
    << "\nNumubar: " << hfluxAllMubar->GetEntries() << " entries with max = " << hfluxAllMubar->GetMaximum()
    << "\nNue:     " << hfluxAllE->GetEntries() << " entries with max = " << hfluxAllE->GetMaximum()
    << "\nNuebar:  " << hfluxAllEbar->GetEntries() << " entries with max = " << hfluxAllEbar->GetMaximum();

  // let's build the mixed flux.
  
  TH1F * spectrumF = (TH1F*) hfluxAllMu->Clone(0);

  if( gOptECoupling == 0.0 ){ // no e coupling
    if( gOptIsMajorana || gOptNHLKind == 2 ){
      spectrumF->Add( hfluxAllMu, 1.0 );
      spectrumF->Add( hfluxAllMubar, 1.0 );
    }
    else if( gOptNHLKind == 0 ){
      spectrumF->Add( hfluxAllMu, 1.0 );
    }
    else{
      spectrumF->Add( hfluxAllMubar, 1.0 );
    }
  }
  else if( gOptMCoupling == 0.0 ){ // no mu coupling
    if( gOptIsMajorana || gOptNHLKind == 2 ){
      spectrumF->Add( hfluxAllE, 1.0 );
      spectrumF->Add( hfluxAllEbar, 1.0 );
    }
    else if( gOptNHLKind == 0 ){
      spectrumF->Add( hfluxAllE, 1.0 );
    }
    else{
      spectrumF->Add( hfluxAllEbar, 1.0 );
    }
  }
  else{ // add larger coupling as 1
    double ratio = gOptECoupling / gOptMCoupling;
    double rE = ( gOptECoupling > gOptMCoupling ) ? 1.0 : ratio;
    double rM = ( gOptMCoupling > gOptECoupling ) ? 1.0 : 1.0 / ratio;
    if( gOptIsMajorana || gOptNHLKind == 2 ){
      spectrumF->Add( hfluxAllMu, rM );
      spectrumF->Add( hfluxAllMubar, rM );
      spectrumF->Add( hfluxAllE, rE );
      spectrumF->Add( hfluxAllEbar, rE );
    }
    else if( gOptNHLKind == 0 ){
      spectrumF->Add( hfluxAllMu, rM );
      spectrumF->Add( hfluxAllE, rE );
    }
    else{
      spectrumF->Add( hfluxAllMubar, rM );
      spectrumF->Add( hfluxAllEbar, rE );
    }
  }

  LOG( "gevgen_nhl", pDEBUG )
    << "\n\n !!! ------------------------------------------------"
    << "\n !!! gOptECoupling, gOptMCoupling, gOptTCoupling = " << gOptECoupling << ", " << gOptMCoupling << ", " << gOptTCoupling
    << "\n !!! gOptNHLKind = " << gOptNHLKind
    << "\n !!! gOptIsMajorana = " << gOptIsMajorana
    << "\n !!! ------------------------------------------------"
    << "\n !!! Flux spectrum has ** " << spectrumF->GetEntries() << " ** entries"
    << "\n !!! Flux spectrum has ** " << spectrumF->GetMaximum() << " ** maximum"
    << "\n !!! ------------------------------------------------ \n";

  // copy into TH1D, *do not use the Copy() function!*
  const int nbins = spectrumF->GetNbinsX();
  spectrum = new TH1D( "s", "s", nbins, spectrumF->GetBinLowEdge(1), 
		       spectrumF->GetBinLowEdge(nbins) + spectrumF->GetBinWidth(nbins) );
  for( Int_t ib = 0; ib <= nbins; ib++ ){
    spectrum->SetBinContent( ib, spectrumF->GetBinContent(ib) );
  }
  
  spectrum->SetNameTitle("spectrum","NHL_flux");
  spectrum->SetDirectory(0);
  for(int ibin = 1; ibin <= hfluxAllMu->GetNbinsX(); ibin++) {
    if(hfluxAllMu->GetBinLowEdge(ibin) + hfluxAllMu->GetBinWidth(ibin) > emax ||
       hfluxAllMu->GetBinLowEdge(ibin) < emin) {
      spectrum->SetBinContent(ibin, 0);
    }
  } // do I want to kill the overflow / underflow bins? Why?
  
  LOG("gevgen_nhl", pINFO) << spectrum->GetEntries() << " entries in spectrum";

  // save input flux

  TFile f("./input-flux.root","RECREATE");
  spectrum->Write();

  // store integrals in histo if not Majorana and mixed flux
  // usual convention: bin 0+1 ==> numu etc
  if( !gOptIsMajorana && gOptNHLKind == 2 ){
    TH1D * hIntegrals = new TH1D( "hIntegrals", "hIntegrals", 4, 0.0, 1.0 );
    hIntegrals->SetBinContent( 1, hfluxAllMu->Integral() );
    hIntegrals->SetBinContent( 2, hfluxAllMubar->Integral() );
    hIntegrals->SetBinContent( 3, hfluxAllE->Integral() );
    hIntegrals->SetBinContent( 4, hfluxAllEbar->Integral() );

    hIntegrals->SetDirectory(0);
    hIntegrals->Write();

    LOG( "gevgen_nhl", pDEBUG )
      << "\n\nIntegrals asked for and stored. Here are their values by type:"
      << "\nNumu: " << hfluxAllMu->Integral()
      << "\nNumubar: " << hfluxAllMubar->Integral()
      << "\nNue: " << hfluxAllE->Integral()
      << "\nNuebar: " << hfluxAllEbar->Integral() << "\n\n";
  }

  f.Close();
  LOG("gevgen_nhl", pDEBUG) 
    << "Written spectrum to ./input-flux.root";

  // keep "beam" == SM-neutrino beam direction at z s.t. cos(theta_z) == 1
  // angular deviation of NHL (which is tiny, if assumption of collimated parents is made) made in main
  
  // Don't use GCylindTH1Flux's in-built methods - yet.

  TVector3 bdir (0.0,0.0,1.0);
  TVector3 bspot(0.0,0.0,1.0);

  flux->SetNuDirection      (bdir);
  flux->SetBeamSpot         (bspot);
  flux->SetTransverseRadius (-1);
  flux->AddEnergySpectrum   (genie::kPdgNHL, spectrum);

  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  LOG("gevgen_nhl", pDEBUG)
    << "Returning flux driver and exiting method.";
  return flux_driver;
}
//............................................................................
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
//_________________________________________________________________________________________
TLorentzVector GeneratePosition( GHepRecord * event )
{
  double weight = 1.0;

  if( gOptUsingRootGeom ){

    Interaction * interaction = event->Summary();

    LOG("gevgen_nhl", pDEBUG)
      << "Set startPoint and momentum";
  
    TVector3 startPoint, momentum, entryPoint, exitPoint;
  
    // sample production vertex
    const TLorentzVector * x4NHL = interaction->InitState().GetTgtP4( kRfLab );
    
    LOG("gevgen_nhl", pDEBUG)
      << "Detected vertex at ( " << x4NHL->Px() << ", " << x4NHL->Py() << ", " << x4NHL->Pz() << ")";
    startPoint.SetXYZ( x4NHL->Px(), x4NHL->Py(), x4NHL->Pz() );
    LOG( "gevgen_nhl", pDEBUG )
      << "Set start point for this trajectory = ( " << startPoint.X() << ", " << startPoint.Y() << ", " << startPoint.Z() << " ) [cm]";

    evProdVtx[0] = 10.0 * x4NHL->Px();
    evProdVtx[1] = 10.0 * x4NHL->Py();
    evProdVtx[2] = 10.0 * x4NHL->Pz();
  
    // get momentum of this channel
    const TLorentzVector * p4NHL = interaction->InitState().GetProbeP4( kRfLab );
    momentum.SetXYZ( p4NHL->Px() / p4NHL->P(), p4NHL->Py() / p4NHL->P(), p4NHL->Pz() / p4NHL->P() );
    LOG( "gevgen_nhl", pDEBUG )
      << "Set momentum for trajectory = ( " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << " ) [GeV]";
  
    TGeoManager * gm = TGeoManager::Import(gOptRootGeom.c_str());
    // RETHERE I only make this for some plots. Must remove
    TGeoVolume * vol = gm->FindVolumeFast("DetectorlvTracker");
    NHLPrimaryVtxGenerator * nhlgen = new NHLPrimaryVtxGenerator();;
    
    int trajIdx = 0;
    bool didIntersectDet = NHLDecayVolume::VolumeEntryAndExitPoints( startPoint, momentum, entryPoint, exitPoint, gm, vol );
    while( !didIntersectDet && trajIdx < 1e+2 ){
      // sample prod vtx and momentum... again
      LOG( "gevgen_nhl", pDEBUG )
	<< "Sampling another trajectory (index = " << trajIdx << ")";
      std::vector< double > * newProdVtx  = nhlgen->GenerateDecayPosition( event );
      //std::vector< double > * newMomentum = nhlgen->GenerateMomentum( event );
      
      startPoint.SetXYZ( newProdVtx->at(0), newProdVtx->at(1), newProdVtx->at(2) );
      LOG( "gevgen_nhl", pDEBUG )
	<< "Set start point for this trajectory = ( " << startPoint.X() << ", " << startPoint.Y() << ", " << startPoint.Z() << " ) [cm]";

      // update the production vertex
      evProdVtx[0] = 10.0 * newProdVtx->at(0);
      evProdVtx[1] = 10.0 * newProdVtx->at(1);
      evProdVtx[2] = 10.0 * newProdVtx->at(2);
      
      trajIdx++;
      didIntersectDet = NHLDecayVolume::VolumeEntryAndExitPoints( startPoint, momentum, entryPoint, exitPoint, gm, vol );
    }
    
    if( trajIdx == 1e+2 && !didIntersectDet ){
      LOG( "gevgen_nhl", pFATAL )
	<< "Unable to make a single good trajectory that intersects the detector after 100 tries! Exiting...";
      exit(1);
    }

    // make sure we have a consistent unit system
    NHLDecayVolume::EnforceUnits( "mm", "rad", "ns" );

    // move CoMLifetime to ns from GeV^{-1}
    CoMLifetime *= 1.0 / ( units::ns * units::GeV );
    LOG( "gevgen_nhl", pDEBUG )
      << "CoMLifetime = " << CoMLifetime << " [ns]";

    double maxDx = exitPoint.X() - entryPoint.X();
    double maxDy = exitPoint.Y() - entryPoint.Y();
    double maxDz = exitPoint.Z() - entryPoint.Z();

    LOG( "gevgen_nhl", pDEBUG )
      << "maxDx, maxDy, maxDz = " << maxDx << ", " << maxDy << ", " << maxDz;

    double maxLength = std::sqrt( std::pow( maxDx , 2.0 ) +
				  std::pow( maxDy , 2.0 ) +
				  std::pow( maxDz , 2.0 ) );

    double betaMag = p4NHL->P() / p4NHL->E();
    double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );

    double elapsed_length = NHLDecayVolume::CalcTravelLength( betaMag, CoMLifetime, maxLength );
    double ratio_length = elapsed_length / maxLength;

    // from these we can also make the weight. It's P( survival ) * P( decay in detector | survival )
    
    double distanceBeforeDet = std::sqrt( std::pow( (entryPoint.X() - startPoint.X()), 2.0 ) + 
					  std::pow( (entryPoint.Y() - startPoint.Y()), 2.0 ) + 
					  std::pow( (entryPoint.Y() - startPoint.Z()), 2.0 ) ); // mm

    double timeBeforeDet = distanceBeforeDet / NHLDecayVolume::kNewSpeedOfLight; // ns lab
    double timeInsideDet = maxLength / NHLDecayVolume::kNewSpeedOfLight; // ns lab
    
    double LabToRestTime = 1.0 / ( betaMag * gamma );
    timeBeforeDet *= LabToRestTime; // ns rest
    timeInsideDet *= LabToRestTime; // ns rest

    weight = std::exp( - timeBeforeDet / CoMLifetime );
    double survProb = weight;
    weight *= ( 1.0 - std::exp( - timeInsideDet / CoMLifetime ) );
    double decayProb = weight / survProb;

    LOG( "gevgen_nhl", pDEBUG )
      << "Decay probability with betaMag, gamma, CoMLifetime, distanceBeforeDet, maxLength = "
      << betaMag << ", " << gamma << ", " << CoMLifetime << ", " << distanceBeforeDet << ", "
      << maxLength << " [mm, ns, ns^{-1}, mm/ns] "
      << "\nand start, entry, exit vertices = ( " << evProdVtx[0] << ", " << evProdVtx[1] << ", " << evProdVtx[2] << " ) , ( " << entryPoint.X() << ", " << entryPoint.Y() << ", " << entryPoint.Z() << " ) , ( " << exitPoint.X() << ", " << exitPoint.Y() << " , " << exitPoint.Z() << " ) [mm] :"
      << "\nLabToRestTime, timeBeforeDet, timeInsideDet = " << LabToRestTime << ", "
      << timeBeforeDet << ", " << timeInsideDet
      << "\nsurvProb, decayProb, weight = " << survProb << ", " << decayProb << ", " << weight;
    evWeight = weight;

    TVector3 decayPoint = NHLDecayVolume::GetDecayPoint( elapsed_length, entryPoint, momentum );

    // RETHERE move to decay vertex!!!
    TLorentzVector x4( decayPoint.X(), decayPoint.Y(), decayPoint.Z(), 0.0 );
    
    return x4;
  }
  else{
    RandomGen * rnd = RandomGen::Instance();
    TRandom3 & rnd_geo = rnd->RndGeom();
    
    double rndx = 2 * (rnd_geo.Rndm() - 0.5); // [-1,1]
    double rndy = 2 * (rnd_geo.Rndm() - 0.5); // [-1,1]
    double rndz = 2 * (rnd_geo.Rndm() - 0.5); // [-1,1]
    
    double t_gen = 0;
    double x_gen = fox + rndx * fdx;
    double y_gen = foy + rndy * fdy;
    double z_gen = foz + rndz * fdz;
    
    TLorentzVector pos(x_gen, y_gen, z_gen, t_gen);
    return pos;
  }
  
  return TLorentzVector( 0, 0, 0, 0);
}
//_________________________________________________________________________________________
const EventRecordVisitorI * NHLGenerator(void)
{
  string sname   = "genie::EventGenerator";
  string sconfig = "NeutralHeavyLepton";
  AlgFactory * algf = AlgFactory::Instance();

  LOG("gevgen_nhl", pINFO)
    << "Instantiating NHL generator.";

  const EventRecordVisitorI * mcgen =
     dynamic_cast<const EventRecordVisitorI *> (algf->GetAlgorithm(sname,sconfig));
  if(!mcgen) {
     LOG("gevgen_nhl", pFATAL) << "Couldn't instantiate the NHL generator";
     gAbortingInErr = true;
     exit(1);
  }

  LOG("gevgen_nhl", pINFO)
    << "NHL generator instantiated successfully.";

  return mcgen;
}
//_________________________________________________________________________________________
int SelectDecayMode( std::vector< NHLDecayMode_t > * intChannels, SimpleNHL sh )
{
  LOG("gevgen_nhl", pDEBUG)
    << " Getting valid channels ";
  const std::map< NHLDecayMode_t, double > gammaMap = sh.GetValidChannels();

  // set CoM lifetime now if unset
  if( CoMLifetime < 0.0 ){
    CoMLifetime = sh.GetCoMLifetime();
    LOG( "gevgen_nhl", pDEBUG )
      << "Rest frame CoMLifetime = " << CoMLifetime << " [GeV^{-1}]";
  }

  for( std::vector< NHLDecayMode_t >::iterator it = intChannels->begin(); it != intChannels->end(); ++it ){
    NHLDecayMode_t mode = *it;
    auto mapG = gammaMap.find( mode );
    double theGamma = mapG->second;
    LOG("gevgen_nhl", pDEBUG)
      << "For mode " << utils::nhl::AsString( mode ) << " gamma = " << theGamma;
  }

  LOG("gevgen_nhl", pDEBUG)
    << " Setting interesting channels map ";
  std::map< NHLDecayMode_t, double > intMap =
    NHLSelector::SetInterestingChannels( (*intChannels), gammaMap );
     
  LOG("gevgen_nhl", pDEBUG)
    << " Telling SimpleNHL about interesting channels ";
  sh.SetInterestingChannels( intMap );

  // get probability that channels in intChannels will be selected
  LOG("gevgen_nhl", pDEBUG)
    << " Building probablilities of interesting channels ";
  std::map< NHLDecayMode_t, double > PMap = 
    NHLSelector::GetProbabilities( intMap );
     
  // now do a random throw
  RandomGen * rnd = RandomGen::Instance();
  double ranThrow = rnd->RndGen().Uniform( 0., 1. ); // NHL's fate is sealed.

  LOG("gevgen_nhl", pDEBUG)
    << "Random throw = " << ranThrow;

  NHLDecayMode_t selectedDecayChan =
    NHLSelector::SelectChannelInclusive( PMap, ranThrow );

  int decay = ( int ) selectedDecayChan;
  return decay;
}
//_________________________________________________________________________________________
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
/// based on the mixing of nu vs nubar in beam, return 1 (nu) or -1 (nubar)
int DecideType(TFile * spectrumFile){
  string intName = "hIntegrals";
  TDirectory * baseDir = spectrumFile->GetDirectory("");
  assert( baseDir->GetListOfKeys()->Contains( intName.c_str() ) );

  // 4 integrals, depending on co-produced lepton pdg. Group mu + e and mubar + ebar
  TH1D * hIntegrals = ( TH1D * ) baseDir->Get( intName.c_str() );
  double muInt    = hIntegrals->GetBinContent(1);
  double eInt     = hIntegrals->GetBinContent(3);
  double mubarInt = hIntegrals->GetBinContent(2);
  double ebarInt  = hIntegrals->GetBinContent(4);

  double nuInt    = muInt + eInt;
  double nubarInt = mubarInt + ebarInt;
  double totInt   = nuInt + nubarInt;

  RandomGen * rnd = RandomGen::Instance();
  double ranthrow = rnd->RndGen().Uniform(0.0, 1.0);

  int typeMod = ( ranthrow <= nuInt / totInt ) ? 1 : -1;
  return typeMod;
}
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_nhl", pINFO) << "Parsing command line arguments";

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
    LOG("gevgen_nhl", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_nhl", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_nhl", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_nhl", pFATAL)
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  // NHL mass
  gOptMassNHL = -1;
  if( parser.OptionExists("mass") ) {
    LOG("gevgen_nhl", pDEBUG)
        << "Reading NHL mass";
    gOptMassNHL = parser.ArgAsDouble("mass");
  } else {
    LOG("gevgen_nhl", pFATAL)
        << "You need to specify the NHL mass";
    PrintSyntax();
    exit(0);
  } //--mass
  PDGLibrary * pdglib = PDGLibrary::Instance();
  pdglib->AddNHL(gOptMassNHL);

  bool isMonoEnergeticFlux = true;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  if( parser.OptionExists('f') ) {
    LOG("gevgen_nhl", pDEBUG)
      << "A flux has been offered. Searching this path: " << parser.ArgAsString('f');
    isMonoEnergeticFlux = false;
    gOptFluxFilePath = parser.ArgAsString('f');
  } else {
    // we need the 'E' option! Log it and pass below
    LOG("gevgen_nhl", pINFO)
      << "No flux file offered. Assuming monoenergetic flux.";
  } //-f
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

  // NHL energy (only relevant if we do not have an input flux)
  gOptEnergyNHL = -1;
  if( isMonoEnergeticFlux ){
    if( parser.OptionExists('E') ) {
      LOG("gevgen_nhl", pDEBUG)
        << "Reading NHL energy";
      gOptEnergyNHL = parser.ArgAsDouble('E');
    } else {
      LOG("gevgen_nhl", pFATAL)
        << "You need to specify the NHL energy";
      PrintSyntax();
      exit(0);
    } //-E
    assert(gOptEnergyNHL > gOptMassNHL);
  }

  gOptIsMonoEnFlux = isMonoEnergeticFlux;

  // NHL decay mode
  int mode = -1;
  if( parser.OptionExists('m') ) {
    LOG("gevgen_nhl", pDEBUG)
        << "Reading NHL decay mode";
    mode = parser.ArgAsInt('m');
  } else {
    LOG("gevgen_nhl", pINFO)
        << "No decay mode specified - will sample from allowed decay modes";
  } //-m
  gOptDecayMode = (NHLDecayMode_t) mode;

  bool allowed = utils::nhl::IsKinematicallyAllowed(gOptDecayMode, gOptMassNHL);
  if(!allowed) {
    LOG("gevgen_nhl", pFATAL)
      << "Specified decay is not allowed kinematically for the given NHL mass";
    PrintSyntax();
    exit(0);
  }

  //
  // geometry
  //

  string geom = "";
  string lunits;
#ifdef __CAN_USE_ROOT_GEOM__
  // string dunits;
  if( parser.OptionExists('g') ) {
    LOG("gevgen_nhl", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');

    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file =
            ! (gSystem->AccessPathName(geom.c_str()));
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    } else {
      LOG("gevgen_nhl", pFATAL)
	<< "Geometry option is not a ROOT file. This is a work in progress; please use ROOT geom.";
      PrintSyntax();
      exit(1);
    }
  } else {
      // LOG("gevgen_nhl", pFATAL)
      //   << "No geometry option specified - Exiting";
      // PrintSyntax();
      // exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // length units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_nhl", pDEBUG)
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_nhl", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // // density units:
     // if( parser.OptionExists('D') ) {
     //    LOG("gevgen_nhl", pDEBUG)
     //       << "Checking for input geometry density units";
     //    dunits = parser.ArgAsString('D');
     // } else {
     //    LOG("gevgen_nhl", pDEBUG) << "Using default geometry density units";
     //    dunits = kDefOptGeomDUnits;
     // } // -D
     gOptGeomLUnits = utils::units::UnitFromString(lunits);
     // gOptGeomDUnits = utils::units::UnitFromString(dunits);

     // check whether an event generation volume name has been
     // specified -- default is the 'top volume'
     if( parser.OptionExists('t') ) {
        LOG("gevgen_nhl", pDEBUG) << "Checking for input volume name";
        gOptRootGeomTopVol = parser.ArgAsString('t');
     } else {
        LOG("gevgen_nhl", pDEBUG) << "Using the <master volume>";
     } // -t

  } // using root geom?
#endif // #ifdef __CAN_USE_ROOT_GEOM__

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_nhl", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_nhl", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_nhl", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_nhl", pINFO) << "Unspecified random number seed - Using default";
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

  LOG("gevgen_nhl", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_nhl job configuration");

  LOG("gevgen_nhl", pNOTICE)
     << "\n @@ Run number    : " << gOptRunNu
     << "\n @@ Random seed   : " << gOptRanSeed
     << "\n @@ NHL mass      : " << gOptMassNHL << " GeV"
     << "\n @@ Decay channel : " << utils::nhl::AsString(gOptDecayMode)
     << "\n @@ Geometry      : " << gminfo.str()
     << "\n @@ Statistics    : " << gOptNev << " events";
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_nhl", pFATAL)
   << "\n **Syntax**"
   << "\n gevgen_nhl [-h] "
   << "\n            [-r run#]"
   << "\n             -n n_of_events"
   << "\n             -f path/to/flux/files"
   << "\n            [-E nhl_energy]"
   << "\n             --mass nhl_mass"
   << "\n            [-m decay_mode]"
   << "\n            [-g geometry (ROOT file)]"
   << "\n            [-t top_volume_name_at_geom]"
   << "\n            [-L length_units_at_geom]"
   << "\n            [-o output_event_file_prefix]"
   << "\n            [--seed random_number_seed]"
   << "\n            [--message-thresholds xml_file]"
   << "\n            [--event-record-print-level level]"
   << "\n            [--mc-job-status-refresh-rate  rate]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << " or look at the source code: $GENIE/src/Apps/gNeutralHeavyLeptonEvGen.cxx"
   << "\n";
}
//_________________________________________________________________________________________
