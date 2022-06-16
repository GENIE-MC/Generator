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
		   [--firstEvent first event for dk2nu flux readin]
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
           -m
              NHL decay mode ID:
           -f
              Input NHL flux.
	   --firstEvent
	      If using dk2nu fluxes, start reading at this entry
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
// TODO: Implement top volume
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

// function prototypes
void  GetCommandLineArgs (int argc, char ** argv);
void  PrintSyntax        (void);

int   SelectDecayMode    (std::vector<NHLDecayMode_t> *intChannels, SimpleNHL sh);
const EventRecordVisitorI * NHLGenerator(void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
void     GenerateEventsUsingFlux (void);
GFluxI * TH1FluxDriver           (void);
int      DecideType              (TFile * spectrumFile);
int      InitialiseTupleFlux     (std::string finpath);
void     MakeNHLFromTuple        (int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf, std::string finpath, int run);
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

string          kDefOptSName   = "genie::EventGenerator";
string          kDefOptSConfig = "NeutralHeavyLepton";

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
map<string,string> gOptFluxShortNames;
bool             gOptIsUsingDk2nu = false;               // using flat dk2nu files?
int              gOptFirstEvent   = 0;                  // skip to this entry in dk2nu
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
bool             gOptIsMonoEnFlux = true;                // do we have monoenergetic flux?

NHLDecayMode_t   gOptDecayMode    = kNHLDcyNull;         // NHL decay mode
std::vector< NHLDecayMode_t > gOptIntChannels;           // decays to un-inhibit

string           gOptEvFilePrefix = kDefOptEvFilePrefix; // event file prefix
bool             gOptUsingRootGeom = false;              // using root geom or target mix?
string           gOptRootGeom;                           // input ROOT file with realistic detector geometry

#ifdef __CAN_USE_ROOT_GEOM__
TGeoManager *    gOptRootGeoManager = 0;                 // the workhorse geometry manager
TGeoVolume  *    gOptRootGeoVolume  = 0;
#endif // #ifdef __CAN_USE_ROOT_GEOM__

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

double NTP_IS_E = 0., NTP_IS_PX = 0., NTP_IS_PY = 0., NTP_IS_PZ = 0.;
double NTP_FS0_E = 0., NTP_FS0_PX = 0., NTP_FS0_PY = 0., NTP_FS0_PZ = 0.;
double NTP_FS1_E = 0., NTP_FS1_PX = 0., NTP_FS1_PY = 0., NTP_FS1_PZ = 0.;
double NTP_FS2_E = 0., NTP_FS2_PX = 0., NTP_FS2_PY = 0., NTP_FS2_PZ = 0.;
int NTP_FS0_PDG = 0, NTP_FS1_PDG = 0, NTP_FS2_PDG = 0;

NHLPrimaryVtxGenerator * nhlgen = 0;
// HNL lifetime in rest frame
double CoMLifetime = -1.0; // GeV^{-1}
// == Gamma( all valid channels ) / Gamma( all interesting channels )
double decayMod = 1.0;
// an array to keep production vertex
double evProdVtx[4] = {0.0, 0.0, 0.0, 0.0}; // x,y,z,t: mm, ns
// event weight
double evWeight = 1.0;

//_________________________________________________________________________________________
int main(int argc, char ** argv)
{
  // suppress ROOT Info messages
  gErrorIgnoreLevel = kWarning;

  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  // Init messenger and random number seed
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);

  __attribute__((unused)) RandomGen * rnd = RandomGen::Instance();

  // Get the NHL generator first to load config
  // config loaded upon instantiation of NHLGenerator algorithm 
  // ==> NHLPrimaryVtxGenerator::LoadConfig()
  const EventRecordVisitorI * mcgen = NHLGenerator();
  if( !nhlgen ){
    nhlgen = new NHLPrimaryVtxGenerator(); // do NOT remove this if( !nhlgen ), it causes a MASSIVE memleak if you do.
  }
  string confString = kDefOptSName + kDefOptSConfig;
  //const double confMass = nhlgen->GetNHLMass( confString );
  //const std::vector< double > confCoups = nhlgen->GetNHLCouplings( confString );

  SimpleNHL confsh = nhlgen->GetNHLInstance( confString );
  const double confMass = confsh.GetMass();
  const std::vector< double > confCoups = confsh.GetCouplings();
  const bool confIsMajorana = confsh.GetIsMajorana();
  const int confType = confsh.GetType();
  const double confAngDev = confsh.GetAngularDeviation();
  const std::vector< double > confT = confsh.GetBeam2UserTranslation();
  const std::vector< double > confR = confsh.GetBeam2UserRotation();
  const std::vector< std::vector< double > > confRM = confsh.GetBeam2UserRotationMatrix();
  const std::vector< NHLDecayMode_t > confIntChan = confsh.GetInterestingChannelsVec();

  LOG( "gevgen_nhl", pDEBUG )
    << "At app stage we see:"
    << "\nMass = " << confMass << " GeV"
    << "\nECoup = " << confCoups.at(0)
    << "\nMCoup = " << confCoups.at(1)
    << "\nTCoup = " << confCoups.at(2)
    << "\nIsMajorana = " << confIsMajorana
    << "\nType = " << confType
    << "\nAngular deviation = " << confAngDev << " deg"
    << "\nBEAM origin is USER ( " << confT.at(0) << ", " << confT.at(1) << ", " << confT.at(2) << " ) [m]"
    << "\nEuler angles (extrinsic x-z-x == 1-2-3) are ( " << confR.at(0) << ", " << confR.at(1) << "," << confR.at(2) << " )"
    << "\nRotation matrix RM = Rx(1) * Rz(2) * Rx(3) : RM * BEAM = USER"
    << "\n = ( " << (confRM.at(0)).at(0) << ", " << (confRM.at(0)).at(1) << ", " << (confRM.at(0)).at(2) << " )"
    << "\n = ( " << (confRM.at(1)).at(0) << ", " << (confRM.at(1)).at(1) << ", " << (confRM.at(1)).at(2) << " )"
    << "\n = ( " << (confRM.at(2)).at(0) << ", " << (confRM.at(2)).at(1) << ", " << (confRM.at(2)).at(2) << " )";

  gOptECoupling = confCoups.at(0);
  gOptMCoupling = confCoups.at(1);
  gOptTCoupling = confCoups.at(2);
  gOptNHLKind = confType; // for mixing
  gOptIsMajorana = confIsMajorana;

  gOptIntChannels = confIntChan;

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

  // let's make NHL-specific FS branches until we get gntpc sorted out
  ntpw.EventTree()->Branch("nhl_IS_E", &NTP_IS_E, "NTP_IS_E/D");
  ntpw.EventTree()->Branch("nhl_IS_PX", &NTP_IS_PX, "NTP_IS_PX/D");
  ntpw.EventTree()->Branch("nhl_IS_PY", &NTP_IS_PY, "NTP_IS_PY/D");
  ntpw.EventTree()->Branch("nhl_IS_PZ", &NTP_IS_PZ, "NTP_IS_PZ/D");
  ntpw.EventTree()->Branch("nhl_FS0_PDG", &NTP_FS0_PDG, "NTP_FS0_PDG/I");
  ntpw.EventTree()->Branch("nhl_FS0_E", &NTP_FS0_E, "NTP_FS0_E/D");
  ntpw.EventTree()->Branch("nhl_FS0_PX", &NTP_FS0_PX, "NTP_FS0_PX/D");
  ntpw.EventTree()->Branch("nhl_FS0_PY", &NTP_FS0_PY, "NTP_FS0_PY/D");
  ntpw.EventTree()->Branch("nhl_FS0_PZ", &NTP_FS0_PZ, "NTP_FS0_PZ/D");
  ntpw.EventTree()->Branch("nhl_FS1_PDG", &NTP_FS1_PDG, "NTP_FS1_PDG/I");
  ntpw.EventTree()->Branch("nhl_FS1_E", &NTP_FS1_E, "NTP_FS1_E/D");
  ntpw.EventTree()->Branch("nhl_FS1_PX", &NTP_FS1_PX, "NTP_FS1_PX/D");
  ntpw.EventTree()->Branch("nhl_FS1_PY", &NTP_FS1_PY, "NTP_FS1_PY/D");
  ntpw.EventTree()->Branch("nhl_FS1_PZ", &NTP_FS1_PZ, "NTP_FS1_PZ/D");
  ntpw.EventTree()->Branch("nhl_FS2_PDG", &NTP_FS2_PDG, "NTP_FS2_PDG/I");
  ntpw.EventTree()->Branch("nhl_FS2_E", &NTP_FS2_E, "NTP_FS2_E/D");
  ntpw.EventTree()->Branch("nhl_FS2_PX", &NTP_FS2_PX, "NTP_FS2_PX/D");
  ntpw.EventTree()->Branch("nhl_FS2_PY", &NTP_FS2_PY, "NTP_FS2_PY/D");
  ntpw.EventTree()->Branch("nhl_FS2_PZ", &NTP_FS2_PZ, "NTP_FS2_PZ/D");

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

  // RETHERE either seek out input flux or loop over some flux tuples
  // WIP
  __attribute__((unused)) GFluxI * ff = 0; // only use this if the flux is not monoenergetic!
  TH1D * spectrum = 0;
  TFile * spectrumFile = 0;
  if( !gOptIsMonoEnFlux ){
    if( !gOptIsUsingDk2nu ){ 
      ff = TH1FluxDriver();

      // read flux from file
      spectrumFile = TFile::Open("./input-flux.root", "READ");
      TDirectory * baseDir = spectrumFile->GetDirectory("");
      std::string fluxName = std::string( "spectrum" );
      assert( baseDir->GetListOfKeys()->Contains( fluxName.c_str() ) );
      spectrum = ( TH1D * ) baseDir->Get( fluxName.c_str() );
      assert( spectrum && spectrum != NULL );

    } else{ // RETHERE gotta make non-flat trees
      LOG( "gevgen_nhl", pWARN )
	<< "Using input flux files. These are *flat dk2nu-like ROOT trees, so far...*";

      int maxFluxEntries = InitialiseTupleFlux( gOptFluxFilePath );
      if( gOptNev > maxFluxEntries ){
	LOG( "gevgen_nhl", pWARN )
	  << "You have asked for " << gOptNev << " events, but only provided "
	  << maxFluxEntries << " flux entries. Truncating events to " << maxFluxEntries << ".";
	gOptNev = maxFluxEntries;
      }
    }
  }

  // Event loop
  int ievent = 0;
  flux::GNuMIFluxPassThroughInfo * gnmf = ( !gOptIsMonoEnFlux && gOptIsUsingDk2nu ) ? 
    new flux::GNuMIFluxPassThroughInfo() : 0;
  
  while (1)
  {
    if( gOptNev >= 10000 ){
      if( ievent % (gOptNev / 1000 ) == 0 ){
	int irat = ievent / ( gOptNev / 1000 );
	std::cerr << 0.1 * irat << " % " << " ( " << ievent
		  << " / " << gOptNev << " ) \r" << std::flush;
      }
    }

    if((ievent-gOptFirstEvent) == gOptNev) break;
      
     LOG("gevgen_nhl", pNOTICE)
          << " *** Generating event............ " << ievent;

     if( !gOptIsMonoEnFlux ){
       if( !gOptIsUsingDk2nu ){
	 LOG( "gevgen_nhl", pDEBUG )
	   << "Getting energy from flux...";
	 gOptEnergyNHL = spectrum->GetRandom();
	 unsigned int ien = 0;
	 while( gOptEnergyNHL <= gOptMassNHL && ien < controls::kRjMaxIterations ){
	   gOptEnergyNHL = spectrum->GetRandom(); // to prevent binning throwing E <= M
	   ien++;
	 }
       } else { // get a full NHL from flux tuples
	 LOG( "gevgen_nhl", pDEBUG )
	   << "Starting reading from event " << gOptFirstEvent;
	 while( ievent < gOptFirstEvent ){
	   ievent++;
	 }
	 LOG( "gevgen_nhl", pDEBUG )
	   << "Making NHL from tuple for event " << ievent;
	 MakeNHLFromTuple( ievent, gnmf, gOptFluxFilePath, gOptRunNu );
       }

       if( gnmf ){ 
	 TLorentzVector gnmfP4 = gnmf->fgP4User;
	 gOptEnergyNHL = (gnmf->fgP4User).E();
	 LOG( "gevgen_nhl", pDEBUG )
	   << "Got TLorentzVector from gnmf: " << utils::print::P4AsString(&gnmfP4);

	 LOG( "gevgen_nhl", pDEBUG )
	   << "\nIn user coordinates, the decay happened at " << utils::print::X4AsString( &(gnmf->fgX4User) ) << " [m]"
	   << "\nAnd the trajectory is pointing towards ( "
	   << (gnmfP4.Px()/gnmfP4.P()) << ", " << (gnmfP4.Py()/gnmfP4.P()) << ", " << (gnmfP4.Pz()/gnmfP4.P());

	 if( gOptEnergyNHL < 0.0 ){
	   ievent++;
	   continue; // hit nonsense & could not generate NHL
	 }
       }
     }
     assert( gOptEnergyNHL > gOptMassNHL );

     int hpdg = genie::kPdgNHL;
     int typeMod = 1;
     
     // if not Majorana, check if we should be doing mixing
     // if yes, get from fluxes
     // if not, enforce nu vs nubar
     if( !gOptIsMajorana && !gOptIsUsingDk2nu ){
       switch( gOptNHLKind ){
       case 0: // always nu
	 break;
       case 1: // always nubar
	 typeMod = -1;
	 break;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
       case 2:
	 if( !gOptIsMonoEnFlux ) typeMod = DecideType( spectrumFile );
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

     const std::map< NHLDecayMode_t, double > gammaMap = sh.GetValidChannels();
     LOG( "gevgen_nhl", pDEBUG )
       << "CoMLifetime = " << sh.GetCoMLifetime();
     CoMLifetime = sh.GetCoMLifetime();

     if( gOptDecayMode == kNHLDcyNull ){ // select from available modes
       LOG("gevgen_nhl", pNOTICE)
	 << "No decay mode specified - sampling from all available modes.";

       LOG("gevgen_nhl", pDEBUG)
	 << "Reading interesting channels vector from config";
       std::vector< NHLDecayMode_t > * intChannels = &gOptIntChannels;

       decay = SelectDecayMode( intChannels, sh );
     }

     Interaction * interaction = Interaction::NHL(typeMod * genie::kPdgNHL, gOptEnergyNHL, decay);

     if( event->Vertex() ){
       LOG( "gevgen_nhl", pDEBUG )
	 << "\nIS p4  = " << utils::print::P4AsString( interaction->InitStatePtr()->GetProbeP4() )
	 << "\nIS vtx = " << utils::print::X4AsString( event->Vertex() );
     } else {
       LOG( "gevgen_nhl", pDEBUG )
	 << "\nIS p4  = " << utils::print::P4AsString( interaction->InitStatePtr()->GetProbeP4() );
     }

     double acceptance = 1.0; // need to weight a spectrum by acceptance and nimpwt as well

     if( gnmf ){ // we have an NHL with definite momentum, so let's set it now
       interaction->InitStatePtr()->SetProbeP4( gnmf->fgP4User );
       LOG( "gevgen_nhl", pDEBUG )
	 << "\ngnmf->fgP4User setting probe p4 = " << utils::print::P4AsString( &gnmf->fgP4User );
       event->SetVertex( gnmf->fgX4User );
       acceptance = gnmf->nimpwt * gnmf->fgXYWgt;
     }

     event->AttachSummary(interaction);

     LOG("gevgen_nhl", pDEBUG)
       << "Note decay mode is " << utils::nhl::AsString(gOptDecayMode);

     // Simulate decay
     mcgen->ProcessEventRecord(event);

     // add the FS 4-momenta to special branches
     // Quite inelegant. Gets the job done, though
     NTP_FS0_PDG = (event->Particle(1))->Pdg();
     NTP_FS0_E  = ((event->Particle(1))->P4())->E();
     NTP_FS0_PX = ((event->Particle(1))->P4())->Px();
     NTP_FS0_PY = ((event->Particle(1))->P4())->Py();
     NTP_FS0_PZ = ((event->Particle(1))->P4())->Pz();
     NTP_FS1_PDG = (event->Particle(2))->Pdg();
     NTP_FS1_E  = ((event->Particle(2))->P4())->E();
     NTP_FS1_PX = ((event->Particle(2))->P4())->Px();
     NTP_FS1_PY = ((event->Particle(2))->P4())->Py();
     NTP_FS1_PZ = ((event->Particle(2))->P4())->Pz();
     if( event->Particle(3) ){
       NTP_FS2_PDG = (event->Particle(3))->Pdg();
       NTP_FS2_E  = ((event->Particle(3))->P4())->E();
       NTP_FS2_PX = ((event->Particle(3))->P4())->Px();
       NTP_FS2_PY = ((event->Particle(3))->P4())->Py();
       NTP_FS2_PZ = ((event->Particle(3))->P4())->Pz();
     }
     else{
       NTP_FS2_PDG = 0;
       NTP_FS2_E = 0.0;
       NTP_FS2_PX = 0.0;
       NTP_FS2_PY = 0.0;
       NTP_FS2_PZ = 0.0;
     }

     // Generate (or read) a position for the decay vertex
     // also currently handles the event weight
     TLorentzVector x4mm = GeneratePosition( event );

     const double mmtom = genie::units::mm / genie::units::m;
     TLorentzVector x4m( x4mm.X() * mmtom, x4mm.Y() * mmtom, x4mm.Z() * mmtom, evProdVtx[3] );
     event->SetVertex(x4m);
     // update weight to scale for couplings, acceptance, inhibited decays + geometry
     LOG( "gevgen_nhl", pDEBUG )
       << "\nWeight modifications:"
       << "\nCouplings^(-1) = " << 1.0 / ( gOptECoupling + gOptMCoupling + gOptTCoupling )
       << "\n(Acceptance * nimpwt)^(-1) = " << acceptance
       << "\nDecays^(-1) = " << decayMod
       << "\nGeometry^(-1) = " << evWeight;
     evWeight *= 1.0 / ( gOptECoupling + gOptMCoupling + gOptTCoupling );
     evWeight *= 1.0 / acceptance;
     evWeight *= 1.0 / decayMod;
     event->SetWeight( evWeight );

     // why does InitState show the wrong p4 here?
     interaction->InitStatePtr()->SetProbeP4( *(event->Particle(0)->P4()) );
     LOG( "gevgen_nhl", pDEBUG ) << 
       "\n!*!*!*! " << utils::print::P4AsString( event->Particle(0)->P4() );
     
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

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 

  // RETHERE implement top volume option from cmd line
  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );
  TGeoShape * ts  = top_volume->GetShape();

  TGeoBBox *  box = (TGeoBBox *)ts;
  
  // pass this box to NHLDecayVolume
  NHLDecayVolume::ImportBoundingBox( box );

  //get box origin and dimensions (in the same units as the geometry)
  fdx = box->GetDX();
  fdy = box->GetDY();
  fdz = box->GetDZ();
  fox = (box->GetOrigin())[0];
  foy = (box->GetOrigin())[1];
  foz = (box->GetOrigin())[2];

  LOG("gevgen_nhl", pINFO)
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

  double emin = 0.0; 
  double emax = utils::nhl::GetCfgDouble( "NHL", "InitialState", "NHL-max-energy" ); 

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
  string prodVtxPath = gOptFluxFilePath; prodVtxPath.append("/NHL_vertex_positions.root");
  __attribute__((unused)) int iset = setenv( "PRODVTXDIR", prodVtxPath.c_str(), 1 );
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
    else if( gOptNHLKind == 1 ){
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
//_________________________________________________________________________________________
int InitialiseTupleFlux( std::string finpath )
{
  LOG( "gevgen_nhl", pDEBUG )
    << "Opening input flux now from path " << finpath.c_str();

  NHLFluxCreator::OpenFluxInput( gOptFluxFilePath );
  //assert( NHLFluxCreator::tree && NHLFluxCreator::meta && NHLFluxCreator::tree->GetEntries() > 0 );
  //return NHLFluxCreator::tree->GetEntries();
  assert( NHLFluxCreator::ctree && NHLFluxCreator::cmeta && NHLFluxCreator::ctree->GetEntries() > 0 );
  return NHLFluxCreator::ctree->GetEntries();
}
//_________________________________________________________________________________________
void MakeNHLFromTuple( int iEntry, flux::GNuMIFluxPassThroughInfo * gnmf, std::string finpath, int run )
{
  // This genereates a full NHL from the flux tuples
  // by interfacing with NHLFluxCreator

  NHLFluxCreator::MakeTupleFluxEntry( iEntry, gnmf, finpath, run );
  
  LOG( "gevgen_nhl", pDEBUG ) << "MakeNHLFromTuple complete.";
}
//............................................................................
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
//_________________________________________________________________________________________
TLorentzVector GeneratePosition( GHepRecord * event )
{
  
  double weight = 1.0;
  double uMult = ( gOptIsUsingDk2nu ) ?
    units::m / units::mm : units::cm / units::mm;
  if( gOptUsingRootGeom ){

    __attribute__((unused)) Interaction * interaction = event->Summary();

    LOG("gevgen_nhl", pDEBUG)
      << "Set startPoint and momentum";
  
    TVector3 startPoint, momentum, entryPoint, exitPoint;
  
    // sample production vertex
    //const TLorentzVector * x4NHL = interaction->InitState().GetTgtP4( kRfLab );
    //if( !nhlgen ) nhlgen = new NHLPrimaryVtxGenerator();
    //const TLorentzVector * x4NHL = nhlgen->GetProdVtxPosition(event);
    TLorentzVector * x4NHL = event->Probe()->GetX4();
    
    std::ostringstream msts;
    if( gOptIsUsingDk2nu ){ msts << "[cm]"; }
    else{ msts << "[cm]"; }
    LOG("gevgen_nhl", pDEBUG)
      << "Detected vertex at ( " << x4NHL->X() << ", " << x4NHL->Y() << ", " << x4NHL->Z() << " )" << msts.str() << ": delay = " << x4NHL->T() << " [ns]";
    double xmult = ( gOptIsUsingDk2nu ) ? 10.0 : 1.0; // cm to mm in dk2nu case
    startPoint.SetXYZ( xmult * x4NHL->X(), xmult * x4NHL->Y(), xmult * x4NHL->Z() );

    evProdVtx[0] = uMult * x4NHL->X();
    evProdVtx[1] = uMult * x4NHL->Y();
    evProdVtx[2] = uMult * x4NHL->Z();
    evProdVtx[3] = x4NHL->T(); // ns

    LOG( "gevgen_nhl", pDEBUG )
      << "Set start point for this trajectory = ( " << startPoint.X() << ", " << startPoint.Y() << ", " << startPoint.Z() << " ) [mm]";
  
    // get momentum of this channel
    //const TLorentzVector * p4NHL = interaction->InitState().GetProbeP4( kRfLab );
    const TLorentzVector * p4NHL = event->Probe()->GetP4();
    
    NTP_IS_E = p4NHL->E(); NTP_IS_PX = p4NHL->Px(); NTP_IS_PY = p4NHL->Py(); NTP_IS_PZ = p4NHL->Pz();
    
    momentum.SetXYZ( p4NHL->Px() / p4NHL->P(), p4NHL->Py() / p4NHL->P(), p4NHL->Pz() / p4NHL->P() );
    LOG( "gevgen_nhl", pDEBUG )
      << "Set momentum for trajectory = ( " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << " ) [GeV]";
  
    if( !gOptRootGeoManager ){
      LOG("gevgen_nhl", pFATAL) << "Importing TGeoManager and doing top-volume business";
      gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str());
    }
    
    int trajIdx = 0; int trajMax = 20; // 1e+2;
    bool didIntersectDet = NHLDecayVolume::VolumeEntryAndExitPoints( startPoint, momentum, entryPoint, exitPoint, gOptRootGeoManager, gOptRootGeoVolume );

    if( gOptIsUsingDk2nu ) assert( didIntersectDet ); // forced to hit detector somewhere!

    std::vector< double > * newProdVtx = new std::vector< double >();
    newProdVtx->emplace_back( startPoint.X() );
    newProdVtx->emplace_back( startPoint.Y() );
    newProdVtx->emplace_back( startPoint.Z() );

    while( !didIntersectDet && trajIdx < trajMax ){
      // sample prod vtx and momentum... again
      LOG( "gevgen_nhl", pDEBUG )
	<< "Sampling another trajectory (index = " << trajIdx << ")";
      newProdVtx  = nhlgen->GenerateDecayPosition( event );
      //std::vector< double > * newMomentum = nhlgen->GenerateMomentum( event );
      
      // RETHERE FIX THIS TO USE MM.
      startPoint.SetXYZ( newProdVtx->at(0), newProdVtx->at(1), newProdVtx->at(2) );
      LOG( "gevgen_nhl", pDEBUG )
	<< "Set start point for this trajectory = ( " << startPoint.X() << ", " << startPoint.Y() << ", " << startPoint.Z() << " ) [cm]";

      // update the production vertex
      evProdVtx[0] = uMult * newProdVtx->at(0);
      evProdVtx[1] = uMult * newProdVtx->at(1);
      evProdVtx[2] = uMult * newProdVtx->at(2);
      
      trajIdx++;
      didIntersectDet = NHLDecayVolume::VolumeEntryAndExitPoints( startPoint, momentum, entryPoint, exitPoint, gOptRootGeoManager, gOptRootGeoVolume );

      newProdVtx->clear();
    }
    LOG("gevgen_nhl", pNOTICE) << "Called NHLDecayVolume::VolumeEntryAndExitPoints " << trajIdx + 1 << " times";
    
    if( trajIdx == trajMax && !didIntersectDet ){
      LOG( "gevgen_nhl", pERROR )
	<< "Unable to make a single good trajectory that intersects the detector after " << trajIdx << " tries! Bailing...";
      return *x4NHL;
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
      << "maxDx, maxDy, maxDz = " << maxDx << ", " << maxDy << ", " << maxDz << "[mm]";

    double maxLength = std::sqrt( std::pow( maxDx , 2.0 ) +
				  std::pow( maxDy , 2.0 ) +
				  std::pow( maxDz , 2.0 ) );

    double betaMag = p4NHL->P() / p4NHL->E();
    double gamma = std::sqrt( 1.0 / ( 1.0 - betaMag * betaMag ) );

    double elapsed_length = NHLDecayVolume::CalcTravelLength( betaMag, CoMLifetime, maxLength ); //mm
    __attribute__((unused)) double ratio_length = elapsed_length / maxLength;

    // from these we can also make the weight. It's P( survival ) * P( decay in detector | survival )
    double distanceBeforeDet = std::sqrt( std::pow( (entryPoint.X() - startPoint.X()), 2.0 ) + 
					  std::pow( (entryPoint.Y() - startPoint.Y()), 2.0 ) + 
					  std::pow( (entryPoint.Y() - startPoint.Z()), 2.0 ) ); // mm

    double timeBeforeDet = distanceBeforeDet / ( betaMag * NHLDecayVolume::kNewSpeedOfLight ); // ns lab
    double timeInsideDet = maxLength / ( betaMag * NHLDecayVolume::kNewSpeedOfLight ); // ns lab
    
    double LabToRestTime = 1.0 / ( gamma );
    timeBeforeDet *= LabToRestTime; // ns rest
    timeInsideDet *= LabToRestTime; // ns rest

    double survProb = std::exp( - timeBeforeDet / CoMLifetime );
    weight *= 1.0 / survProb;
    double decayProb = 1.0 - std::exp( - timeInsideDet / CoMLifetime );
    weight *= 1.0 / decayProb;

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

    TLorentzVector x4( decayPoint.X(), decayPoint.Y(), decayPoint.Z(), x4NHL->T() );

    delete x4NHL;
    delete p4NHL;

    newProdVtx->clear();
    
    return x4;
  }
  else{
    __attribute__((unused)) RandomGen * rnd = RandomGen::Instance();
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
  //string sname   = "genie::EventGenerator";
  //string sconfig = "NeutralHeavyLepton";
  AlgFactory * algf = AlgFactory::Instance();

  LOG("gevgen_nhl", pINFO)
    << "Instantiating NHL generator.";

  const Algorithm * algmcgen = algf->GetAlgorithm(kDefOptSName, kDefOptSConfig);
  LOG("gevgen_nhl", pDEBUG)
    << "Got algorithm " << kDefOptSName.c_str() << "/" << kDefOptSConfig.c_str();;

  const EventRecordVisitorI * mcgen = 
    dynamic_cast< const EventRecordVisitorI * >( algmcgen );
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

  // update fraction of total decay width that is not in inhibited channels
  double gammaAll = 0.0, gammaInt = 0.0;
  for( std::map< NHLDecayMode_t, double >::const_iterator itall = gammaMap.begin() ;
       itall != gammaMap.end() ; ++itall ){
    gammaAll += (*itall).second;
  }
  for( std::map< NHLDecayMode_t, double >::iterator itint = intMap.begin() ;
       itint != intMap.end() ; ++itint ){
    gammaInt += (*itint).second;
  }
  assert( gammaInt > 0.0 && gammaAll >= gammaInt );
  decayMod = gammaInt / gammaAll;
  

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

  /*
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
  //pdglib->AddNHL(gOptMassNHL);
  */

  // get NHL mass directly from config
  gOptMassNHL = genie::utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Mass" );
  // RETHERE check to see if mass is (not) given from config!

  bool isMonoEnergeticFlux = true;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  if( parser.OptionExists('f') ) {
    LOG("gevgen_nhl", pDEBUG)
      << "A flux has been offered. Searching this path: " << parser.ArgAsString('f');
    isMonoEnergeticFlux = false;
    gOptFluxFilePath = parser.ArgAsString('f');
    
    // check if this is dk2nu
    if( gOptFluxFilePath.find( "dk2nu" ) != string::npos ){
      gOptIsUsingDk2nu = true;
      LOG("gevgen_nhl", pDEBUG)
	<< "dk2nu flux files detected. Will create flux spectrum dynamically.";
    }
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

  // first flux entry to read
  if( parser.OptionExists("firstEvent") ) {
    gOptFirstEvent = parser.ArgAsInt("firstEvent");
    LOG( "gevgen_nhl", pINFO )
      << "Starting flux readin at first event = " << gOptFirstEvent;
  } // --firstEvent

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
   << "\n            [--firstEvent first_event_for_dk2nu_readin]"  
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
