//________________________________________________________________________________________
/*!

\program gevgen_hnl

\brief   A GENIE-based neutral heavy lepton event generation application.

         *** Synopsis :

         gevgen_hnl [-h]
                   [-r run#]
                    -n n_of_events
		    -f path/to/flux/files
                   [-E hnl_energy]
		   [--firstEvent first event for dk2nu flux readin]
                   [-m decay_mode]
		   [-g geometry (ROOT file)]
                   [-L geometry_length_units]
                   [-o output_event_file_prefix]
                   [--seed random_number_seed]
                   [--message-thresholds xml_file]
                   [--event-record-print-level level]
                   [--mc-job-status-refresh-rate  rate]

         *** Options :

           [] Denotes an optional argument

           -h
              Prints out the gevgen_hnl syntax and exits.
           -r
              Specifies the MC run number [default: 1000].
           -n
              Specifies how many events to generate.
           -m
              HNL decay mode ID:
           -f
              Input HNL flux.
	   --firstEvent
	      If using dk2nu fluxes, start reading at this entry
           -g
              Input detector geometry.
              If a geometry is specified, HNL decay vertices will be distributed
              in the desired detector volume.
              Using this argument, you can pass a ROOT file containing your
              detector geometry description.
           -L
              Input geometry length units, eg 'm', 'cm', 'mm', ...
              [default: 'mm']
           -o
              Sets the prefix of the output event file.
              The output filename is built as:
              [prefix].[run_number].[event_tree_format].[file_format]
              The default output filename is:
              gntp.[run_number].ghep.root
              This cmd line arguments lets you override 'gntp'
           --seed
              Random number seed.

\author  John Plows <komninos-john.plows \at cern.ch>
         University of Oxford

         Costas Andreopoulos <c.andreopoulos \at cern.ch>
         University of Liverpool

\created February 11, 2020

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
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
#include "Physics/BeamHNL/HNLVertexGenerator.h"
#include "Physics/BeamHNL/HNLFluxContainer.h"
#include "Physics/BeamHNL/HNLFluxCreator.h"
#include "Physics/BeamHNL/HNLDecayer.h"
#include "Physics/BeamHNL/SimpleHNL.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/TuneId.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/CmdLnArgParser.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;
using namespace genie::hnl;
using namespace genie::hnl::enums;

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX__
//#include "Tools/Flux/GNuMIFlux.h"
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
void   GetCommandLineArgs (int argc, char ** argv);
void   PrintSyntax        (void);

int    SelectDecayMode    (std::vector<HNLDecayMode_t> *intChannels, SimpleHNL sh);
const  EventRecordVisitorI * HNLGenerator(void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
void     GenerateEventsUsingFlux (void);
void     FillFluxNonsense        (FluxContainer &ggn);
void     FillFlux                (FluxContainer &ggn, FluxContainer &tgn);
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
string          kDefOptSConfig = "BeamHNL";
string          kDefOptSTune   = "GHNL20_00a_00_000";

//
Long_t           gOptRunNu        = 1000;                // run number
int              gOptNev          = 10;                  // number of events to generate

double           gOptEnergyHNL    = -1;                  // HNL energy
double           gOptMassHNL      = -1;                  // HNL mass
double           gOptECoupling    = -1;                  // |U_e4|^2
double           gOptMCoupling    = -1;                  // |U_m4|^2
double           gOptTCoupling    = -1;                  // |U_t4|^2

bool             gOptIsMajorana   = false;               // is this Majorana?

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
string           gOptFluxFilePath = kDefOptFluxFilePath; // where flux files live
map<string,string> gOptFluxShortNames;
bool             gOptIsUsingDk2nu = false;               // using flat dk2nu files?
int              gOptFirstEvent   = 0;                  // skip to this entry in dk2nu
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
bool             gOptIsMonoEnFlux = true;                // do we have monoenergetic flux?

HNLDecayMode_t   gOptDecayMode    = kHNLDcyNull;         // HNL decay mode
std::vector< HNLDecayMode_t > gOptIntChannels;           // decays to un-inhibit

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

// HNL lifetime in rest frame
double CoMLifetime = -1.0; // GeV^{-1}
// == Gamma( all valid channels ) / Gamma( all interesting channels )
double decayMod = 1.0;
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

  // Get the HNL generator first to load config
  // config loaded upon instantiation of HNLGenerator algorithm 
  // calls LoadConfig() of each sub-alg
  __attribute__((unused)) const EventRecordVisitorI * mcgen = HNLGenerator();
  const Algorithm * algFluxCreator = AlgFactory::Instance()->GetAlgorithm("genie::hnl::FluxCreator", "Default");
  const Algorithm * algHNLGen = AlgFactory::Instance()->GetAlgorithm("genie::hnl::Decayer", "Default");
  const Algorithm * algVtxGen = AlgFactory::Instance()->GetAlgorithm("genie::hnl::VertexGenerator", "Default");

  const FluxCreator * fluxCreator = dynamic_cast< const FluxCreator * >( algFluxCreator );
  const Decayer * hnlgen = dynamic_cast< const Decayer * >( algHNLGen );
  const VertexGenerator * vtxGen = dynamic_cast< const VertexGenerator * >( algVtxGen );
  
  //string confString = kDefOptSName + "/" + kDefOptSConfig;
  string confString = kDefOptSConfig;

  CoMLifetime = hnlgen->GetHNLLifetime(); // GeV^{-1}

  std::vector< double > U4l2s = hnlgen->GetHNLCouplings();
  gOptECoupling  = U4l2s.at(0);
  gOptMCoupling  = U4l2s.at(1);
  gOptTCoupling  = U4l2s.at(2);
  gOptIsMajorana = hnlgen->IsHNLMajorana();

  std::string stIntChannels = hnlgen->GetHNLInterestingChannels(); int iChan = -1;
  if( gOptIntChannels.size() > 0 ) gOptIntChannels.clear();
  while( stIntChannels.size() > 0 ){ // read channels from right (lowest mass) to left (highest mass)
    iChan++;
    HNLDecayMode_t md = static_cast< HNLDecayMode_t >( iChan );
    std::string tmpSt = stIntChannels.substr( stIntChannels.size()-1, stIntChannels.size() );
    if( std::strcmp( tmpSt.c_str(), "1" ) == 0 )
      gOptIntChannels.emplace_back( md );

    stIntChannels.erase( stIntChannels.end()-1, stIntChannels.end() );
  }

  //gOptIntChannels = confIntChan;

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu, gOptRanSeed);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // add flux info to the tree
  FluxContainer gnmf;
  if( gOptIsUsingDk2nu ) {
    // fill the flux object with nonsense to start with
    FluxContainer * ptGnmf = new FluxContainer();
    gnmf = *ptGnmf;
    delete ptGnmf;
    FillFluxNonsense( gnmf );
    TBranch * flux = ntpw.EventTree()->Branch( "flux",
					       "genie::hnl::FluxContainer",
					       &gnmf, 32000, 1 );
    flux->SetAutoDelete(kFALSE);
  }

  LOG("gevgen_hnl", pNOTICE)
    << "Initialised Ntuple Writer";

  // add another few branches to tree.
  ntpw.EventTree()->Branch("hnl_mass", &gOptMassHNL, "gOptMassHNL/D");
  ntpw.EventTree()->Branch("hnl_coup_e", &gOptECoupling, "gOptECoupling/D");
  ntpw.EventTree()->Branch("hnl_coup_m", &gOptMCoupling, "gOptMCoupling/D");
  ntpw.EventTree()->Branch("hnl_coup_t", &gOptTCoupling, "gOptTCoupling/D");
  ntpw.EventTree()->Branch("hnl_ismaj", &gOptIsMajorana, "gOptIsMajorana/I");

  // Create a MC job monitor for a periodically updated status file
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  LOG("gevgen_hnl", pNOTICE)
    << "Initialised MC job monitor";

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

#ifdef __CAN_USE_ROOT_GEOM__
  // Read geometry bounding box - for vertex position generation
  if( gOptUsingRootGeom ){
    InitBoundingBox();
  }
#endif // #ifdef __CAN_USE_ROOT_GEOM__

  if( !gOptIsMonoEnFlux ){
    LOG( "gevgen_hnl", pWARN )
      << "Using input flux files. These are *flat dk2nu-like ROOT trees, so far...*";

  }
  // Event loop
  int iflux = (gOptFirstEvent < 0) ? 0 : gOptFirstEvent; int ievent = iflux;
  int maxFluxEntries = -1;
  fluxCreator->SetInputFluxPath( gOptFluxFilePath );
  fluxCreator->SetGeomFile( gOptRootGeom );
  fluxCreator->SetFirstFluxEntry( iflux );
  vtxGen->SetGeomFile( gOptRootGeom );
  
  bool tooManyEntries = false;
  while (1)
  {
    if( tooManyEntries ){
      if( gOptNev >= 10000 ){
	if( (ievent-gOptFirstEvent) % (gOptNev / 1000) == 0 ){
	  int irat = (iflux-gOptFirstEvent) / ( gOptNev / 1000 );
	  std::cerr << 0.1 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev  << " processed ) \r" << std::flush;
	}
      } else if( gOptNev >= 100 ) {
	if( (ievent-gOptFirstEvent) % (gOptNev / 10) == 0 ){
	  int irat = (iflux-gOptFirstEvent) / ( gOptNev / 10 );
	  std::cerr << 10.0 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev  << " processed ) \r" << std::flush;
	}
      }
    } else {
      if( gOptNev >= 10000 ){
	if( (ievent-gOptFirstEvent) % (gOptNev / 1000) == 0 ){
	  int irat = (ievent-gOptFirstEvent) / ( gOptNev / 1000 );
	  std::cerr << 0.1 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev <<  " processed ) \r" << std::flush;
	}
      } else if( gOptNev >= 100 ) {
	if( (ievent-gOptFirstEvent) % (gOptNev / 10) == 0 ){
	  int irat = (ievent-gOptFirstEvent) / ( gOptNev / 10 );
	  std::cerr << 10.0 * irat << " % " << " ( " << (iflux-gOptFirstEvent)
		    << " seen ), ( " << (ievent-gOptFirstEvent) << " / " << gOptNev  << " processed ) \r" << std::flush;
	}
      }
    }

    if( tooManyEntries && ((iflux-gOptFirstEvent) == gOptNev) ) break;
    else if( (ievent-gOptFirstEvent) == gOptNev ) break;
    
    if( ievent < gOptFirstEvent ){ ievent++; continue; }
    
    assert( ievent >= gOptFirstEvent && gOptFirstEvent >= 0 );
    
    LOG("gevgen_hnl", pNOTICE)
      << " *** Generating event............ " << (ievent-gOptFirstEvent);
    
    EventRecord * event = new EventRecord;
    event->SetWeight(1.0);
    event->SetProbability( CoMLifetime );
    event->SetXSec( iflux ); // will be overridden, use as handy container
    evWeight = 1.0;
    
    if( !gOptIsMonoEnFlux ){
       fluxCreator->ProcessEventRecord( event );

       // fluxCreator->ProcessEventRecord now tells us how many entries there are
       maxFluxEntries = fluxCreator->GetNFluxEntries();
       if( gOptNev > maxFluxEntries - gOptFirstEvent ){
	 LOG( "gevgen_hnl", pWARN )
	   << "You have asked for " << gOptNev << " events, but only provided "
	   << maxFluxEntries - gOptFirstEvent << " flux entries. Truncating events to " << maxFluxEntries - gOptFirstEvent << ".";
	 gOptNev = maxFluxEntries - gOptFirstEvent;
	 tooManyEntries = true;
       }
       if( iflux >= maxFluxEntries - 1 ){
	 LOG( "gevgen_hnl", pWARN )
	   << "Reached end of flux input (iflux = " << iflux << "), about to stop.";
	 break;
       }
       
       FluxContainer retGnmf = fluxCreator->RetrieveFluxInfo();
       FillFlux( gnmf, retGnmf );
       
       // check to see if this was nonsense
       if( ! event->Particle(0) ){ iflux++; delete event; continue; }
       
       gOptEnergyHNL = event->Particle(0)->GetP4()->E();
       iflux++;
     } else { // monoenergetic HNL. Add it with energy and momentum pointing on z axis
       
       assert( gOptEnergyHNL > gOptMassHNL );
       double HNLP = std::sqrt( gOptEnergyHNL*gOptEnergyHNL - gOptMassHNL*gOptMassHNL );
       TLorentzVector probeP4( 0.0, 0.0, HNLP, gOptEnergyHNL );
       TLorentzVector v4( 0.0, 0.0, 0.0, 0.0 );
       GHepParticle ptHNL( kPdgHNL, kIStInitialState, -1, -1, -1, -1, probeP4, v4 );
       event->AddParticle( ptHNL );

     }
     assert( gOptEnergyHNL > gOptMassHNL );
     
     int hpdg = genie::kPdgHNL;
     int typeMod = 1;
     if( gOptIsMajorana ) typeMod = ( rnd->RndGen().Uniform( 0.0, 1.0 ) >= 0.5 ) ? -1 : 1;
     else if( event->Particle(0)->Pdg() > 0 ) typeMod = 1;
     else if( event->Particle(0)->Pdg() < 0 ) typeMod = -1;

     // int target = SelectInitState();
     int decay  = (int) gOptDecayMode;
     
     assert( gOptECoupling >= 0.0 && gOptMCoupling >= 0.0 && gOptTCoupling >= 0.0 );
     
     // RETHERE assuming all these HNL have K+- parent. This is wrong 
     // (but not very wrong for interesting masses)
     SimpleHNL sh( "HNL", ievent, hpdg, genie::kPdgKP, 
		   gOptMassHNL, gOptECoupling, gOptMCoupling, gOptTCoupling, false );

     const std::map< HNLDecayMode_t, double > gammaMap = sh.GetValidChannels();
     CoMLifetime = sh.GetCoMLifetime();

     if( gOptDecayMode == kHNLDcyNull ){ // select from available modes
       LOG("gevgen_hnl", pNOTICE)
	 << "No decay mode specified - sampling from all available modes.";

       std::vector< HNLDecayMode_t > * intChannels = &gOptIntChannels;

       decay = SelectDecayMode( intChannels, sh );
     }

     Interaction * interaction = Interaction::HNL(typeMod * genie::kPdgHNL, gOptEnergyHNL, decay);

     if( event->Particle(0) ){ // we have an HNL with definite momentum, so let's set it now
       interaction->InitStatePtr()->SetProbeP4( *(event->Particle(0)->P4()) );
       interaction->InitStatePtr()->SetProbePdg( event->Particle(0)->Pdg() );
       LOG( "gevgen_hnl", pDEBUG )
	 << "\nsetting probe p4 = " << utils::print::P4AsString( event->Particle(0)->P4() );
     }

     event->AttachSummary(interaction);

     // Simulate decay
     hnlgen->ProcessEventRecord(event);

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
     // also currently handles the geometrical weight
     TLorentzVector x4mm;
     if( gOptUsingRootGeom ){
       TLorentzVector * p4HNL = event->Particle(0)->GetP4();
       NTP_IS_E = p4HNL->E(); NTP_IS_PX = p4HNL->Px(); NTP_IS_PY = p4HNL->Py(); NTP_IS_PZ = p4HNL->Pz();
       vtxGen->ProcessEventRecord( event );
       x4mm = *(event->Vertex());
     } else {
       x4mm = GeneratePosition( event );
     }

     // update weight to scale for couplings, inhibited decays
     // acceptance is already handled in FluxCreator
     // geometry handled in VertexGenerator
     evWeight = event->Weight();
     evWeight *= 1.0 / ( gOptECoupling + gOptMCoupling + gOptTCoupling );
     evWeight *= 1.0 / decayMod;
     event->SetWeight( evWeight / 1.0e+20 ); // in units of 1e+20 POT

     // store the decayMod in the second daughter's polarisation
     event->Particle(2)->SetPolarization( 1.0, decayMod );

     // why does InitState show the wrong p4 here?
     interaction->InitStatePtr()->SetProbeP4( *(event->Particle(0)->P4()) );
     
     LOG("gevgen_hnl", pDEBUG) << "Weight = " << evWeight;

     LOG("gevgen_hnl", pINFO)
         << "Generated event: " << *event;

     // Add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);

     delete event;

     ievent++;
  } // event loop

  // Save the generated event tree & close the output file
  ntpw.Save();

  LOG("gevgen_hnl", pNOTICE) << "Done!";

  return 0;
}
//_________________________________________________________________________________________
//............................................................................
#ifdef __CAN_USE_ROOT_GEOM__
void InitBoundingBox(void)
{
// Initialise geometry bounding box, used for generating HNL vertex positions

  LOG("gevgen_hnl", pINFO)
    << "Initialising geometry bounding box.";

  fdx = 0; // half-length - x
  fdy = 0; // half-length - y
  fdz = 0; // half-length - z
  fox = 0; // origin - x
  foy = 0; // origin - y
  foz = 0; // origin - z

  if(!gOptUsingRootGeom){ // make a unit-m sided box
    LOG("gevgen_hnl", pINFO)
      << "No geometry file input detected, making a unit-m side box volume.";

    TGeoManager * geom = new TGeoManager( "box1", "A simple box detector" );

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

    //--- make the top container volume
    //const double boxSideX = 2.5, boxSideY = 2.5, boxSideZ = 2.5; // m
    //const double bigBoxSide = 2.0 * std::max( boxSideX, std::max( boxSideY, boxSideZ ) ); // m
    //const double worldLen = 1.01 * bigBoxSide; // m

    TGeoVolume * topvol = geom->MakeBox( "TOP", Vacuum, 101.0, 101.0, 101.0 );
    geom->SetTopVolume( topvol );

    //--- make the detector box container
    TGeoVolume * boxvol = geom->MakeBox( "VOL", Vacuum, 100.5, 100.5, 100.5 );
    boxvol->SetVisibility(kFALSE);

    //--- origin is at centre of the box
    TGeoVolume * box = geom->MakeBox( "BOX", Al, 100.0, 100.0, 100.0 );
    //TGeoTranslation * tr0 = new TGeoTranslation( 0.0, 0.0, 0.0 );
    TGeoRotation * rot0 = new TGeoRotation( "rot0", 90.0, 0.0, 90.0, 90.0, 0.0, 0.0 );

    //--- add directly to top volume
    topvol->AddNode( box, 1, rot0 );
    
    gOptRootGeoManager = geom;

    return;
  } 

  bool geom_is_accessible = ! (gSystem->AccessPathName(gOptRootGeom.c_str()));
  if (!geom_is_accessible) {
    LOG("gevgen_hnl", pFATAL)
      << "The specified ROOT geometry doesn't exist! Initialization failed!";
    exit(1);
  }

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 

  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );
  TGeoShape * ts  = top_volume->GetShape();

  TGeoBBox *  box = (TGeoBBox *)ts;

  //get box origin and dimensions (in the same units as the geometry)
  fdx = box->GetDX();
  fdy = box->GetDY();
  fdz = box->GetDZ();
  fox = (box->GetOrigin())[0];
  foy = (box->GetOrigin())[1];
  foz = (box->GetOrigin())[2];

  LOG("gevgen_hnl", pDEBUG)
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

  LOG("gevgen_hnl", pINFO)
    << "Initialised bounding box successfully.";

}
#endif // #ifdef __CAN_USE_ROOT_GEOM__
//............................................................................
//_________________________________________________________________________________________
//............................................................................
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
void FillFluxNonsense( FluxContainer &ggn )
{
  ggn.evtno = -9999;

  ggn.pdg = -9999;
  ggn.parPdg = -9999;
  ggn.lepPdg = -9999;
  ggn.nuPdg = -9999;

  ggn.prodChan = -9999;
  ggn.nuProdChan = -9999;

  ggn.startPoint.SetXYZ(-9999.9, -9999.9, -9999.9);
  ggn.targetPoint.SetXYZ(-9999.9, -9999.9, -9999.9);
  ggn.startPointUser.SetXYZ(-9999.9, -9999.9, -9999.9);
  ggn.targetPointUser.SetXYZ(-9999.9, -9999.9, -9999.9);
  ggn.delay = -9999.9;

  ggn.polz.SetXYZ(-9999.9, -9999.9, -9999.9);

  ggn.p4.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);
  ggn.parp4.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);
  ggn.p4User.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);
  ggn.parp4User.SetPxPyPzE(-9999.9, -9999.9, -9999.9, -9999.9);

  ggn.Ecm = -9999.9;
  ggn.nuEcm = -9999.9;
  ggn.XYWgt = -9999.9;
  ggn.boostCorr = -9999.9;
  ggn.accCorr = -9999.9;
  ggn.zetaMinus = -9999.9;
  ggn.zetaPlus = -9999.9;
  ggn.acceptance = -9999.9;

  ggn.nimpwt = -9999.9;

  return;
}
//_________________________________________________________________________________________
void FillFlux( FluxContainer &ggn, FluxContainer &tgn )
{
  ggn.evtno = tgn.evtno;

  ggn.pdg = tgn.pdg;
  ggn.parPdg = tgn.parPdg;
  ggn.lepPdg = tgn.lepPdg;
  ggn.nuPdg = tgn.nuPdg;

  ggn.prodChan = tgn.prodChan;
  ggn.nuProdChan = tgn.nuProdChan;

  ggn.startPoint = tgn.startPoint;
  ggn.targetPoint = tgn.targetPoint;
  ggn.startPointUser = tgn.startPointUser;
  ggn.targetPointUser = tgn.targetPointUser;
  ggn.delay = tgn.delay;

  ggn.polz = tgn.polz;

  ggn.p4 = tgn.p4;
  ggn.parp4 = tgn.parp4;
  ggn.p4User = tgn.p4User;
  ggn.parp4User = tgn.parp4User;

  ggn.Ecm = tgn.Ecm;
  ggn.nuEcm = tgn.nuEcm;
  ggn.XYWgt = tgn.XYWgt;
  ggn.boostCorr = tgn.boostCorr;
  ggn.accCorr = tgn.accCorr;
  ggn.zetaMinus = tgn.zetaMinus;
  ggn.zetaPlus = tgn.zetaPlus;
  ggn.acceptance = tgn.acceptance;

  ggn.nimpwt = tgn.nimpwt;
}
//............................................................................
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
//_________________________________________________________________________________________
TLorentzVector GeneratePosition( GHepRecord * event )
{ 
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
//_________________________________________________________________________________________
const EventRecordVisitorI * HNLGenerator(void)
{
  //string sname   = "genie::EventGenerator";
  //string sconfig = "BeamHNL";
  AlgFactory * algf = AlgFactory::Instance();

  LOG("gevgen_hnl", pINFO)
    << "Instantiating HNL generator.";

  const Algorithm * algmcgen = algf->GetAlgorithm(kDefOptSName, kDefOptSConfig);
  LOG("gevgen_hnl", pDEBUG)
    << "Got algorithm " << kDefOptSName.c_str() << "/" << kDefOptSConfig.c_str();;

  const EventRecordVisitorI * mcgen = 
    dynamic_cast< const EventRecordVisitorI * >( algmcgen );
  if(!mcgen) {
     LOG("gevgen_hnl", pFATAL) << "Couldn't instantiate the HNL generator";
     gAbortingInErr = true;
     exit(1);
  }

  LOG("gevgen_hnl", pINFO)
    << "HNL generator instantiated successfully.";

  return mcgen;
}
//_________________________________________________________________________________________
int SelectDecayMode( std::vector< HNLDecayMode_t > * intChannels, SimpleHNL sh )
{
  const std::map< HNLDecayMode_t, double > gammaMap = sh.GetValidChannels();

  // set CoM lifetime now if unset
  if( CoMLifetime < 0.0 ){
    CoMLifetime = sh.GetCoMLifetime();
  }

  std::vector< HNLDecayMode_t > intAndValidChannels;
  for( std::vector< HNLDecayMode_t >::iterator it = intChannels->begin(); it != intChannels->end(); ++it ){
    HNLDecayMode_t mode = *it;
    
    //check if this is a valid mode
    if( !utils::hnl::IsKinematicallyAllowed( mode, gOptMassHNL ) ) continue;

    intAndValidChannels.emplace_back( mode );
    auto mapG = gammaMap.find( mode );
    double theGamma = mapG->second;
    LOG("gevgen_hnl", pDEBUG)
      << "For mode " << utils::hnl::AsString( mode ) << " gamma = " << theGamma;
  }

  if( intAndValidChannels.size() == 0 ){ // all the modes picked by user are too heavy. Abort.
    LOG( "gevgen_hnl", pFATAL )
      << "None of the channels specified as interesting are kinematically allowed. Please either increase the HNL mass or change interesting channels in config.";
    exit(1);
  }
  
  std::map< HNLDecayMode_t, double > intMap =
    selector::SetInterestingChannels( intAndValidChannels, gammaMap );
     
  sh.SetInterestingChannels( intMap );

  // update fraction of total decay width that is not in inhibited channels
  double gammaAll = 0.0, gammaInt = 0.0;
  for( std::map< HNLDecayMode_t, double >::const_iterator itall = gammaMap.begin() ;
       itall != gammaMap.end() ; ++itall ){
    gammaAll += (*itall).second;
  }
  for( std::map< HNLDecayMode_t, double >::iterator itint = intMap.begin() ;
       itint != intMap.end() ; ++itint ){
    gammaInt += (*itint).second;
  }
  assert( gammaInt > 0.0 && gammaAll >= gammaInt );
  decayMod = gammaInt / gammaAll;

  // get probability that channels in intAndValidChannels will be selected
  std::map< HNLDecayMode_t, double > PMap = 
    selector::GetProbabilities( intMap );
     
  // now do a random throw
  RandomGen * rnd = RandomGen::Instance();
  double ranThrow = rnd->RndGen().Uniform( 0., 1. ); // HNL's fate is sealed.

  HNLDecayMode_t selectedDecayChan =
    selector::SelectChannelInclusive( PMap, ranThrow );

  int decay = ( int ) selectedDecayChan;
  return decay;
}
//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_hnl", pINFO) << "Parsing command line arguments";

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);
  
  // force the app to look at HNL tune by default
  // if user passes --tune argument, look at the user input tune instead
  char * expargv[ argc + 2 ];
  bool didFindUserInputTune = false;
  std::string stExtraTuneBit = kDefOptSTune;
  
  if( parser.OptionExists("tune") ){
    didFindUserInputTune = true;
    stExtraTuneBit = parser.ArgAsString("tune");
    LOG( "gevgen_hnl", pWARN )
      << "Using input HNL tune " << parser.ArgAsString("tune");
  } else {
    LOG( "gevgen_hnl", pWARN )
      << "Using default HNL tune " << kDefOptSTune;
  }
  // append this to argv
  for( int iArg = 0; iArg < argc; iArg++ ){
    expargv[iArg] = argv[iArg];
  }
  if( !didFindUserInputTune ){
    char * chBit = const_cast< char * >( stExtraTuneBit.c_str() ); // ugh. Ugly.
    std::string stune("--tune"); char * tBit = const_cast< char * >( stune.c_str() );
    expargv[argc] = tBit;
    expargv[argc+1] = chBit;
  }

  // Common run options.
  int expargc = ( didFindUserInputTune ) ? argc : argc+2;
  std::string stnull(""); char * nBit = const_cast< char * >( stnull.c_str() );
  expargv[expargc] = nBit;

  RunOpt::Instance()->ReadFromCommandLine(expargc,expargv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
    PrintSyntax();
    exit(0);
  }

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevgen_hnl", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_hnl", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_hnl", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_hnl", pFATAL)
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  // get HNL mass 
  const Algorithm * algHNLGen = AlgFactory::Instance()->GetAlgorithm("genie::hnl::Decayer", "Default");
  const Decayer * hnlgen = dynamic_cast< const Decayer * >( algHNLGen );
  gOptMassHNL = hnlgen->GetHNLMass();

  bool isMonoEnergeticFlux = true;
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
  if( parser.OptionExists('f') ) {
    LOG("gevgen_hnl", pDEBUG)
      << "A flux has been offered. Searching this path: " << parser.ArgAsString('f');
    isMonoEnergeticFlux = false;
    gOptFluxFilePath = parser.ArgAsString('f');
    
    // check if this is valid path (assume these are dk2nu files)
    if( gSystem->OpenDirectory( gOptFluxFilePath.c_str() ) != NULL ){
      gOptIsUsingDk2nu = true;
      LOG("gevgen_hnl", pDEBUG)
	<< "dk2nu flux files detected. Will create flux spectrum dynamically.";
    } else {
      LOG("gevgen_hnl", pFATAL)
	<< "Invalid flux file path " << gOptFluxFilePath;
      exit(1);
    }
  } else {
    // we need the 'E' option! Log it and pass below
    LOG("gevgen_hnl", pINFO)
      << "No flux file offered. Assuming monoenergetic flux.";
  } //-f
#endif // #ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

  // HNL energy (only relevant if we do not have an input flux)
  gOptEnergyHNL = -1;
  if( isMonoEnergeticFlux ){
    if( parser.OptionExists('E') ) {
      LOG("gevgen_hnl", pDEBUG)
        << "Reading HNL energy";
      gOptEnergyHNL = parser.ArgAsDouble('E');
    } else {
      LOG("gevgen_hnl", pFATAL)
        << "You need to specify the HNL energy";
      PrintSyntax();
      exit(0);
    } //-E
    assert(gOptEnergyHNL > gOptMassHNL);
  }

  gOptIsMonoEnFlux = isMonoEnergeticFlux;

  // first flux entry to read
  if( parser.OptionExists("firstEvent") ) {
    gOptFirstEvent = parser.ArgAsInt("firstEvent");
    LOG( "gevgen_hnl", pINFO )
      << "Starting flux readin at first event = " << gOptFirstEvent;
  } // --firstEvent

  // HNL decay mode
  int mode = -1;
  if( parser.OptionExists('m') ) {
    LOG("gevgen_hnl", pDEBUG)
        << "Reading HNL decay mode";
    mode = parser.ArgAsInt('m');
  } else {
    LOG("gevgen_hnl", pINFO)
        << "No decay mode specified - will sample from allowed decay modes";
  } //-m
  gOptDecayMode = (HNLDecayMode_t) mode;

  bool allowed = utils::hnl::IsKinematicallyAllowed(gOptDecayMode, gOptMassHNL);
  if(!allowed) {
    LOG("gevgen_hnl", pFATAL)
      << "Specified decay is not allowed kinematically for the given HNL mass";
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
    LOG("gevgen_hnl", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');

    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file =
            ! (gSystem->AccessPathName(geom.c_str()));
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    } else {
      LOG("gevgen_hnl", pFATAL)
	<< "Geometry option is not a ROOT file. Please use ROOT geom.";
      PrintSyntax();
      exit(1);
    }
  } else {
      // LOG("gevgen_hnl", pFATAL)
      //   << "No geometry option specified - Exiting";
      // PrintSyntax();
      // exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // length units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_hnl", pDEBUG)
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_hnl", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // // density units:
     // if( parser.OptionExists('D') ) {
     //    LOG("gevgen_hnl", pDEBUG)
     //       << "Checking for input geometry density units";
     //    dunits = parser.ArgAsString('D');
     // } else {
     //    LOG("gevgen_hnl", pDEBUG) << "Using default geometry density units";
     //    dunits = kDefOptGeomDUnits;
     // } // -D
     gOptGeomLUnits = utils::units::UnitFromString(lunits);
     // gOptGeomDUnits = utils::units::UnitFromString(dunits);

  } // using root geom?
#endif // #ifdef __CAN_USE_ROOT_GEOM__

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_hnl", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_hnl", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_hnl", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_hnl", pINFO) << "Unspecified random number seed - Using default";
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

  LOG("gevgen_hnl", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_hnl job configuration");

  LOG("gevgen_hnl", pNOTICE)
     << "\n @@ Run number    : " << gOptRunNu
     << "\n @@ Random seed   : " << gOptRanSeed
     << "\n @@ HNL mass      : " << gOptMassHNL << " GeV"
     << "\n @@ Decay channel : " << utils::hnl::AsString(gOptDecayMode)
     << "\n @@ Geometry      : " << gminfo.str()
     << "\n @@ Statistics    : " << gOptNev << " events";
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_hnl", pFATAL)
   << "\n **Syntax**"
   << "\n gevgen_hnl [-h] "
   << "\n            [-r run#]"
   << "\n             -n n_of_events"
   << "\n             -f path/to/flux/files"
   << "\n            [-E hnl_energy]"
   << "\n            [--firstEvent first_event_for_dk2nu_readin]"  
   << "\n            [-m decay_mode]"
   << "\n            [-g geometry (ROOT file)]"
   << "\n            [-L length_units_at_geom]"
   << "\n            [-o output_event_file_prefix]"
   << "\n            [--seed random_number_seed]"
   << "\n            [--message-thresholds xml_file]"
   << "\n            [--event-record-print-level level]"
   << "\n            [--mc-job-status-refresh-rate  rate]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << " or look at the source code: $GENIE/src/Apps/gBeamHNLEvGen.cxx"
   << "\n";
}
//_________________________________________________________________________________________
