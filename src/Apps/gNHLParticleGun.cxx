//________________________________________________________________________________________
/*!

\program gevgen_pgnhl

\brief   A particle gun for NHL.

         *** Synopsis :

         gevgen_pgnhl [-h]
                     [-r run#]
                      -n n_of_events
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
	 John Plows <komninos-john.plows \at physics.ox.ac.uk>
	 University of Oxford

\created May 27th, 2022

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
//#include "Physics/NeutralHeavyLepton/NHLFluxReader.h"
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
using namespace genie::NHL::NHLenums;

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

TLorentzVector * GenerateOriginPosition( GHepRecord * event );
TLorentzVector * GenerateOriginMomentum( GHepRecord * event );

TLorentzVector GeneratePosition( GHepRecord * event );
#ifdef __CAN_USE_ROOT_GEOM__
void  InitBoundingBox    (void);
#endif // #ifdef __CAN_USE_ROOT_GEOM__

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format
string          kDefOptEvFilePrefix = "gntp";

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

//NHLPrimaryVtxGenerator * nhlgen = 0;
// HNL lifetime in rest frame
double CoMLifetime = -1.0; // GeV^{-1}
// an array to keep production vertex
double evProdVtx[3] = {0.0, 0.0, 0.0}; // mm
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
  const Algorithm * algNHLGen = AlgFactory::Instance()->GetAlgorithm("genie::NHLPrimaryVtxGenerator", "Default");
  const Algorithm * algDkVol = AlgFactory::Instance()->GetAlgorithm("genie::NHL::NHLDecayVolume", "Default");
  
  const NHLPrimaryVtxGenerator * nhlgen = dynamic_cast< const NHLPrimaryVtxGenerator * >( algNHLGen );
  const NHLDecayVolume * dkVol = dynamic_cast< const NHLDecayVolume * >( algDkVol );

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 

  // RETHERE implement top volume option from cmd line
  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );
  TGeoShape * ts  = top_volume->GetShape();

  TGeoBBox *  box = (TGeoBBox *)ts;

  // pass this box to NHLDecayVolume
  dkVol->ImportBoundingBox( box );

  string confString = kDefOptSName + kDefOptSConfig;
  //const double confMass = nhlgen->GetNHLMass( confString );
  //const std::vector< double > confCoups = nhlgen->GetNHLCouplings( confString );

  SimpleNHL confsh = nhlgen->GetNHLInstance( confString );
  const double confMass = confsh.GetMass();
  const std::vector< double > confCoups = confsh.GetCouplings();
  const bool confIsMajorana = confsh.GetIsMajorana();
  const int confType = confsh.GetType();
  const std::vector< NHLDecayMode_t > confIntChan = confsh.GetInterestingChannelsVec();

  LOG( "gevgen_pgnhl", pDEBUG )
    << "At app stage we see:"
    << "\nMass = " << confMass << " GeV"
    << "\nECoup = " << confCoups.at(0)
    << "\nMCoup = " << confCoups.at(1)
    << "\nTCoup = " << confCoups.at(2)
    << "\nIsMajorana = " << confIsMajorana
    << "\nType = " << confType;

  gOptECoupling = confCoups.at(0);
  gOptMCoupling = confCoups.at(1);
  gOptTCoupling = confCoups.at(2);
  gOptNHLKind = confType; // for mixing
  gOptIsMajorana = confIsMajorana;

  gOptIntChannels = confIntChan;

  assert( gOptECoupling >= 0.0 && gOptMCoupling >= 0.0 && gOptTCoupling >= 0.0 );

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu, gOptRanSeed);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  LOG("gevgen_pgnhl", pNOTICE)
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

  LOG("gevgen_pgnhl", pNOTICE)
    << "Initialised MC job monitor";

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

#ifdef __CAN_USE_ROOT_GEOM__
  // Read geometry bounding box - for vertex position generation
  if( gOptUsingRootGeom ){
    InitBoundingBox();
  }
#endif // #ifdef __CAN_USE_ROOT_GEOM__

  // read in energy. This is always constant
  gOptEnergyNHL = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-Energy" );
  if( gOptEnergyNHL <= gOptMassNHL ){
    LOG( "gevgen_pgnhl", pFATAL )
      << "\nInsufficient energy " << gOptEnergyNHL << " GeV for NHL of mass " << gOptMassNHL << " GeV. Exiting."
      << "\nPlease check ${GENIE}/config/CommonNHL.xml sections \"ParamSpace\" and \"ParticleGun\"";
    exit(1);
  }
  assert( gOptEnergyNHL > gOptMassNHL );

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

    if((ievent) == gOptNev) break;
      
     LOG("gevgen_pgnhl", pNOTICE)
          << " *** Generating event............ " << ievent;

     int hpdg = genie::kPdgNHL;
     EventRecord * event = new EventRecord;

     int decay  = (int) gOptDecayMode;

     LOG("gevgen_pgnhl", pDEBUG)
       << " Building SimpleNHL object ";
     SimpleNHL sh( "NHL", ievent, hpdg, genie::kPdgKP, 
		   gOptMassNHL, gOptECoupling, gOptMCoupling, gOptTCoupling, false );
     
     const std::map< NHLDecayMode_t, double > gammaMap = sh.GetValidChannels();
     LOG( "gevgen_pgnhl", pDEBUG )
       << "CoMLifetime = " << sh.GetCoMLifetime();
     CoMLifetime = sh.GetCoMLifetime();
     
     if( gOptDecayMode == kNHLDcyNull ){ // select from available modes
       LOG("gevgen_pgnhl", pNOTICE)
	 << "No decay mode specified - sampling from all available modes.";

       LOG("gevgen_pgnhl", pDEBUG)
	 << "Reading interesting channels vector from config";
       std::vector< NHLDecayMode_t > * intChannels = &gOptIntChannels;
       
       decay = SelectDecayMode( intChannels, sh );
     }

     Interaction * interaction = Interaction::NHL(genie::kPdgNHL, gOptEnergyNHL, decay);
     event->AttachSummary(interaction);

     // generate origin position and momentum direction for this NHL
     TLorentzVector * x4orig = GenerateOriginPosition( event );
     TLorentzVector * p4NHL  = GenerateOriginMomentum( event );

     LOG("gevgen_pgnhl", pDEBUG)
       << "Origin information for NHL at event number " << ievent << ":"
       << "\n Origin x4 = " << utils::print::X4AsString( x4orig )
       << "\n Momentum  = " << utils::print::P4AsString( p4NHL );
     
     if( event->Vertex() ){
       LOG( "gevgen_pgnhl", pDEBUG )
	 << "\nIS p4  = " << utils::print::P4AsString( interaction->InitStatePtr()->GetProbeP4() )
	 << "\nIS vtx = " << utils::print::X4AsString( event->Vertex() );
     } else {
       LOG( "gevgen_pgnhl", pDEBUG )
	 << "\nIS p4  = " << utils::print::P4AsString( interaction->InitStatePtr()->GetProbeP4() );
     }

     LOG("gevgen_pgnhl", pDEBUG)
       << "Note decay mode is " << utils::nhl::AsString(gOptDecayMode);

     // Simulate decay
     nhlgen->ProcessEventRecord(event);
     dkVol->SetStartingParameters( event, CoMLifetime, false, gOptUsingRootGeom, gOptRootGeom.c_str() );
     dkVol->ProcessEventRecord(event);

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
     event->SetWeight( evWeight );

     // why does InitState show the wrong p4 here?
     interaction->InitStatePtr()->SetProbeP4( *(event->Particle(0)->P4()) );
     LOG( "gevgen_pgnhl", pDEBUG ) << 
       "\n!*!*!*! " << utils::print::P4AsString( event->Particle(0)->P4() );
     
     LOG("gevgen_pgnhl", pDEBUG) << "Weight = " << evWeight;

     LOG("gevgen_pgnhl", pINFO)
         << "Generated event: " << *event;

     // Add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     delete event;

     ievent++;
  } // event loop

  // Save the generated event tree & close the output file
  ntpw.Save();

  LOG("gevgen_pgnhl", pNOTICE) << "Done!";

  return 0;
}
//_________________________________________________________________________________________
//............................................................................
#ifdef __CAN_USE_ROOT_GEOM__
void InitBoundingBox(void)
{
// Initialise geometry bounding box, used for generating NHL vertex positions

  LOG("gevgen_pgnhl", pINFO)
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
    LOG("gevgen_pgnhl", pFATAL)
      << "The specified ROOT geometry doesn't exist! Initialization failed!";
    exit(1);
  }

  if( !gOptRootGeoManager ) gOptRootGeoManager = TGeoManager::Import(gOptRootGeom.c_str()); 

  // RETHERE implement top volume option from cmd line
  TGeoVolume * top_volume = gOptRootGeoManager->GetTopVolume();
  assert( top_volume );
  TGeoShape * ts  = top_volume->GetShape();

  TGeoBBox *  box = (TGeoBBox *)ts;
  
  const Algorithm * algDkVol = AlgFactory::Instance()->GetAlgorithm("genie::NHL::NHLDecayVolume", "Default");

  const NHLDecayVolume * dkVol = dynamic_cast< const NHLDecayVolume * >( algDkVol );
  
  // pass this box to NHLDecayVolume
  dkVol->ImportBoundingBox( box );

  //get box origin and dimensions (in the same units as the geometry)
  fdx = box->GetDX();
  fdy = box->GetDY();
  fdz = box->GetDZ();
  fox = (box->GetOrigin())[0];
  foy = (box->GetOrigin())[1];
  foz = (box->GetOrigin())[2];

  LOG("gevgen_pgnhl", pINFO)
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

  LOG("gevgen_pgnhl", pINFO)
    << "Initialised bounding box successfully.";

}
#endif // #ifdef __CAN_USE_ROOT_GEOM__
//............................................................................
//_________________________________________________________________________________________
TLorentzVector * GenerateOriginMomentum( GHepRecord * event )
{
  double E = gOptEnergyNHL;
  double M = gOptMassNHL;
  double pMag = std::sqrt( E*E - M*M );

  // first map directional cosines to unit sphere
  double cx = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-cx" );
  double cy = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-cy" );
  double cz = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-cz" );
  double c2 = std::sqrt( cx*cx + cy*cy + cz*cz );
  assert( c2 > 0.0 );

  cx *= 1.0 / c2; cy *= 1.0 / c2; cz *= 1.0 / c2;

  double theta = TMath::ACos( cz / c2 );
  double phi   = ( TMath::Sin( theta ) != 0.0 ) ? TMath::ACos( cx / ( c2 * TMath::Sin( theta ) ) ) : 0.0;

  // apply uniform random deviation
  double dthetaDeg = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-DTheta" );
  double dphiDeg   = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-DPhi" );
  double dThetaMax = dthetaDeg * constants::kPi / 180.0;
  double dPhiMax   = dphiDeg * constants::kPi / 180.0;

  RandomGen * rng = RandomGen::Instance();
  double dTheta = rng->RndGen().Uniform( -dThetaMax, dThetaMax );
  double dPhi   = rng->RndGen().Uniform( -dPhiMax, dPhiMax );
  
  std::ostringstream asts;
  asts
    << "Output details for the momentum:"
    << "\nDirectional cosines: ( " << cx << ", " << cy << ", " << cz << " )"
    << "\nMapped to ( " << theta * 180.0 / constants::kPi 
    << ", " << phi * 180.0 / constants::kPi << " ) [deg] on unit sphere"
    << "\nApplied deviation: dTheta = " << dTheta * 180.0 / constants::kPi
    << ", dPhi = " << dPhi * 180.0 / constants::kPi << " [deg]";

  // make momentum
  theta += dTheta; phi += dPhi;
  TLorentzVector * p4NHL = new TLorentzVector();
  p4NHL->SetPxPyPzE( pMag * TMath::Sin( theta ) * TMath::Cos( phi ),
		     pMag * TMath::Sin( theta ) * TMath::Sin( phi ),
		     pMag * TMath::Cos( theta ), E );
  
  asts << "\nMomentum = " << utils::print::P4AsString( p4NHL );

  Interaction * interaction = event->Summary();
  interaction->InitStatePtr()->SetProbeP4( *p4NHL );

  LOG( "gevgen_pgnhl", pDEBUG ) << asts.str();

  //delete rng;
  return p4NHL;
}
//_________________________________________________________________________________________
TLorentzVector * GenerateOriginPosition( GHepRecord * event )
{
  // get centre of box from config
  double ox = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-OriginX" );
  double oy = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-OriginY" );
  double oz = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-OriginZ" );

  // allow uniform deviation
  double Dxmax = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-OriginDX" );
  double Dymax = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-OriginDY" );
  double Dzmax = utils::nhl::GetCfgDouble( "NHL", "ParticleGun", "PG-OriginDZ" );

  RandomGen * rng = RandomGen::Instance();
  double dx = rng->RndGen().Uniform( -Dxmax, Dxmax );
  double dy = rng->RndGen().Uniform( -Dymax, Dymax );
  double dz = rng->RndGen().Uniform( -Dzmax, Dzmax );

  std::ostringstream asts;
  asts
    << "Output details for the origin position:"
    << "\nCentred around user ( " << ox << ", " << oy << ", " << oz << " ) [m]"
    << "\nApplied deviation: ( " << dx << ", " << dy << ", " << dz << " ) [m]";

  // make position
  ox += dx; oy += dy; oz += dz;
  TLorentzVector * x4NHL = new TLorentzVector();
  x4NHL->SetXYZT( ox, oy, oz, 0.0 );

  asts << "\nx4NHL = " << utils::print::X4AsString( x4NHL );

  LOG( "gevgen_pgnhl", pDEBUG ) << asts.str();

  event->SetVertex( *x4NHL );
  //delete rng;
  return x4NHL;
}
//_________________________________________________________________________________________
TLorentzVector GeneratePosition( GHepRecord * event )
{
  if( gOptUsingRootGeom ){
  
    const Algorithm * algDkVol = AlgFactory::Instance()->GetAlgorithm("genie::NHL::NHLDecayVolume", "Default");
    
    const NHLDecayVolume * dkVol = dynamic_cast< const NHLDecayVolume * >( algDkVol );
    dkVol->SetStartingParameters( event, CoMLifetime, true, gOptUsingRootGeom, gOptRootGeom.c_str() );
    
    dkVol->ProcessEventRecord( event );
    
    TLorentzVector x4 = *(event->Vertex());
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

  LOG("gevgen_pgnhl", pINFO)
    << "Instantiating NHL generator.";

  const Algorithm * algmcgen = algf->GetAlgorithm(kDefOptSName, kDefOptSConfig);
  LOG("gevgen_pgnhl", pDEBUG)
    << "Got algorithm " << kDefOptSName.c_str() << "/" << kDefOptSConfig.c_str();;

  const EventRecordVisitorI * mcgen = 
    dynamic_cast< const EventRecordVisitorI * >( algmcgen );
  if(!mcgen) {
     LOG("gevgen_pgnhl", pFATAL) << "Couldn't instantiate the NHL generator";
     gAbortingInErr = true;
     exit(1);
  }

  LOG("gevgen_pgnhl", pINFO)
    << "NHL generator instantiated successfully.";

  return mcgen;
}
//_________________________________________________________________________________________
int SelectDecayMode( std::vector< NHLDecayMode_t > * intChannels, SimpleNHL sh )
{
  LOG("gevgen_pgnhl", pDEBUG)
    << " Getting valid channels ";
  const std::map< NHLDecayMode_t, double > gammaMap = sh.GetValidChannels();

  // set CoM lifetime now if unset
  if( CoMLifetime < 0.0 ){
    CoMLifetime = sh.GetCoMLifetime();
    LOG( "gevgen_pgnhl", pDEBUG )
      << "Rest frame CoMLifetime = " << CoMLifetime << " [GeV^{-1}]";
  }

  for( std::vector< NHLDecayMode_t >::iterator it = intChannels->begin(); it != intChannels->end(); ++it ){
    NHLDecayMode_t mode = *it;
    auto mapG = gammaMap.find( mode );
    double theGamma = mapG->second;
    LOG("gevgen_pgnhl", pDEBUG)
      << "For mode " << utils::nhl::AsString( mode ) << " gamma = " << theGamma;
  }

  LOG("gevgen_pgnhl", pDEBUG)
    << " Setting interesting channels map ";
  std::map< NHLDecayMode_t, double > intMap =
    NHLSelector::SetInterestingChannels( (*intChannels), gammaMap );
     
  LOG("gevgen_pgnhl", pDEBUG)
    << " Telling SimpleNHL about interesting channels ";
  sh.SetInterestingChannels( intMap );

  // get probability that channels in intChannels will be selected
  LOG("gevgen_pgnhl", pDEBUG)
    << " Building probablilities of interesting channels ";
  std::map< NHLDecayMode_t, double > PMap = 
    NHLSelector::GetProbabilities( intMap );
     
  // now do a random throw
  RandomGen * rnd = RandomGen::Instance();
  double ranThrow = rnd->RndGen().Uniform( 0., 1. ); // NHL's fate is sealed.

  LOG("gevgen_pgnhl", pDEBUG)
    << "Random throw = " << ranThrow;

  NHLDecayMode_t selectedDecayChan =
    NHLSelector::SelectChannelInclusive( PMap, ranThrow );

  int decay = ( int ) selectedDecayChan;
  return decay;
}
//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_pgnhl", pINFO) << "Parsing command line arguments";

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
    LOG("gevgen_pgnhl", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_pgnhl", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_pgnhl", pDEBUG)
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_pgnhl", pFATAL)
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  // get NHL mass directly from config
  gOptMassNHL = genie::utils::nhl::GetCfgDouble( "NHL", "ParameterSpace", "NHL-Mass" );

  // NHL decay mode
  int mode = -1;
  if( parser.OptionExists('m') ) {
    LOG("gevgen_pgnhl", pDEBUG)
        << "Reading NHL decay mode";
    mode = parser.ArgAsInt('m');
  } else {
    LOG("gevgen_pgnhl", pINFO)
        << "No decay mode specified - will sample from allowed decay modes";
  } //-m
  gOptDecayMode = (NHLDecayMode_t) mode;

  bool allowed = utils::nhl::IsKinematicallyAllowed(gOptDecayMode, gOptMassNHL);
  if(!allowed) {
    LOG("gevgen_pgnhl", pFATAL)
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
    LOG("gevgen_pgnhl", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');

    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file =
            ! (gSystem->AccessPathName(geom.c_str()));
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    } else {
      LOG("gevgen_pgnhl", pFATAL)
	<< "Geometry option is not a ROOT file. This is a work in progress; please use ROOT geom.";
      PrintSyntax();
      exit(1);
    }
  } else {
      // LOG("gevgen_pgnhl", pFATAL)
      //   << "No geometry option specified - Exiting";
      // PrintSyntax();
      // exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // length units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_pgnhl", pDEBUG)
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_pgnhl", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // // density units:
     // if( parser.OptionExists('D') ) {
     //    LOG("gevgen_pgnhl", pDEBUG)
     //       << "Checking for input geometry density units";
     //    dunits = parser.ArgAsString('D');
     // } else {
     //    LOG("gevgen_pgnhl", pDEBUG) << "Using default geometry density units";
     //    dunits = kDefOptGeomDUnits;
     // } // -D
     gOptGeomLUnits = utils::units::UnitFromString(lunits);
     // gOptGeomDUnits = utils::units::UnitFromString(dunits);

     // check whether an event generation volume name has been
     // specified -- default is the 'top volume'
     if( parser.OptionExists('t') ) {
        LOG("gevgen_pgnhl", pDEBUG) << "Checking for input volume name";
        gOptRootGeomTopVol = parser.ArgAsString('t');
     } else {
        LOG("gevgen_pgnhl", pDEBUG) << "Using the <master volume>";
     } // -t

  } // using root geom?
#endif // #ifdef __CAN_USE_ROOT_GEOM__

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_pgnhl", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_pgnhl", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_pgnhl", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_pgnhl", pINFO) << "Unspecified random number seed - Using default";
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

  LOG("gevgen_pgnhl", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_pgnhl job configuration");

  LOG("gevgen_pgnhl", pNOTICE)
     << "\n @@ Run number    : " << gOptRunNu
     << "\n @@ Random seed   : " << gOptRanSeed
     << "\n @@ Decay channel : " << utils::nhl::AsString(gOptDecayMode)
     << "\n @@ Geometry      : " << gminfo.str()
     << "\n @@ Statistics    : " << gOptNev << " events";
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_pgnhl", pFATAL)
   << "\n **Syntax**"
   << "\n gevgen_pgnhl [-h] "
   << "\n              [-r run#]"
   << "\n               -n n_of_events"
   << "\n              [-m decay_mode]"
   << "\n              [-g geometry (ROOT file)]"
   << "\n              [-t top_volume_name_at_geom]"
   << "\n              [-L length_units_at_geom]"
   << "\n              [-o output_event_file_prefix]"
   << "\n              [--seed random_number_seed]"
   << "\n              [--message-thresholds xml_file]"
   << "\n              [--event-record-print-level level]"
   << "\n              [--mc-job-status-refresh-rate  rate]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << " or look at the source code: $GENIE/src/Apps/gNHLParticleGun.cxx"
   << "\n";
}
//_________________________________________________________________________________________
