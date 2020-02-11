//________________________________________________________________________________________
/*!

\program gevgen_nhl

\brief   A GENIE-based neutral heavy lepton event generation application.

         *** Synopsis :

         gevgen_nhl [-h]
                   [-r run#]
                    -n n_of_events
                    -E nhl_energy (temporary)
                    --mass nhl_mass
                    -m decay_mode
                    -f nhl_flux
	                 [-g geometry]
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
              *** not implemented for gevgen_nhl yet ***
           -g
              Input detector geometry.
              *** not implemented for gevgen_nhl yet ***
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

\cpright Copyright (c) 2003-2020, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org

*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>

#include <TSystem.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <TGeoShape.h>
#include <TGeoBBox.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"
#include "Physics/NeutralHeavyLepton/NHLDecayUtils.h"
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

// function prototypes
void  GetCommandLineArgs (int argc, char ** argv);
void  PrintSyntax        (void);
void  InitBoundingBox    (void);
TLorentzVector GeneratePosition(void);
const EventRecordVisitorI * NHLGenerator(void);

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format
string          kDefOptEvFilePrefix = "gntp";

//
Long_t           gOptRunNu        = 1000;                // run number
int              gOptNev          = 10;                  // number of events to generate
double           gOptEnergyNHL    = -1;                  // NHL mass
double           gOptMassNHL      = -1;                  // NHL mass
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

//_________________________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  // Init messenger and random number seed
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // Create a MC job monitor for a periodically updated status file
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  // Read geometry bounding box - for vertex position generation
  InitBoundingBox();

  // Get the nucleon decay generator
  const EventRecordVisitorI * mcgen = NHLGenerator();

  // Event loop
  int ievent = 0;
  while (1)
  {
     if(ievent == gOptNev) break;

     LOG("gevgen_nhl", pNOTICE)
          << " *** Generating event............ " << ievent;

     EventRecord * event = new EventRecord;
     // int target = SelectInitState();
     int decay  = (int)gOptDecayMode;
     Interaction * interaction = Interaction::NHL(gOptEnergyNHL, decay);
     event->AttachSummary(interaction);

     // Simulate decay
     mcgen->ProcessEventRecord(event);

     // Generate a position within the geometry bounding box
     TLorentzVector x4 = GeneratePosition();
     event->SetVertex(x4);

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
void InitBoundingBox(void)
{
// Initialise geometry bounding box, used for generating NHL vertex positions

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
  TGeoShape * ts  = top_volume->GetShape();
  TGeoBBox *  box = (TGeoBBox *)ts;
  //get box origin and dimensions (in the same units as the geometry)
  fdx = box->GetDX(); // half-length
  fdy = box->GetDY(); // half-length
  fdz = box->GetDZ(); // half-length
  fox = (box->GetOrigin())[0];
  foy = (box->GetOrigin())[1];
  foz = (box->GetOrigin())[2];

  // Convert from local to SI units
  fdx *= gOptGeomLUnits;
  fdy *= gOptGeomLUnits;
  fdz *= gOptGeomLUnits;
  fox *= gOptGeomLUnits;
  foy *= gOptGeomLUnits;
  foz *= gOptGeomLUnits;
}
//_________________________________________________________________________________________
TLorentzVector GeneratePosition(void)
{
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
//_________________________________________________________________________________________
const EventRecordVisitorI * NHLGenerator(void)
{
  string sname   = "genie::EventGenerator";
  string sconfig = "NeutralHeavyLepton";
  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * mcgen =
     dynamic_cast<const EventRecordVisitorI *> (algf->GetAlgorithm(sname,sconfig));
  if(!mcgen) {
     LOG("gevgen_nhl", pFATAL) << "Couldn't instantiate the NHL generator";
     gAbortingInErr = true;
     exit(1);
  }
  return mcgen;
}
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

  // NHL energy (temporary - will disappear once we add an option to read a flux)
  gOptEnergyNHL = -1;
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

  // NHL decay mode
  int mode = -1;
  if( parser.OptionExists('m') ) {
    LOG("gevgen_nhl", pDEBUG)
        << "Reading NHL decay mode";
    mode = parser.ArgAsInt('m');
  } else {
    LOG("gevgen_nhl", pFATAL)
        << "You need to specify the decay mode";
    PrintSyntax();
    exit(0);
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
    }
  } else {
      // LOG("gevgen_nhl", pFATAL)
      //   << "No geometry option specified - Exiting";
      // PrintSyntax();
      // exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // legth units:
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
   << "\n             -E nhl_energy (temporary)"
   << "\n             --mass nhl_mass"
   << "\n             -m decay_mode"
   << "\n             -f flux ** not installed yet **"
   << "\n            [-g geometry]"
   << "\n            [-t top_volume_name_at_geom]"
   << "\n            [-L length_units_at_geom]"
   << "\n             -n n_of_events "
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
