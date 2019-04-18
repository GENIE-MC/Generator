//________________________________________________________________________________________
/*!

\program gevgen_nosc

\brief   A GENIE-based n-nbar oscillation event generation application.

         *** Synopsis :

         gevgen_nosc [-h] 
                     [-r run#] 
                      -n n_of_events
                     [-m decay_mode]
	              -g geometry
                     [-L geometry_length_units] 
                     [-D geometry_density_units]
                     [-t geometry_top_volume_name]
                     [-o output_event_file_prefix]
                     [--seed random_number_seed]
                     [--message-thresholds xml_file]
                     [--event-record-print-level level]
                     [--mc-job-status-refresh-rate  rate]

         *** Options :

           [] Denotes an optional argument

           -h 
              Prints out the gevgen_nosc syntax and exits.
           -r 
              Specifies the MC run number [default: 1000].
           -n  
              Specifies how many events to generate.
           -m 
              Nucleon decay mode ID:
             ---------------------------------------------------------
              ID |   Decay Mode                     
                 |                                  
             ---------------------------------------------------------
               0 |    Random decay mode
               1 |    p + nbar --> \pi^{+} + \pi^{0}
               2 |    p + nbar --> \pi^{+} + 2\pi^{0}
               3 |    p + nbar --> \pi^{+} + 3\pi^{0}
               4 |    p + nbar --> 2\pi^{+} + \pi^{-} + \pi^{0}
               5 |    p + nbar --> 2\pi^{+} + \pi^{-} + 2\pi^{0}
               6 |    p + nbar --> 2\pi^{+} + \pi^{-} + 2\omega^{0}
               7 |    p + nbar --> 3\pi^{+} + 2\pi^{-} + \pi^{0}
               8 |    n + nbar --> \pi^{+} + \pi^{-}
               9 |    n + nbar --> 2\pi^{0}
              10 |    n + nbar --> \pi^{+} + \pi^{-} + \pi^{0}
              11 |    n + nbar --> \pi^{+} + \pi^{-} + 2\pi^{0}
              12 |    n + nbar --> \pi^{+} + \pi^{-} + 3\pi^{0}
              13 |    n + nbar --> 2\pi^{+} + 2\pi^{-}
              14 |    n + nbar --> 2\pi^{+} + 2\pi^{-} + \pi^{0}
              15 |    n + nbar --> \pi^{+} + \pi^{-} + \omega^{0}
              16 |    n + nbar --> 2\pi^{+} + 2\pi^{-} + 2\pi^{0}
             ---------------------------------------------------------

           -g 
              Input 'geometry'.
              This option can be used to specify any of:
              1 > A ROOT file containing a ROOT/GEANT geometry description
                  [Examples] 
                  - To use the master volume from the ROOT geometry stored 
                    in the laguna-lbno.root file, type:
                    '-g /some/path/laguna-lbno.root'
              2 > A mix of target materials, each with its corresponding weight,
                  typed as a comma-separated list of nuclear PDG codes (in the
                  std PDG2006 convention: 10LZZZAAAI) with the weight fractions
                  in brackets, eg code1[fraction1],code2[fraction2],...
                  If that option is used (no detailed input geometry description) 
                  then the interaction vertices are distributed in the detector
                  by the detector MC.
                  [Examples] 
                  - To use a target mix of 88.9% O16 and 11.1% Hydrogen type:
                    '-g 1000080160[0.889],1000010010[0.111]'
           -L 
              Input geometry length units, eg 'm', 'cm', 'mm', ...
              [default: 'mm']
           -D 
              Input geometry density units, eg 'g_cm3', 'clhep_def_density_unit',... 
              [default: 'g_cm3']
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

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
              
\created November 03, 2011
             
\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE

*/
//_________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <string> 
#include <vector>
#include <sstream>

#include <TSystem.h> 

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Physics/NNBarOscillation/NNBarOscMode.h"
#include "Physics/NNBarOscillation/NNBarOscUtils.h"
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
int   SelectAnnihilationMode (int pdg_code);
int   SelectInitState    (void);
const EventRecordVisitorI * NeutronOscGenerator(void);

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format   
string          kDefOptEvFilePrefix = "gntp";

//
Long_t             gOptRunNu        = 1000;                // run number
int                gOptNev          = 10;                  // number of events to generate
NNBarOscMode_t     gOptDecayMode    = kNONull;             // neutron oscillation mode
string             gOptEvFilePrefix = kDefOptEvFilePrefix; // event file prefix
bool               gOptUsingRootGeom = false;              // using root geom or target mix?
map<int,double>    gOptTgtMix;                             // target mix  (tgt pdg -> wght frac) / if not using detailed root geom
string             gOptRootGeom;                           // input ROOT file with realistic detector geometry
string             gOptRootGeomTopVol = "";                // input geometry top event generation volume 
double             gOptGeomLUnits = 0;                     // input geometry length units 
double             gOptGeomDUnits = 0;                     // input geometry density units 
long int           gOptRanSeed = -1;                       // random number seed

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

  // Get the nucleon decay generator
  const EventRecordVisitorI * mcgen = NeutronOscGenerator();

  // Event loop
  int ievent = 0;
  while (1)
  {
     if(ievent == gOptNev) break;

     LOG("gevgen_nnbar_osc", pNOTICE)
          << " *** Generating event............ " << ievent;

     EventRecord * event = new EventRecord;
     int target = SelectInitState();
     int decay = SelectAnnihilationMode(target);
     Interaction * interaction = Interaction::NOsc(target,decay);
     event->AttachSummary(interaction);

     // Simulate decay     
     mcgen->ProcessEventRecord(event);

     LOG("gevgen_nnbar_osc", pINFO)
         << "Generated event: " << *event;

     // Add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     delete event;

     ievent++;
  } // event loop

  // Save the generated event tree & close the output file
  ntpw.Save();

  LOG("gevgen_nnbar_osc", pNOTICE) << "Done!";

  return 0;
}
//_________________________________________________________________________________________
int SelectAnnihilationMode(int pdg_code)
{
  // if the mode is set to 'random' (the default), pick one at random!
  if (gOptDecayMode == kNORandom) {
    int mode;

    std::string pdg_string = std::to_string(static_cast<long long>(pdg_code));
    if (pdg_string.size() != 10) {
      LOG("gevgen_nnbar_osc", pERROR)
        << "Expecting PDG code to be a 10-digit integer; instead, it's the following: " << pdg_string;
      gAbortingInErr = true;
      exit(1);
    }

    // count number of protons & neutrons
    int n_nucleons = std::stoi(pdg_string.substr(6,3)) - 1;
    int n_protons  = std::stoi(pdg_string.substr(3,3));

    // factor proton / neutron ratio into branching ratios
    double proton_frac  = ((double)n_protons) / ((double)n_nucleons);
    double neutron_frac = 1 - proton_frac;

    // set branching ratios, taken from bubble chamber data
    const int n_modes = 16;
    double br [n_modes] = { 0.010, 0.080, 0.100, 0.220,
                            0.360, 0.160, 0.070, 0.020,
                            0.015, 0.065, 0.110, 0.280,
                            0.070, 0.240, 0.100, 0.100 };

    for (int i = 0; i < n_modes; i++) {
      if (i < 7)
        br[i] *= proton_frac;
      else
        br[i] *= neutron_frac;
    }

    // randomly generate a number between 1 and 0
    RandomGen * rnd = RandomGen::Instance();
    rnd->SetSeed(0);
    double p = rnd->RndNum().Rndm();

    // loop through all modes, figure out which one our random number corresponds to
    double threshold = 0;
    for (int i = 0; i < n_modes; i++) {
      threshold += br[i];
      if (p < threshold) {
        // once we've found our mode, return it!
        mode = i + 1;
        return mode;
      }
    }

    // error message, in case the random number selection fails
    LOG("gevgen_nnbar_osc", pFATAL) << "Random selection of final state failed!";
    gAbortingInErr = true;
    exit(1);
  }

  // if specific annihilation mode specified, just use that
  else {
    int mode = (int) gOptDecayMode;
    return mode;
  }
}
//_________________________________________________________________________________________
int SelectInitState(void)
{
  if (gOptTgtMix.size() > 1) {
    LOG("gevgen_nnbar_osc", pERROR)
      << "Target mix not currently supported. You must specify a single target nucleus!";
    gAbortingInErr = true;
    exit(1);
  }

  int pdg_code = gOptTgtMix.begin()->first;

  return pdg_code;
}
//_________________________________________________________________________________________
const EventRecordVisitorI * NeutronOscGenerator(void)
{
  string sname   = "genie::EventGenerator";
  string sconfig = "NNBarOsc";
  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * mcgen =
     dynamic_cast<const EventRecordVisitorI *> (algf->GetAlgorithm(sname,sconfig));
  if(!mcgen) {
     LOG("gevgen_nnbar_osc", pFATAL)
       << "Couldn't instantiate the neutron oscillation generator";
     gAbortingInErr = true;
     exit(1);
  }
  return mcgen;
}
//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_nnbar_osc", pINFO) << "Parsing command line arguments";

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
    LOG("gevgen_nnbar_osc", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_nnbar_osc", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r


  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_nnbar_osc", pDEBUG) 
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_nnbar_osc", pFATAL) 
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  // decay mode
  int mode = 0;
  if( parser.OptionExists('m') ) {
    LOG("gevgen_nnbar_osc", pDEBUG) 
        << "Reading annihilation mode";
    mode = parser.ArgAsInt('m');
  }
  gOptDecayMode = (NNBarOscMode_t) mode;
  bool valid_mode = utils::nnbar_osc::IsValidMode(gOptDecayMode);
  if(!valid_mode) {
    LOG("gevgen_nnbar_osc", pFATAL) 
        << "You need to specify a valid annihilation mode";
    PrintSyntax();
    exit(0);
  } //-m

  //
  // geometry
  //

  string geom = "";
  string lunits, dunits;
  if( parser.OptionExists('g') ) {
    LOG("gevgen_nnbar_osc", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');

    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file = 
            ! (gSystem->AccessPathName(geom.c_str()));
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    }                 
  } else {
      LOG("gevgen_nnbar_osc", pFATAL) 
        << "No geometry option specified - Exiting";
      PrintSyntax();
      exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // legth units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_nnbar_osc", pDEBUG) 
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_nnbar_osc", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // density units:
     if( parser.OptionExists('D') ) {
        LOG("gevgen_nnbar_osc", pDEBUG) 
           << "Checking for input geometry density units";
        dunits = parser.ArgAsString('D');
     } else {
        LOG("gevgen_nnbar_osc", pDEBUG) << "Using default geometry density units";
        dunits = kDefOptGeomDUnits;
     } // -D 
     gOptGeomLUnits = utils::units::UnitFromString(lunits);
     gOptGeomDUnits = utils::units::UnitFromString(dunits);

     // check whether an event generation volume name has been 
     // specified -- default is the 'top volume'
     if( parser.OptionExists('t') ) {
        LOG("gevgen_nnbar_osc", pDEBUG) << "Checking for input volume name";
        gOptRootGeomTopVol = parser.ArgAsString('t');
     } else {
        LOG("gevgen_nnbar_osc", pDEBUG) << "Using the <master volume>";
     } // -t 

  } // using root geom?

  else {
    // User has specified a target mix.
    // Decode the list of target pdf codes & their corresponding weight fraction
    // (specified as 'pdg_code_1[fraction_1],pdg_code_2[fraction_2],...')
    // See documentation on top section of this file.
    //
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(geom,",");
    if(tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));    
    } else {
      vector<string>::const_iterator tgtmix_iter = tgtmix.begin();
      for( ; tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
         string tgt_with_wgt = *tgtmix_iter;
         string::size_type open_bracket  = tgt_with_wgt.find("[");
         string::size_type close_bracket = tgt_with_wgt.find("]");
         if (open_bracket ==string::npos || 
             close_bracket==string::npos) 
         {
             LOG("gevgen_nnbar_osc", pFATAL) 
                << "You made an error in specifying the target mix"; 
             PrintSyntax();
             exit(1);
         }
         string::size_type ibeg = 0;
         string::size_type iend = open_bracket;
         string::size_type jbeg = open_bracket+1;
         string::size_type jend = close_bracket;
         int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend-ibeg).c_str());
         double wgt = atof(tgt_with_wgt.substr(jbeg,jend-jbeg).c_str());
         LOG("gevgen_nnbar_osc", pDEBUG) 
            << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));

      }// tgtmix_iter
    } // >1 materials in mix
  } // using tgt mix?

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_nnbar_osc", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_nnbar_osc", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o


  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_nnbar_osc", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_nnbar_osc", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  //
  // >>> print the command line options
  //

  PDGLibrary * pdglib = PDGLibrary::Instance();

  ostringstream gminfo;
  if (gOptUsingRootGeom) {
    gminfo << "Using ROOT geometry - file: " << gOptRootGeom
           << ", top volume: "
           << ((gOptRootGeomTopVol.size()==0) ? "<master volume>" : gOptRootGeomTopVol)
           << ", length  units: " << lunits
           << ", density units: " << dunits;
  } else {
    gminfo << "Using target mix - ";
    map<int,double>::const_iterator iter;
    for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
          int    pdg_code = iter->first;
          double wgt      = iter->second;
          TParticlePDG * p = pdglib->Find(pdg_code);
          if(p) {
            string name = p->GetName();
            gminfo << "(" << name << ") -> " << 100*wgt << "% / ";
          }//p?
    }
  }

  LOG("gevgen_nnbar_osc", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_nosc job configuration");

  LOG("gevgen_nnbar_osc", pNOTICE) 
     << "\n @@ Run number: " << gOptRunNu
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Decay channel $ " << utils::nnbar_osc::AsString(gOptDecayMode)
     << "\n @@ Geometry      $ " << gminfo.str()
     << "\n @@ Statistics    $ " << gOptNev << " events";

  //
  // Temporary warnings...
  //
  if(gOptUsingRootGeom) {
     LOG("gevgen_nnbar_osc", pWARN) 
        << "\n ** ROOT geometries not supported yet in neutron oscillation mode"
        << "\n ** (they will be in the very near future)"
        << "\n ** Please specify a `target mix' instead.";
     gAbortingInErr = true;
     exit(1);
  }
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_nnbar_osc", pFATAL) 
   << "\n **Syntax**"
   << "\n gevgen_nnbarosc [-h] "
   << "\n             [-r run#]"
   << "\n              -m decay_mode"
   << "\n              -g geometry"
   << "\n             [-t top_volume_name_at_geom]"
   << "\n             [-L length_units_at_geom]"
   << "\n             [-D density_units_at_geom]"
   << "\n              -n n_of_events "
   << "\n             [-o output_event_file_prefix]"
   << "\n             [--seed random_number_seed]"
   << "\n             [--message-thresholds xml_file]"
   << "\n             [--event-record-print-level level]"
   << "\n             [--mc-job-status-refresh-rate  rate]"
   << "\n"
   << " Please also read the detailed documentation at http://www.genie-mc.org"
   << " or look at the source code: $GENIE/src/Apps/gNNBarOscEvGen.cxx"
   << "\n";
}
//_________________________________________________________________________________________

