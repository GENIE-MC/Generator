//________________________________________________________________________________________
/*!

\program gevgen_ndcy

\brief   A GENIE-based nucleon decay event generation application.

         *** Synopsis :

         gevgen_ndcy [-h] 
                     [-r run#] 
                      -n n_of_events
                      -m decay_mode
		     [-N decayed_nucleon_pdg]
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
              Prints out the gevgen_ndcy syntax and exits.
           -r 
              Specifies the MC run number [default: 1000].
           -n  
              Specifies how many events to generate.
           -m 
              Nucleon decay mode ID:
	      see http://www-pdg.lbl.gov/2016/listings/rpp2016-list-p.pdf
	      for nucleon decay mode numbering convention
	      Example:
	      m = 1: N -> e+ pi
	      ...
	      m = 60: n -> 5nu
            

	   -N
	      Decayed nucleon PDG code.
	      Either N=2212 (proton) or N=2112 (neutron) 
	      Example:
	      m = 1 and N = 2112: n -> e+ pi-
	      m = 1 and N = 2212: p -> e+ pi0

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
#include "Physics/NucleonDecay/NucleonDecayMode.h"
#include "Physics/NucleonDecay/NucleonDecayUtils.h"
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
int   SelectInitState    (void);
const EventRecordVisitorI * NucleonDecayGenerator(void);

//
string          kDefOptGeomLUnits   = "mm";    // default geometry length units
string          kDefOptGeomDUnits   = "g_cm3"; // default geometry density units
NtpMCFormat_t   kDefOptNtpFormat    = kNFGHEP; // default event tree format   
string          kDefOptEvFilePrefix = "gntp";

//
Long_t             gOptRunNu        = 1000;                // run number
int                gOptNev          = 10;                  // number of events to generate
NucleonDecayMode_t gOptDecayMode    = kNDNull;             // nucleon decay mode
int                gOptDecayedNucleon = kNDNull;           // decayed nucleon PDG
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
  const EventRecordVisitorI * mcgen = NucleonDecayGenerator();

  // Event loop
  int ievent = 0;
  int dpdg = 0;
  if (gOptDecayedNucleon > 0) {
    dpdg = gOptDecayedNucleon;
  } else {
    dpdg = utils::nucleon_decay::DecayedNucleonPdgCode(gOptDecayMode);
  }  

  while (1)
  {
     if(ievent == gOptNev) break;

     LOG("gevgen_ndcy", pNOTICE)
          << " *** Generating event............ " << ievent;

     EventRecord * event = new EventRecord;
     int target = SelectInitState();
     int decay  = (int)gOptDecayMode;
     Interaction * interaction = Interaction::NDecay(target,decay,dpdg);
     event->AttachSummary(interaction);

     // Simulate decay     
     mcgen->ProcessEventRecord(event);

     LOG("gevgen_ndcy", pINFO)
         << "Generated event: " << *event;

     // Add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     delete event;

     ievent++;
  } // event loop

  // Save the generated event tree & close the output file
  ntpw.Save();

  LOG("gevgen_ndcy", pNOTICE) << "Done!";

  return 0;
}
//_________________________________________________________________________________________
int SelectInitState(void)
{
  if (!utils::nucleon_decay::IsValidMode(gOptDecayMode, gOptDecayedNucleon)) {
    LOG("gevgen_ndcy", pFATAL) << "Not a valid decay mode and/or decayed nucleon...";
    gAbortingInErr = true;
    exit(1);
  }

  int dpdg = 0;
  if (gOptDecayedNucleon > 0) {
    dpdg = gOptDecayedNucleon;
  } else {
     dpdg = utils::nucleon_decay::DecayedNucleonPdgCode(gOptDecayMode);
  }

  map<int,double> cprob; // cumulative probability 
  map<int,double>::const_iterator iter;
 
  double sum_prob = 0;
  for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
     int pdg_code = iter->first;
     int A = pdg::IonPdgCodeToA(pdg_code);
     int Z = pdg::IonPdgCodeToZ(pdg_code);

     double nucleon_decay_fraction = 0.;
     if (dpdg == kPdgProton ) { nucleon_decay_fraction = (double)Z / (double)A; }
     else if (dpdg == kPdgNeutron ) { nucleon_decay_fraction = (double)(A-Z) / (double)A; }

     double wgt  = iter->second;
     double prob = wgt*nucleon_decay_fraction;

     sum_prob += prob;
     cprob.insert(map<int, double>::value_type(pdg_code, sum_prob));
  }

  assert(sum_prob > 0.);

  RandomGen * rnd = RandomGen::Instance();
  double r = sum_prob * rnd->RndEvg().Rndm();

  for(iter = cprob.begin(); iter != cprob.end(); ++iter) {
     int pdg_code = iter->first;
     double prob  = iter->second;
     if(r < prob) {
       LOG("gevgen_ndcy", pNOTICE) << "Selected initial state = " << pdg_code;
       return pdg_code;
     }
  }  

  LOG("gevgen_ndcy", pFATAL) << "Couldn't select an initial state...";
  gAbortingInErr = true;
  exit(1);
}
//_________________________________________________________________________________________
const EventRecordVisitorI * NucleonDecayGenerator(void)
{
  string sname   = "genie::EventGenerator";
  string sconfig = "NucleonDecay";
  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * mcgen =
     dynamic_cast<const EventRecordVisitorI *> (algf->GetAlgorithm(sname,sconfig));
  if(!mcgen) {
     LOG("gevgen_ndcy", pFATAL) << "Couldn't instantiate the nucleon decay generator";
     gAbortingInErr = true;
     exit(1);
  }
  return mcgen;
}
//_________________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_ndcy", pINFO) << "Parsing command line arguments";

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
    LOG("gevgen_ndcy", pDEBUG) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_ndcy", pDEBUG) << "Unspecified run number - Using default";
    gOptRunNu = 1000;
  } //-r


  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_ndcy", pDEBUG) 
        << "Reading number of events to generate";
    gOptNev = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_ndcy", pFATAL) 
        << "You need to specify the number of events";
    PrintSyntax();
    exit(0);
  } //-n

  // decay mode
  int mode = -1;
  if( parser.OptionExists('m') ) {
    LOG("gevgen_ndcy", pDEBUG) 
        << "Reading decay mode";
    mode = parser.ArgAsInt('m');
  } else {
    LOG("gevgen_ndcy", pFATAL) 
        << "You need to specify the decay mode";
    PrintSyntax();
    exit(0);
  } //-m
  gOptDecayMode = (NucleonDecayMode_t) mode;

  // decayed nucleon PDG
  if( parser.OptionExists('N') ) {
    LOG("gevgen_ndcy", pINFO) << "Reading decayed nucleon PDG";
    gOptDecayedNucleon = parser.ArgAsInt('N');
  } else {
    LOG("gevgen_ndcy", pINFO) << "Unspecified decayed nucleon PDG - Using default";
    gOptDecayedNucleon = 0;
  }  

  bool valid_mode = utils::nucleon_decay::IsValidMode(gOptDecayMode, gOptDecayedNucleon);
  if(!valid_mode) {
    LOG("gevgen_ndcy", pFATAL) 
        << "You need to specify a valid decay mode / decayed nucleon PDG combination";
    PrintSyntax();
    exit(0);
  }

  //
  // geometry
  //

  string geom = "";
  string lunits, dunits;
  if( parser.OptionExists('g') ) {
    LOG("gevgen_ndcy", pDEBUG) << "Getting input geometry";
    geom = parser.ArgAsString('g');

    // is it a ROOT file that contains a ROOT geometry?
    bool accessible_geom_file = 
            ! (gSystem->AccessPathName(geom.c_str()));
    if (accessible_geom_file) {
      gOptRootGeom      = geom;
      gOptUsingRootGeom = true;
    }                 
  } else {
      LOG("gevgen_ndcy", pFATAL) 
        << "No geometry option specified - Exiting";
      PrintSyntax();
      exit(1);
  } //-g

  if(gOptUsingRootGeom) {
     // using a ROOT geometry - get requested geometry units

     // legth units:
     if( parser.OptionExists('L') ) {
        LOG("gevgen_ndcy", pDEBUG) 
           << "Checking for input geometry length units";
        lunits = parser.ArgAsString('L');
     } else {
        LOG("gevgen_ndcy", pDEBUG) << "Using default geometry length units";
        lunits = kDefOptGeomLUnits;
     } // -L
     // density units:
     if( parser.OptionExists('D') ) {
        LOG("gevgen_ndcy", pDEBUG) 
           << "Checking for input geometry density units";
        dunits = parser.ArgAsString('D');
     } else {
        LOG("gevgen_ndcy", pDEBUG) << "Using default geometry density units";
        dunits = kDefOptGeomDUnits;
     } // -D 
     gOptGeomLUnits = utils::units::UnitFromString(lunits);
     gOptGeomDUnits = utils::units::UnitFromString(dunits);

     // check whether an event generation volume name has been 
     // specified -- default is the 'top volume'
     if( parser.OptionExists('t') ) {
        LOG("gevgen_ndcy", pDEBUG) << "Checking for input volume name";
        gOptRootGeomTopVol = parser.ArgAsString('t');
     } else {
        LOG("gevgen_ndcy", pDEBUG) << "Using the <master volume>";
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
             LOG("gevgen_ndcy", pFATAL) 
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
         LOG("gevgen_ndcy", pDEBUG) 
            << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));

      }// tgtmix_iter
    } // >1 materials in mix
  } // using tgt mix?

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_ndcy", pDEBUG) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_ndcy", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o


  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_ndcy", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_ndcy", pINFO) << "Unspecified random number seed - Using default";
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

  LOG("gevgen_ndcy", pNOTICE)
     << "\n\n"
     << utils::print::PrintFramedMesg("gevgen_ndcy job configuration");

  LOG("gevgen_ndcy", pNOTICE) 
     << "\n @@ Run number: " << gOptRunNu
     << "\n @@ Random number seed: " << gOptRanSeed
     << "\n @@ Decay channel $ " << utils::nucleon_decay::AsString(gOptDecayMode, gOptDecayedNucleon)
     << "\n @@ Geometry      $ " << gminfo.str()
     << "\n @@ Statistics    $ " << gOptNev << " events";

  //
  // Temporary warnings...
  //
  if(gOptUsingRootGeom) {
     LOG("gevgen_ndcy", pWARN) 
        << "\n ** ROOT geometries not supported yet in the nucleon decay mode"
        << "\n ** (they will be in the very near future)"
        << "\n ** Please specify a `target mix' instead.";
     gAbortingInErr = true;
     exit(1);
  }
}
//_________________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_ndcy", pFATAL) 
   << "\n **Syntax**"
   << "\n gevgen_ndcy [-h] "
   << "\n             [-r run#]"
   << "\n              -m decay_mode"
   << "\n             [-N decayed_nucleon]"
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
   << " or look at the source code: $GENIE/src/Apps/gNucleonDecayEvGen.cxx"
   << "\n";
}
//_________________________________________________________________________________________

