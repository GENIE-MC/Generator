//____________________________________________________________________________
/*!

\program gevgen_hadron

\brief   Generates hadron + nucleus interactions using GENIE's INTRANUKE
         Similar to NEUGEN's pitest (S.Dytman & H.Gallagher)

         Syntax :
           gevgen_hadron [-n nev] -p probe -t tgt [-r run#] -k KE
                         [-f flux] [-o prefix] [-m mode]
                         [--seed random_number_seed]
                         [--message-thresholds xml_file]
                         [--event-record-print-level level]
                         [--mc-job-status-refresh-rate  rate]

         Options :
           [] Denotes an optional argument
           -n
              Specifies the number of events to generate (default: 10000)
           -p
              Specifies the incoming hadron PDG code
           -t
              Specifies the nuclear target PDG code (10LZZZAAAI)
           -r
              Specifies the MC run number (default: 0)
           -k
              Specifies the incoming hadron's kinetic energy (in GeV).
              The same option can be use to specify a kinetic energy range as
              a comma-separated set of numbers (eg 0.1,1.2).
              If no flux is specified then the hadrons will be fired with a
              uniform kinetic energy distribution within the specified range.
              If a kinetic energy spectrum is supplied then the hadrons will be
              fired with a kinetic energy following the input spectrum within
              the specified range.
           -f
              Specifies the incoming hadron's kinetic energy spectrum -
              it can be either a function, eg 'x*x+4*exp(-x)' or a text file
              containing 2 columns corresponding to (kinetic energy {GeV}, 'flux').
           -o
              Output filename prefix
           -m
              INTRANUKE mode <hA, hN> (default: hA)
           --seed
              Random number seed.
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
           --event-record-print-level
              Allows users to set the level of information shown when the event
              record is printed in the screen. See GHepRecord::Print().
           --mc-job-status-refresh-rate
              Allows users to customize the refresh rate of the status file.

         Examples:

         (1) Generate 100k pi^{+}+Fe56 events with a pi^{+} kinetic energy
             of 165 MeV:
             % gevgen_hadron -n 100000 -p 211 -t 1000260560 -k 0.165

         (2) Generate 100k pi^{+}+Fe56 events with the pi^{+} kinetic energy
             distributed uniformly in the [165 MeV, 1200 MeV] range:
             % gevgen_hadron -n 100000 -p 211 -t 1000260560 -k 0.165,1.200

         (3) Generate 100k pi^{+}+Fe56 events with the pi^{+} kinetic energy
             distributed as f(KE) = 1/KE in the [165 MeV, 1200 MeV] range:
             % ghAevgen -n gevgen_hadron -p 211 -t 1000260560 -k 0.165,1.200 -f '1/x'

\authors  Steve Dytman, Minsuk Kim and Aaron Meyer
          University of Pittsburgh

          Costas Andreopoulos,
          University of Liverpool & STFC Rutherford Appleton Lab

\version 1.3

\created May 1, 2007

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <cstdlib>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/CmdLnArgParser.h"

#include "Physics/HadronTransport/INukeHadroFates.h"
#include "Physics/HadronTransport/INukeUtils.h"

using namespace genie;
using namespace genie::controls;

using namespace genie::utils::intranuke;

// Function prototypes
void                        GetCommandLineArgs    (int argc, char ** argv);
const EventRecordVisitorI * GetIntranuke          (void);
double                      GenProbeKineticEnergy (void);
EventRecord *               InitializeEvent       (void);
void                        BuildSpectrum         (void);
void                        PrintSyntax           (void);

// Default options
int     kDefOptNevents      = 10000;   // n-events to generate
Long_t  kDefOptRunNu        = 0;       // default run number
string  kDefOptEvFilePrefix = "gntp.inuke";  // default output file prefix
string  kDefOptMode         = "hA";    // default mode

// User-specified options:
string   gOptMode;             // mode variable
Long_t   gOptRunNu;            // run number
int      gOptNevents;          // n-events to generate
int      gOptProbePdgCode;     // probe  PDG code
int      gOptTgtPdgCode;       // target PDG code
double   gOptProbeKE;          // incoming hadron kinetic enegy (GeV) - for monoenergetic probes
double   gOptProbeKEmin;       // incoming hadron kinetic enegy (GeV) - if using flux
double   gOptProbeKEmax;       // incoming hadron kinetic enegy (GeV) - if using flux
string   gOptFlux;             // input flux (function or flux file)
string   gOptEvFilePrefix;     // event file prefix
bool     gOptUsingFlux=false;  // using kinetic energy distribution?
long int gOptRanSeed ;         // random number seed

TH1D * gSpectrum  = 0;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gmkspl", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // Init random number generator generator with user-specified seed number,
  // set user-specified mesg thresholds, set user-specified GHEP print-level
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  // Build the incident hadron kinetic energy spectrum, if required
  BuildSpectrum();

  LOG("gevgen_hadron",pNOTICE) << "finish setup";

  // Get the specified INTRANUKE model
  const EventRecordVisitorI * intranuke = GetIntranuke();

  // Initialize an Ntuple Writer to save GHEP records into a ROOT tree
  NtpWriter ntpw(kNFGHEP, gOptRunNu);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // Create an MC job monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  LOG("gevgen_hadron",pNOTICE) << "ready to generate events";

  //
  // Generate events
  //

  int ievent = 0;
  while (ievent < gOptNevents) {
      LOG("gevgen_hadron", pNOTICE)
         << " *** Generating event............ " << ievent;

      // initialize
      EventRecord * evrec = InitializeEvent();

      // generate full h+A event
      intranuke->ProcessEventRecord(evrec);

      // print generated events
      LOG("gevgen_hadron", pNOTICE ) << *evrec;

      // add event at the output ntuple
      ntpw.AddEventRecord(ievent, evrec);

      // refresh the mc job monitor
      mcjmonitor.Update(ievent,evrec);

      ievent++;
      delete evrec;

  } // end loop events

  // Save the generated MC events
  ntpw.Save();

  // Clean-up
  if(gSpectrum) {
    delete gSpectrum;
    gSpectrum = 0;
  }

  return 0;
}
//____________________________________________________________________________
const EventRecordVisitorI * GetIntranuke(void)
{
// get the requested INTRANUKE module

  string sname = "";
  string sconf = "";

  if(gOptMode.compare("hA")==0) {
     sname = "genie::HAIntranuke";
     sconf = "Default";
  }
  else
  if(gOptMode.compare("hN")==0) {
     sname = "genie::HNIntranuke";
     sconf = "Default";
  }
  else
  if(gOptMode.compare("hA2019")==0) {
     sname = "genie::HAIntranuke2019";
     sconf = "Default";
  }
  else
  if(gOptMode.compare("hN2019")==0) {
     sname = "genie::HNIntranuke2019";
     sconf = "Default";
  }
  else
  if(gOptMode.compare("hA2018")==0) {
     sname = "genie::HAIntranuke2018";
     sconf = "Default";
  }
  else
  if(gOptMode.compare("hN2018")==0) {
     sname = "genie::HNIntranuke2018";
     sconf = "Default";
  }
  else {
    LOG("gevgen_hadron", pFATAL) << "Invalid Intranuke mode - Exiting";
    gAbortingInErr = true;
    exit(1);
  }

  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * intranuke =
   dynamic_cast<const EventRecordVisitorI *> (algf->GetAlgorithm(sname,sconf));
  assert(intranuke);

  return intranuke;
}
//____________________________________________________________________________
EventRecord * InitializeEvent(void)
{
// Initialize event record. Inserting the probe and target particles.

  EventRecord * evrec = new EventRecord();
  Interaction * interaction = new Interaction;
  evrec->AttachSummary(interaction);

  // dummy vertex position
  TLorentzVector x4null(0.,0.,0.,0.);

  // incident hadron & target nucleon masses
  PDGLibrary * pdglib = PDGLibrary::Instance();
  double mh  = pdglib -> Find (gOptProbePdgCode) -> Mass();
  double M   = pdglib -> Find (gOptTgtPdgCode  ) -> Mass();

  // incident hadron kinetic energy
  double ke = GenProbeKineticEnergy();

  // form  incident hadron and target 4-momenta
  double Eh  = mh + ke;
  double pzh = TMath::Sqrt(TMath::Max(0.,Eh*Eh-mh*mh));
  TLorentzVector p4h   (0.,0.,pzh,Eh);
  TLorentzVector p4tgt (0.,0.,0., M);

  // insert probe and target entries
  GHepStatus_t ist = kIStInitialState;
  evrec->AddParticle(gOptProbePdgCode, ist, -1,-1,-1,-1, p4h,   x4null);
  evrec->AddParticle(gOptTgtPdgCode,   ist, -1,-1,-1,-1, p4tgt, x4null);

  return evrec;
}
//____________________________________________________________________________
double GenProbeKineticEnergy(void)
{
  if(gOptUsingFlux) return gSpectrum->GetRandom(); // spectrum
  else              return gOptProbeKE;            // mono-energetic
}
//____________________________________________________________________________
void BuildSpectrum(void)
{
// Create kinetic energy spectrum from input function
//

  if(!gOptUsingFlux) return;

  if(gSpectrum) {
    delete gSpectrum;
    gSpectrum = 0;
  }

  LOG("gevgen_hadron", pNOTICE)
      << "Generating a flux histogram ... ";

  int    flux_bins    = 300;
  int    flux_entries = 1000000;
  double ke_min       = gOptProbeKEmin;
  double ke_max       = gOptProbeKEmax;
  double dke          = ke_max - ke_min;
  assert(dke>0 && ke_min>=0.);

  // kinetic energy spectrum
  gSpectrum  = new TH1D(
    "spectrum","hadron kinetic energy spectrum", flux_bins, ke_min, ke_max);

  // check whether the input flux is a file or a functional form
  bool input_is_file = ! gSystem->AccessPathName(gOptFlux.c_str());
  if(input_is_file) {
    Spline * input_flux = new Spline(gOptFlux.c_str());

    // generate the flux hisogram from the input flux file
    int    n       = 100;
    double ke_step = (ke_max-ke_min)/(n-1);
    double ymax  = -1, ry = -1, gy = -1, ke = -1;
    for(int i=0; i<n; i++) {
      ke   = ke_min + i*ke_step;
      ymax = TMath::Max(ymax, input_flux->Evaluate(ke));
    }
    ymax *= 1.3;

    RandomGen * r = RandomGen::Instance();

    for(int ientry=0; ientry<flux_entries; ientry++) {
      bool accept = false;
      unsigned int iter=0;
      while(!accept) {
        iter++;
        if(iter > kRjMaxIterations) {
           LOG("gevgen_hadron", pFATAL)
             << "Couldn't generate a flux histogram";
           gAbortingInErr = true;
           exit(1);
        }
        ke = ke_min + dke * r->RndGen().Rndm();
        gy = ymax * r->RndGen().Rndm();
        ry = input_flux->Evaluate(ke);
        accept = gy < ry;
        if(accept) gSpectrum->Fill(ke);
      }
    }
    delete input_flux;

  } else {
    // generate the flux histogram from the input functional form
    TF1 *  input_func = new TF1("input_func", gOptFlux.c_str(), ke_min, ke_max);
    gSpectrum->FillRandom("input_func", flux_entries);
    delete input_func;
  }
  TFile f("./input-hadron-flux.root","recreate");
  gSpectrum->Write();
  f.Close();
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_hadron", pINFO) << "Parsing command line arguments";

  // Common run options.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_hadron", pINFO) << "Reading number of events to generate";
    gOptNevents = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_hadron", pINFO)
       << "Unspecified number of events to generate - Using default";
    gOptNevents = kDefOptNevents;
  }

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevgen_hadron", pINFO) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_hadron", pINFO) << "Unspecified run number - Using default";
    gOptRunNu = kDefOptRunNu;
  }

  // incoming hadron PDG code
  if( parser.OptionExists('p') ) {
    LOG("gevgen_hadron", pINFO) << "Reading rescattering particle PDG code";
    gOptProbePdgCode = parser.ArgAsInt('p');
  } else {
    LOG("gevgen_hadron", pFATAL) << "Unspecified PDG code - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // target PDG code
  if( parser.OptionExists('t') ) {
    LOG("gevgen_hadron", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = parser.ArgAsInt('t');
  } else {
    LOG("gevgen_hadron", pFATAL) << "Unspecified target PDG code - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // target PDG code
  if( parser.OptionExists('m') ) {
    LOG("gevgen_hadron", pINFO) << "Reading mode";
    gOptMode = parser.ArgAsString('m');
  } else {
    LOG("gevgen_hadron", pFATAL) << "Unspecified mode - Using default";
    gOptMode = kDefOptMode;
  }

  // flux functional form or flux file
  if( parser.OptionExists('f') ) {
    LOG("gevgen_hadron", pINFO) << "Reading hadron's kinetic energy spectrum";
    gOptFlux = parser.ArgAsString('f');
    gOptUsingFlux = true;
  }

  // incoming hadron kinetic energy (or kinetic energy range, if using flux)
  if( parser.OptionExists('k') ) {
    LOG("gevgen_hadron", pINFO) << "Reading probe kinetic energy";
    string ke = parser.ArgAsString('k');
    // is it just a value or a range (comma separated set of values)
    if(ke.find(",") != string::npos) {
       // split the comma separated list
       vector<string> kerange = utils::str::Split(ke, ",");
       assert(kerange.size() == 2);
       double kemin = atof(kerange[0].c_str());
       double kemax = atof(kerange[1].c_str());
       assert(kemax>kemin && kemin>0);
       gOptProbeKE    = -1;
       gOptProbeKEmin = kemin;
       gOptProbeKEmax = kemax;
       // if no flux was specified, generate uniformly within given range
       if(!gOptUsingFlux) {
          gOptFlux = "1";
          gOptUsingFlux = true;
       }
    } else {
       gOptProbeKE    = atof(ke.c_str());
       gOptProbeKEmin = -1;
       gOptProbeKEmax = -1;
       if(gOptUsingFlux) {
          LOG("gevgen_hadron", pFATAL)
            << "You specified an input flux without a kinetic energy range";
          PrintSyntax();
          gAbortingInErr = true;
          exit(1);
       }
    }
  } else {
    LOG("gevgen_hadron", pFATAL) << "Unspecified kinetic energy - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_hadron", pINFO) << "Reading the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_hadron", pDEBUG)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // INTRANUKE mode
  if( parser.OptionExists('m') ) {
    LOG("gevgen_hadron", pINFO) << "Reading mode";
    gOptMode = parser.ArgAsString('m');
  } else {
    LOG("gevgen_hadron", pDEBUG)
       << "Unspecified mode - Using default";
    gOptMode = kDefOptMode;
  }

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_hadron", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_hadron", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }


  LOG("gevgen_hadron", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gevgen_hadron job configuration");

  LOG("gevgen_hadron", pNOTICE) << "MC Run Number      = " << gOptRunNu;
  LOG("gevgen_hadron", pNOTICE) << "Random number seed = " << gOptRanSeed;
  LOG("gevgen_hadron", pNOTICE) << "Mode               = " << gOptMode;
  LOG("gevgen_hadron", pNOTICE) << "Number of events   = " << gOptNevents;
  LOG("gevgen_hadron", pNOTICE) << "Probe PDG code     = " << gOptProbePdgCode;
  LOG("gevgen_hadron", pNOTICE) << "Target PDG code    = " << gOptTgtPdgCode;
  if(gOptProbeKEmin<0 && gOptProbeKEmax<0) {
    LOG("gevgen_hadron", pNOTICE)
        << "Hadron input KE    = " << gOptProbeKE;
  } else {
    LOG("gevgen_hadron", pNOTICE)
        << "Hadron input KE range = ["
        << gOptProbeKEmin << ", " << gOptProbeKEmax << "]";
  }
  if(gOptUsingFlux) {
    LOG("gevgen_hadron", pNOTICE)
        << "Input flux            = "
        << gOptFlux;
  }

  LOG("gevgen_hadron", pNOTICE) << "\n";
  LOG("gevgen_hadron", pNOTICE) << *RunOpt::Instance();
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_hadron", pNOTICE)
    << "\n\n"
    << "Syntax:" << "\n"
    << "   gevgen_hadron [-r run] [-n nev] -p hadron_pdg -t tgt_pdg -k KE [-m mode] "
    << "                 [-f flux] "
    << "                 [--seed random_number_seed]"
    << "                 [--message-thresholds xml_file]"
    << "                 [--event-record-print-level level]"
    << "                 [--mc-job-status-refresh-rate rate]"
    << "\n";
}
//____________________________________________________________________________
