//____________________________________________________________________________
/*!

\program gevgen_hadron

\brief   Generates hadron + nucleus interactions using GENIE's INTRANUKE
	 Similar to NEUGEN's pitest (S.Dytman & H.Gallagher)

         Syntax :
           gevgen_hadron [-n nev] -p probe -t tgt [-r run#] -k KE 
                    [-f flux] [-o prefix] [-m mode]

         Options :
           [] Denotes an optional argument
           -n Specifies the number of events to generate (default: 10000)
           -p Specifies the incoming hadron PDG code 
           -t Specifies the nuclear target PDG code (10LZZZAAAI)
           -r Specifies the MC run number (default: 0)
           -k Specifies the incoming hadron's kinetic energy (in GeV).
              The same option can be use to specify a kinetic energy range as
	      a comma-separated set of numbers (eg 0.1,1.2).
              If no flux is specified then the hadrons will be fired with a 
              uniform kinetic energy distribution within the specified range.
	      If a kinetic energy spectrum is supplied then the hadrons will be
              fired with a kinetic energy following the input spectrum within 
              the specified range.
           -f Specifies the incoming hadron's kinetic energy spectrum - 
              it can be either a function, eg 'x*x+4*exp(-x)' or a text file 
              containing 2 columns corresponding to (kinetic energy {GeV}, 'flux').
           -o output filename prefix
           -m INTRANUKE mode <hA, hN> (default: hA)
	      
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
          STFC, Rutherford Appleton Laboratory
	  
\version 1.3

\created May 1, 2007

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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

#include "Algorithm/AlgFactory.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepStatus.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLnArgParser.h"

using namespace genie;
using namespace genie::controls;

// Function prototypes
void BuildKineticEnergySpectrum (void);

void                        GetCommandLineArgs    (int argc, char ** argv);
const EventRecordVisitorI * GetIntranuke          (void);
double                      GenProbeKineticEnergy (void);
EventRecord *               InitializeEvent       (void);
void                        PrintSyntax           (void);

// Default options 
int     kDefOptNevents      = 10000;   // n-events to generate
Long_t  kDefOptRunNu        = 0;       // default run number
string  kDefOptEvFilePrefix = "gntp";

// User-specified options:
string gOptMode;             // mode variable
Long_t gOptRunNu;            // run number
int    gOptNevents;          // n-events to generate
int    gOptProbePdgCode;     // probe  PDG code
int    gOptTgtPdgCode;       // target PDG code
double gOptProbeKE;          // incoming hadron kinetic enegy (GeV) - for monoenergetic probes
double gOptProbeKEmin;       // incoming hadron kinetic enegy (GeV) - if using flux
double gOptProbeKEmax;       // incoming hadron kinetic enegy (GeV) - if using flux
string gOptFlux;             // input flux (function or flux file)
string gOptEvFilePrefix;     // event file prefix
bool   gOptUsingFlux=false;  // using kinetic energy distribution?

TH1D * gSpectrum  = 0;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // parse command line arguments
  GetCommandLineArgs(argc,argv);

  // build the incident hadron kinetic energy spectrum, if required
  BuildKineticEnergySpectrum();

  // get the specified INTRANUKE model
  const EventRecordVisitorI * intranuke = GetIntranuke();

  // initialize an Ntuple Writer to save GHEP records into a ROOT tree
  NtpWriter ntpw(kNFGHEP, gOptRunNu);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // create an MC job monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //
  // generate events 
  //

  int ievent = 0;
  while (ievent < gOptNevents) {
      LOG("gevgen_hadron", pINFO) 
         << " *** Generating event............ " << ievent;
      
      // initialize
      EventRecord * evrec = InitializeEvent();

      // generate full h+A event
      intranuke->ProcessEventRecord(evrec);
  
      // print n first generated events (then continue printing out 
      // with debug priority level)
      if(ievent > 100) { LOG("gevgen_hadron", pDEBUG ) << *evrec; }
      else             { LOG("gevgen_hadron", pNOTICE) << *evrec; }
  
      // add event at the output ntuple
      ntpw.AddEventRecord(ievent, evrec);
      
      // refresh the mc job monitor
      mcjmonitor.Update(ievent,evrec);
      
      ievent++;
      delete evrec;
      
  } // end loop events
  
  //-- save the generated MC events
  ntpw.Save();

  //-- clean-up
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

  string sname   = "";
  string sconfig = "";

  // comment-out this block after Aaron's and Steve's intranuke upgrade
  sname   = "genie::Intranuke"; 
  sconfig = "hA";  

/*
  // uncomment this block after Aaron's and Steve's intranuke upgrade
  
  if(gOptMode.compare("hA")==0) {
     sname   = "genie::HAIntranuke"; 
     sconfig = "hA";
  }
  else 
  if(gOptMode.compare("hN")==0) {
     sname   = "genie::HNIntranuke"; 
     sconfig = "hN";
  }
  else {
    LOG("gevgen_hadron", pFATAL) << "Invalid Intranuke mode - Exiting";
    gAbortingInErr = true;
    exit(1);
  }
*/

  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * intranuke = 
            dynamic_cast<const EventRecordVisitorI *> (
                     algf->GetAlgorithm(sname,sconfig + "-TestMode"));

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
void BuildKineticEnergySpectrum(void)
{
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
  LOG("gevgen_hadron", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // number of events:
  if( parser.OptionExists('n') ) {
    LOG("gevgen_hadron", pINFO) << "Reading number of events to generate";
    gOptNevents = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_hadron", pINFO)
       << "Unspecified number of events to generate - Using default";
    gOptNevents = kDefOptNevents;
  }

  // run number:
  if( parser.OptionExists('r') ) {
    LOG("gevgen_hadron", pINFO) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_hadron", pINFO) << "Unspecified run number - Using default";
    gOptRunNu = kDefOptRunNu;
  }

  // incoming hadron PDG code:
  if( parser.OptionExists('p') ) {
    LOG("gevgen_hadron", pINFO) << "Reading rescattering particle PDG code";
    gOptProbePdgCode = parser.ArgAsInt('p');
  } else {
    LOG("gevgen_hadron", pFATAL) << "Unspecified PDG code - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // target PDG code:
  if( parser.OptionExists('t') ) {
    LOG("gevgen_hadron", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = parser.ArgAsInt('t');
  } else {
    LOG("gevgen_hadron", pFATAL) << "Unspecified target PDG code - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // flux functional form or flux file
  if( parser.OptionExists('f') ) {
    LOG("gevgen_hadron", pINFO) << "Reading hadron's kinetic energy spectrum";
    gOptFlux = parser.ArgAsString('f');
    gOptUsingFlux = true;
  }

  // incoming hadron kinetic energy (or kinetic energy range, if using flux):
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

  LOG("gevgen_hadron", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("gevgen_hadron", pINFO) << "MC Run Number              = " << gOptRunNu;
  LOG("gevgen_hadron", pINFO) << "Probe PDG code             = " << gOptProbePdgCode;
  LOG("gevgen_hadron", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  if(gOptProbeKEmin<0 && gOptProbeKEmax<0) {
    LOG("gevgen_hadron", pINFO) 
        << "Hadron input KE            = " << gOptProbeKE;
  } else {
    LOG("gevgen_hadron", pINFO) 
        << "Hadron input KE range      = [" 
        << gOptProbeKEmin << ", " << gOptProbeKEmax << "]";
  }
  if(gOptUsingFlux) {
    LOG("gevgen_hadron", pINFO) 
        << "Input flux                 = " 
        << gOptFlux;
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_hadron", pNOTICE)
    << "\n\n" 
    << "Syntax:" << "\n"
    << "   ghAevgen [-n nev] -p hadron_pdg -t tgt_pdg [-r run] "
    << "             -k KE [-f flux] [-m mode]"
    << "\n";
}
//____________________________________________________________________________
