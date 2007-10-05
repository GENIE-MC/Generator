//____________________________________________________________________________
/*!

\program ghAevgen

\brief   Generates hadron + nucleus interactions using GENIE's INTRANUKE
	 Similar to NEUGEN's pitest (S.Dytman & H.Gallagher)

         Syntax :
           ghAevgen [-n nev] -p hadron_pdg -t tgt_pdg [-r run#] -k KE [-f flux]

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
	      
         Examples:

           shell% ghAevgen -n 100000 -p 211 -t 1000260560 -k 0.165
           * will generate 100k pi^{+}+Fe56 events with a pi^{+} kinetic
             energy of 165 MeV

           shell% ghAevgen -n 100000 -p 211 -t 1000260560 -k 0.165,1.200
           * will generate 100k pi^{+}+Fe56 events with the pi^{+} kinetic
             energy distributed uniformly in the [165 MeV, 1200 MeV] range

           shell% ghAevgen -n 100000 -p 211 -t 1000260560 -k 0.165,1.200 -f '1/x'
           * will generate 100k pi^{+}+Fe56 events with the pi^{+} kinetic
             energy distributed as f(KE) = 1/KE in the [165 MeV, 1200 MeV] range

\authors  Minsuk Kim and Steve Dytman
          University of Pittsburgh

	  Costas Andreopoulos,
          STFC, Rutherford Appleton Laboratory
	  
\version 1.2

\created May 1, 2007

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>

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
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using namespace genie;
using namespace genie::controls;

void GetCommandLineArgs         (int argc, char ** argv);
void BuildKineticEnergySpectrum (void);
void PrintSyntax                (void);

//Default options 
int     kDefOptNevents   = 10000;   // n-events to generate
Long_t  kDefOptRunNu     = 0;       // default run number

//User-specified options:
int           gOptNevents;          // n-events to generate
int           gOptTgtPdgCode;       // target PDG code
Long_t        gOptRunNu;            // run number
int           gOptInpHadPdgCode;    // incoming hadron PDG code
double        gOptInpHadKE;         // incoming hadron kinetic enegy (GeV) or min energy in specified range
double        gOptInpHadDeltaKE;    // incoming hadron kinetic energy range (GeV)
string        gOptFlux;             // input flux (function or flux file)

// globals
bool          gOptUsingFlux = false;
TH1D *        gSpectrum     = 0;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- build the incident hadron kinetic energy spectrum, if required
  if(gOptUsingFlux) 
	BuildKineticEnergySpectrum();

  //-- get the intranuke module
  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * intranuke = 
            dynamic_cast<const EventRecordVisitorI *> (
                     algf->GetAlgorithm("genie::Intranuke","hA-TestMode"));

  //-- initialize an Ntuple Writer to save GHEP records into a ROOT tree
  NtpWriter ntpw(kNFEventRecord, gOptRunNu);
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- get the pdg library
  PDGLibrary * pdglib = PDGLibrary::Instance();

  //-- dummy vertex position
  TLorentzVector x4null(0.,0.,0.,0.);

  //-- incident hadron & target nucleon masses
  double mh  = pdglib->Find(gOptInpHadPdgCode)->Mass();
  double M   = pdglib->Find(gOptTgtPdgCode)->Mass();

  //-- generate events / print the GHEP record / add it to the event tree
  int ievent = 0;
  while (ievent < gOptNevents) {
      LOG("ghAevgen", pINFO) << " *** Generating event............ " << ievent;
      
      //-- generate a kinetic energy & form the input 4-momenta
      double ke = 0;
      if(gOptUsingFlux) {
        ke = gSpectrum->GetRandom();
      } else {
        if(gOptInpHadDeltaKE<0) { ke = gOptInpHadKE; }
        else                    { ke = gOptInpHadKE + rnd->RndGen().Rndm() * gOptInpHadDeltaKE; }
      }

      double Eh  = mh + ke;
      double pzh = TMath::Sqrt(TMath::Max(0.,Eh*Eh-mh*mh));
      TLorentzVector p4h   (0.,0.,pzh,Eh);
      TLorentzVector p4tgt (0.,0.,0.,M);

      // insert initial state
      EventRecord * evrec = new EventRecord();
      Interaction * interaction = new Interaction;
      evrec->AttachSummary(interaction);
            
      evrec->AddParticle(gOptInpHadPdgCode,kIStInitialState,-1,-1,-1,-1,p4h,  x4null);
      evrec->AddParticle(gOptTgtPdgCode,   kIStInitialState,-1,-1,-1,-1,p4tgt,x4null);
      
      // generate h+A eventw
      intranuke->ProcessEventRecord(evrec);
  
      // print generated event    
      LOG("ghAevgen", pINFO) << *evrec;
  
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
void BuildKineticEnergySpectrum(void)
{
  if(gSpectrum) {
    delete gSpectrum;
    gSpectrum = 0;
  }

  int    flux_bins    = 300;
  int    flux_entries = 100000;
  double ke_min       = gOptInpHadKE; 
  double ke_max       = gOptInpHadKE + gOptInpHadDeltaKE; 
  double dke          = gOptInpHadDeltaKE; 

  gSpectrum  = new TH1D(
    "spectrum","hadron kinetic energy spectrum", flux_bins, ke_min, ke_max);

  // check whether the input flux is a file or a functional form
  //
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
           LOG("ghAevgen", pFATAL) << "Couldn't generate a flux histogram";
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
  LOG("ghAevgen", pNOTICE) << "Parsing command line arguments";

  //number of events:
  try {
    LOG("ghAevgen", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAevgen", pINFO)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }

  //run number:
  try {
    LOG("ghAevgen", pINFO) << "Reading MC run number";
    gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAevgen", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }

  // incoming hadron PDG code:
  try {
    LOG("ghAevgen", pINFO) << "Reading rescattering particle PDG code";
    gOptInpHadPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
      if(!e.ArgumentFound()) {
      LOG("ghAevgen", pFATAL) << "Unspecified PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //target PDG code:
  try {
    LOG("ghevAgen", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAevgen", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // incoming hadron kinetic energy:
  try {
    LOG("ghAevgen", pINFO) << "Reading rescattering particle KE energy";
    string ke = genie::utils::clap::CmdLineArgAsString(argc,argv,'k');

    // is it just a value or a range (comma separated set of values)   
    if(ke.find(",") != string::npos) {
       // split the comma separated list
       vector<string> kerange = utils::str::Split(ke, ",");
       assert(kerange.size() == 2);
       double kemin = atof(kerange[0].c_str());
       double kemax = atof(kerange[1].c_str());
       assert(kemax>kemin && kemin>0);
       gOptInpHadKE      = kemin; 
       gOptInpHadDeltaKE = kemax-kemin; 
    } else {
       gOptInpHadKE      = atof(ke.c_str());
       gOptInpHadDeltaKE = -1;
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("ghAevgen", pFATAL) << "Unspecified KE - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // flux functional form or flux file
  try {
    LOG("ghAevgen", pINFO) << "Reading hadron's kinetic energy spectrum";
    gOptFlux = genie::utils::clap::CmdLineArgAsString(argc,argv,'f');
    gOptUsingFlux = true;
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
    }
  }


  LOG("ghAevgen", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("ghAevgen", pINFO) << "MC Run Number              = " << gOptRunNu;
  LOG("ghAevgen", pINFO) << "Incoming hadron PDG code   = " << gOptInpHadPdgCode;
  LOG("ghAevgen", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  if(gOptInpHadDeltaKE<0) {
    LOG("ghAevgen", pINFO) 
        << "Hadron input KE            = " << gOptInpHadKE;
  } else {
    LOG("ghAevgen", pINFO) 
        << "Hadron input KE range      = [" 
        << gOptInpHadKE << ", " << gOptInpHadKE+gOptInpHadDeltaKE << "]";
  }
  if(gOptUsingFlux) {
    LOG("ghAevgen", pINFO) << "Input flux                 = " 
                           << gOptFlux;
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("ghAevgen", pNOTICE)
    << "\n\n" 
    << "Syntax:" << "\n"
    << "   ghAevgen [-n nev] -p hadron_pdg -t tgt_pdg [-r run] "
    << "[-a R0] -k KE [-f flux]"
    << "\n";
}
//____________________________________________________________________________
