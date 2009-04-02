//____________________________________________________________________________
/*!

\program testIntranuke

\brief   Generates hadron + nucleus interactions using GENIE's INTRANUKE
	 Similar to NEUGEN's pitest (S.Dytman & H.Gallagher)

         Syntax :
           gtestIntranuke [-n nev] -p hadron_pdg -t tgt_pdg [-r run#] -k KE -o output_file [-m mode]

         Options :
           [] denotes an optional argument
           -n specifies the number of events to generate (default: 10000)
           -p specifies the incoming hadron PDG code 
           -t specifies the nuclear target PDG code (10LZZZAAAI)
           -r specifies the MC run number (default: 0)
           -k specifies the incoming hadron's kinetic energy (in GeV)
           -o specifies the output file for the x-sec data
           -m specifies the Intranuke mode (this functionality is not working yet)

         Example:
           gtestIntranuke -n 100000 -p 211 -t 1000260560 -k 0.165 -o testIntranukeXsec
           will generate 100k pi^{+} + Fe56 events at E(pi^{+})=165 MeV, and
             write the x-sec data output to a file named testIntranukeXsec
             (in the directory from which gtestIntranuke was run)

\authors  Minsuk Kim and Steve Dytman
          University of Pittsburgh

	  Costas Andreopoulos,
          STFC, Rutherford Appleton Laboratory
	  
\version 1.6

\created January 17, 2008

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <time.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "EVGCore/EventRecordVisitorI.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepStatus.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/Intranuke.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Registry/Registry.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
#include "Utils/INukeUtils.h"

using std::endl;
using std::setw;
using std::setfill;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//Default options 
int     kDefOptNevents   = 10000;   // n-events to generate
Long_t  kDefOptRunNu     = 0;       // default run number
int           kDefOptTgtPdgCode   = 1000260560;
double        kDefOptR0           = 1.4;
int           kDefOptNstep        = 1; 
double        kDefOptLstep        = 0.2;        // some amount of KE (GeV)
int           gOptInputPdgCode    = 211;      // probe particle PDG code
string        kDefOptMode           ("hN");   // default mode, not working yet
string        kDefOptXsecFName      ("gtestIntranuke_xsec");  // default xsec output  filename 
//User-specified options:
int           gOptNevents;          // n-events to generate
int           gOptTgtPdgCode;       // target PDG code
Long_t        gOptRunNu;            // run number
int           gOptInpHadPdgCode;    // incoming hadron PDG code
double        gOptInpHadKE;         // incoming hadron kinetic enegy (GeV)
double        gOptR0;               // nuclear radius param 
string        gOptXsecFName;        // filename for output x-secs
string        gOptMode;             // mode variable, not working yet
ofstream      gOptXsecFile;         // output stream
bool          gOptDoOutput = true;  // whether or not xsec is output
//Analysis variables
bool   isPiProbe = false;  //defines flavor of probe
bool   isNuclProbe = false;  //defines flavor of probe
bool   isKProbe = false;  //defines flavor of probe (future)
bool   isGamProbe = false;  //defines flavor of probe (future)
int    fate = 0;  // fate for each event
const int    nfates = 9;  //total number of possible fates
int    countfate[nfates]; // total no. of events with given fate
double    sigma[nfates]; // cross sections
double    sigma_err[nfates]; // cross section errors
string    fatestr[nfates] = "  ";
INukeFateHA_t fatetype[nfates];
double    nuclear_radius = 0., area = 0.; //target params
int       A, Z;  // target params
clock_t   start,end;    // time variables
//____________________________________________________________________________
int main(int argc, char ** argv)
{
//-- test logger
LOG("testInstranuke", pFATAL) << "Logger is working fine.";

//-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- set seed for random number generator
  RandomGen * r = RandomGen::Instance();
  time_t t; (void) time(&t);
  r->SetSeed(t);

  //-- set output stream flags
  if(gOptDoOutput)
  {
    gOptXsecFile.exceptions(ofstream::failbit | ofstream::eofbit | ofstream::badbit);
  }

  AlgFactory * algf = AlgFactory::Instance();
  const EventRecordVisitorI * intranuke = 
            dynamic_cast<const EventRecordVisitorI *> (
                     algf->GetAlgorithm("genie::Intranuke","hA-TestMode"));

  //-- Code to implement mode; not working yet
  if(gOptMode.compare("hA")==0)
  {
    AlgConfigPool * confpool = AlgConfigPool::Instance();
    Registry * config = confpool->FindRegistry(intranuke);

    config->UnLock();
    config->Set("mode","hA");
  }

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(kNFGHEP, gOptRunNu);   //kNFEventRecord,
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- get the pdg library
  PDGLibrary * pdglib = PDGLibrary::Instance();

  // initialize
  fatetype[0]= kIHAFtUndefined;   
  fatetype[1]= kIHAFtNoInteraction;
  fatetype[2]= kIHAFtCEx;
  fatetype[3]= kIHAFtElas;
  fatetype[4]= kIHAFtInelas;
  fatetype[5]= kIHAFtAbs;
  fatetype[6]= kIHAFtKo;
  fatetype[7]= kIHAFtPiProd;
  fatetype[8]= kIHAFtDCEx;
 
  for (int k=0; k<nfates; k++) {
    countfate[k] = 0; 
    sigma[k] = 0.; 
    sigma_err[k] = 0.;
    fatestr[k] = INukeHadroFates::AsString(fatetype[k]);
  }

  //-- input 4-momenta & dummy vtx
  double mh  = pdglib->Find(gOptInpHadPdgCode)->Mass();
  double M   = pdglib->Find(gOptTgtPdgCode)->Mass();
  int iprobe =  gOptInpHadPdgCode;
  int itarget = gOptTgtPdgCode;
  string probe_name = "", tgt_name = "";
  probe_name = PDGLibrary::Instance()->Find(iprobe)->GetName();
  tgt_name = PDGLibrary::Instance()->Find(itarget)->GetName();
  A = pdg::IonPdgCodeToA(itarget);
  Z = pdg::IonPdgCodeToZ(itarget);
  if (fabs(iprobe)==211 || iprobe==111) isPiProbe = true;
  else if (iprobe==2112 || iprobe==2212) isNuclProbe = true;
  else LOG("testIntranuke", pWARN) << " probe not implemented yet"; //isKProbe = true;
  gOptR0 = kDefOptR0;
  nuclear_radius = 3.0* gOptR0*TMath::Power(A,1./3.);
  area = TMath::Pi()*TMath::Power(nuclear_radius,2);

  double Eh  = mh + gOptInpHadKE;
  double pzh = TMath::Sqrt(TMath::Max(0.,Eh*Eh-mh*mh));

  TLorentzVector p4h   (0.,0.,pzh,Eh);
  TLorentzVector p4tgt (0.,0.,0.,M);
  TLorentzVector x4null(0.,0.,0.,0.);

  //-- generate events / print the GHEP record / add it to the event tree
  int ievent = 0;
  start = clock(); // starting time
  while (ievent < gOptNevents) {
      LOG("testIntranuke", pDEBUG) << " *** Generating event............ " << ievent;
      
      // insert initial state
      EventRecord * evrec = new EventRecord();
      Interaction * interaction = new Interaction;
      evrec->AttachSummary(interaction);
            
      evrec->AddParticle(gOptInpHadPdgCode,kIStInitialState,-1,-1,-1,-1,p4h,  x4null);
      evrec->AddParticle(gOptTgtPdgCode,   kIStInitialState,-1,-1,-1,-1,p4tgt,x4null);
      
      // generate h+A eventw
      intranuke->ProcessEventRecord(evrec);
  
      // print generated event    
      if (ievent<50) LOG("testIntranuke", pINFO) << *evrec;
  
      // add event at the output ntuple
      ntpw.AddEventRecord(ievent, evrec);

	// analyze
      INukeFateHA_t fate = utils::inuke::FindhAFate(evrec);
      if (ievent<50) LOG("testIntranuke", pINFO) << "fate = " << INukeHadroFates::AsString(fate);

      // We don't want the specific fate data, just the main (8) fate types
      // if fate num is greater than 99, then use the first digit + 1
      // -- otherwise, fate is either 0 or 1, so just use that fate
      countfate[(((int) fate)>99)?((100 + (int) fate)/100):((int) fate)]++;

      // refresh the mc job monitor
      mcjmonitor.Update(ievent,evrec);
      
      ievent++;
      delete evrec;
      
      } // end loop events
  
  end=clock();  // ending time

  //-- save the generated MC events
  ntpw.Save();

  // output section
  LOG("testIntranuke", pINFO) << " Probe = " << probe_name << " (KE=" << (Eh-mh) << ")" << ENDL;
  LOG("testIntranuke", pINFO) << " Target = " << tgt_name << " (Z=" << Z << ") nuclear radius, area are: " << nuclear_radius << " fm, " << area << " fm**2 " << ENDL;
 
  double fm2tomb = 10.; //(units::fm2 / units::mb);
 
  int cnttot = 0;
  int nullint = countfate[1]; // no interactions
  double sigtot = 0, sigtoterr = 0;
  double sigtotScat = 0;
  double sigtotAbs  = 0;

  if(gOptDoOutput) gOptXsecFile << gOptInpHadKE;

  for(int k=0; k<nfates; k++) {
    if(k!=1) {
      cnttot += countfate[k];
      double ratio = countfate[k]/(double)ievent;
      sigma[k] = fm2tomb * area * ratio;
      sigma_err[k]   = fm2tomb * area * TMath::Sqrt(ratio*(1-ratio)/(double)ievent);
      if(sigma_err[k]==0) sigma_err[k] = fm2tomb * area * TMath::Sqrt(countfate[k])/(double)ievent;
    
      if(countfate[k]>0) {
        LOG("testIntranuke", pINFO) << " --> " << setw(26) << fatestr[k] << ": " << setw(7) << countfate[k] << " events -> " << setw(7) << sigma[k] << " +- " << sigma_err[k] << " (mb)" << ENDL;
      }
    
      if(k==5) sigtotAbs += sigma[k];
      else if (k!=1) sigtotScat += sigma[k];

      if (gOptDoOutput) gOptXsecFile << "    " << sigma[k] << "    " << sigma_err[k];
    }
  }

  sigtot    = fm2tomb * area * cnttot/(double)ievent;
  sigtoterr = fm2tomb * area * TMath::Sqrt(cnttot)/(double)ievent;

  double sigtot_noElas    = fm2tomb * area * (cnttot-countfate[3])/(double)ievent;
  double sigtoterr_noElas = fm2tomb * area * TMath::Sqrt(cnttot-countfate[3])/(double)ievent;

  if(gOptDoOutput) 
  {
    gOptXsecFile << "    " << sigtot_noElas << "    " << sigtoterr_noElas;
    gOptXsecFile << "    " << sigtot << "    " << sigtoterr << endl;
  }

  double ratioAS = (sigtotScat==0) ? 0 : sigtotAbs/(double)sigtotScat;
  
  LOG("testIntranuke", pINFO) << " -------------------------------------------------------------------------- " << ENDL;
  LOG("testIntranuke", pINFO) << " ==> " << setw(28) << " Runtime: " << setw(7) 
 	<< (((double)(end-start))/((double) CLOCKS_PER_SEC)) << ENDL;
  LOG("testIntranuke", pINFO) << " ==> " << setw(28) << " Total: " << setw(7) << cnttot << " events -> " << setw(7) << sigtot << " +- " << sigtoterr << " (mb)" << ENDL;
  LOG("testIntranuke", pINFO) << " (-> " << setw(28) << " Hadron escaped nucleus: " << setw(7) << nullint << " events) " << ENDL;
  LOG("testIntranuke", pINFO) << " ==> " << setw(28) << " Ratio (abs/scat) = " << setw(7) << ratioAS << ENDL;
  LOG("testIntranuke", pINFO) << " ==> " << setw(28) << " avg. num of int. = " << setw(7) << cnttot/(double)ievent << ENDL;
  LOG("testIntranuke", pINFO) << " ==> " << setw(28) << " no interaction   = " << setw(7) << (ievent-cnttot)/(double)ievent << ENDL; 
  LOG("testIntranuke", pINFO) << " -------------------------------------------------------------------------- " << ENDL;
  //  LOG("testIntranuke", pINFO) << " ===> " << setw(2) << " Nevt w/o hadrons scatter= " << setw(7) << nscatter[0] << " events " << ENDL;
  //  LOG("testIntranuke", pINFO) << " ===> " << setw(2) << " Nevt w/ a hadron scatter= " << setw(7) << nscatter[1] << " events " << ENDL;
  //  LOG("testIntranuke", pINFO) << " ===> " << setw(2) << " Nevt w/ mult-had scatter= " << setw(7) << nscatter[2] << " events " << ENDL;
  //  LOG("testIntranuke", pINFO) << " -------------------------------------------------------------------------- " << ENDL;
  //  LOG("testIntranuke", pINFO) << " ===> " << setw(2) << " Nevt generated          = " << setw(7) << nscatter[0]+nscatter[1]+nscatter[2] << " events " << ENDL;
  LOG("testIntranuke", pINFO) << ENDL;

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("testIntranuke", pNOTICE) << "Parsing command line arguments";

  //number of events:
  try {
    LOG("testIntranuke", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }

  //run number:
  try {
    LOG("testIntranuke", pINFO) << "Reading MC run number";
    gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }

  // incoming hadron PDG code:
  try {
    LOG("testIntranuke", pINFO) << "Reading rescattering particle PDG code";
    gOptInpHadPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
      if(!e.ArgumentFound()) {
      LOG("testIntranuke", pFATAL) << "Unspecified PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //target PDG code:
  try {
    LOG("testIntranuke", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO) << "Unspecified target - Using default";
      gOptTgtPdgCode = kDefOptTgtPdgCode;
    }
  }

  // incoming hadron kinetic energy:
  try {
    LOG("testIntranuke", pINFO) << "Reading rescattering particle KE energy";
    string ke = genie::utils::clap::CmdLineArgAsString(argc,argv,'k');
    gOptInpHadKE = atof(ke.c_str());
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pFATAL) << "Unspecified KE - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // xsec output filename:
  try {
    LOG("testIntranuke", pINFO) << "Reading xsec filename";
    gOptXsecFName = genie::utils::clap::CmdLineArgAsString(argc,argv,'o');
    gOptXsecFile.open(gOptXsecFName.c_str(), ofstream::app);
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptDoOutput = false;
    }
  }
  catch(ofstream::failure e) {
    LOG("testIntranuke", pFATAL) << "Error opening output file: ";
    LOG("testIntranuke", pFATAL) << "--> " << e.what();
  }

  // mode
  try {
    LOG("testIntranuke", pINFO) << "Reading Intranuke mode";
    gOptMode = genie::utils::clap::CmdLineArgAsString(argc,argv,'m');
    if(gOptMode.compare("hA")!=0 && gOptMode.compare("hN")!=0) {
      LOG("testIntranuke", pINFO) << "Invalid specified mode - Using default";
      gOptMode = kDefOptMode;
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO) << "Unspecified mode - Using default";
      gOptMode = kDefOptMode;
    }
  }

  LOG("testIntranuke", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("testIntranuke", pINFO) << "MC Run Number              = " << gOptRunNu;
  LOG("testIntranuke", pINFO) << "Incoming hadron PDG code   = " << gOptInpHadPdgCode;
  LOG("testIntranuke", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  LOG("testIntranuke", pINFO) << "Hadron input KE            = " << gOptInpHadKE;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("testIntranuke", pNOTICE)
    << "\n\n" 
    << "Syntax:" << "\n"
    << "   gtestIntranuke [-n nev] -p hadron_pdg -t tgt_pdg [-r run] [-a R0] -k KE -o output_file [-m mode]"
    << "\n";
}
//____________________________________________________________________________
