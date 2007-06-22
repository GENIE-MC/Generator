//____________________________________________________________________________
/*!

\program testIntranuke

\brief   testIntranuke from hadtest.F

         Syntax :
           testIntranuke [-n nev] [-t tgtpdg] [-f format] [-r run#] [-a R0] -i inputpdg -k KE

         Options :
           [] denotes an optional argument
           -n specifies the number of events to generate
           -t specifies the target PDG code (std format: 1aaazzz000)
           -f specifies the output TTree format. If set to 0 will create a
              single-branch TTree with NtpMCPlainRecord objects in its leaves,
              while if set to 1 it will have NtpMCEventRecord objects in
              its leaves (see the Ntuple package for descriptions of the ntuple
              records and their intended usage). Default options is 1.
           -r specifies the MC run number
           -i specifies the rescattering particle PDG code as a input
           -k specifies the rescattering particle kinetic energy
              (if what follows the -k option is a comma separated pair of
               values it will be interpreted as an energy range and input
               energies will be generated uniformly in this range)
	      *Note that this is the kinetic energy of an "scattering" particle, not
	      the energy of an "interacting" or a "flux" neutrino -- the particle is forced to
	      scatter anyway (see gEvGen for how to input neutrinos)
	   -l length of step corresponding to some amount of KE (step size)
	   -m number of steps with some amount of KE
           -a specifies an effective nucleus size given in fm (R=Ro*A^1/3)

         Example:
           gtestIntranuke -n 300 -t 1056026000 -a 1.2 -i 211 -k 0.165

           will generate 300 events hadron-nucleus interaction and show cross section.
	   A nucleus is Iron (A=56,Z=26) as default target.
	   A hadron is pi+ (pdg = 211) with R0 = 1.2 as defalt option, KE = .165 GeV, 
	   and a vertex as a point on the outer edge using an effective nucleus size.

\author  Minsuk Kim and Steve Dytman
         University of Pittsburgh

\version 1.1

\created May 1, 2007

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TNtuple.h>

#include <iomanip>

#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Conventions/Constants.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodeList.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepStatus.h"
#include "HadronTransport/Intranuke.h"
#include "HadronTransport/INukeHadroData.h"
#include "HadronTransport/INukeHadroFates.h"
#include "EVGCore/EventRecordVisitorI.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

using std::cout;
using std::endl;
using std::ios;
using std::setw;

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 100;            // n-events to generate
Long_t        kDefOptRunNu     = 0;              // default run number

int           kDefOptTgtPdgCode   = 1056026000;
double        kDefOptR0           = 1.2;
int           kDefOptNstep        = 1; 
double        kDefOptLstep        = 0.02;        // some amount of KE (GeV)

//User-specified options:
int           gOptNevents;           // n-events to generate
int           gOptTgtPdgCode;        // target PDG code
NtpMCFormat_t gOptNtpFormat;         // ntuple format
Long_t        gOptRunNu;             // run number

int           gOptInputPdgCode;      // rescattering particle PDG code as a input
double        gOptInputKE;           // This is KE = E - M. So E = M + KE
double        gOptRangeKE;           // max-min KE
double        gOptR0;                // R0
double        gOptLstep;             // Length of step corresponding to some amount of KE (size) 
int           gOptNstep;             // Number of steps with some amount of KE

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- print the options you got from command line arguments

  string fmts = NtpMCFormat::AsString(gOptNtpFormat);
  if(gOptNstep>1) gOptRangeKE = gOptLstep*(gOptNstep-1);

  LOG("testIntranuke", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("testIntranuke", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  LOG("testIntranuke", pINFO) << "MC Run Number              = " << gOptRunNu;
  LOG("testIntranuke", pINFO) << "Output ntuple format       = " << fmts;
  LOG("testIntranuke", pINFO) << "Effective nucleus size, R0 = " << gOptR0;
  LOG("testIntranuke", pINFO) << "Length of step             = " << gOptLstep;
  LOG("testIntranuke", pINFO) << "Number of steps            = " << gOptNstep;
  if(gOptRangeKE>0) {
    LOG("testIntranuke", pINFO) << "Hadron input KE            = ["
        << gOptInputKE << ", " << gOptInputKE+gOptRangeKE << "]";
  } else {
    LOG("testIntranuke", pINFO) << "Hadron input KE            = " << gOptInputKE;
  }

  //-- Get an instance of PDG library to access PDG particle data (masses etc)
  PDGLibrary * pdglib = PDGLibrary::Instance();

  //-- Get a handle to the algorithm factory and ask it to give you an algorithm
  //-- called 'genie::Intranuke' configured with a parameter set named 'Default'
  //-- and implementing an algorithm interface named 'EventRecordVisitorI'.
  AlgFactory * algf = AlgFactory::Instance();

  const EventRecordVisitorI * intranuke = 
                 dynamic_cast<const EventRecordVisitorI *> (
                         algf->GetAlgorithm("genie::Intranuke","Default"));

  //-- Access the algorithm's configuration
  AlgConfigPool * confpool = AlgConfigPool::Instance();
  Registry * config = confpool->FindRegistry(intranuke);

  //-- Set the 'test-mode' parameter to true (same as fortran intranuke's pitest mode)
  config->UnLock();
  bool test=true;
  config->Set("test-mode", test);
  if(gOptR0!=1.2) config->Set("R0",gOptR0);
  ////config->Set("Kpt2",1);
  ////config->Set("ct0","2");
  //config->Set("mode","hN");
  //config->Set("nuc-removal-energy",0.007); //config->Set("nucleon-removal-energy",0.007); 

  //-- Read modified configurations
  algf->ForceReconfiguration();

  int NNUCLN = pdg::IonPdgCodeToA(gOptTgtPdgCode); // from $GENIE/PDG/PDGUtils.h
  //double nuclear_density = GV_RHONUC;
  ////gOptR0 = gOptR0 + 0.00000005; (in neugen)
  double nuclear_radius = gOptR0*TMath::Power(NNUCLN,1./3.); //exactly rsiz*TMath::Power(NNUCLN,1/3.);

  double area = TMath::Pi()*TMath::Power(nuclear_radius,2);
  double mass = pdglib->Find(gOptInputPdgCode)->Mass();
  double nuclear_vertex[3];
  for(int k=0; k<3; k++) nuclear_vertex[k] = 0;
  const int nfates = 12;
  const int nsteps = gOptNstep; double stepKE[nsteps];
  int countfate[nsteps][nfates]; string FateType[nsteps][nfates];
  double sigma[nfates], err[nfates];
  char fname[80]; sprintf(fname,"xsec-%d.root",static_cast<Int_t>(gOptRunNu));
  TFile f(fname,"RECREATE");
  TNtuple* nt = new TNtuple("nt","nt","probe:ke:step:sigtot:sigcex:sigelas:siginelas:sigabs:sigprod:sigtoterr:sigcexerr:sigelaserr:siginelaserr:sigabserr:sigproderr");
  TNtuple* nt2 = new TNtuple("nt2","nt2","probe:ke:step:energy:fate:fspid:fske");
  int bin = 100; if(gOptNstep>1) bin = gOptNstep;
  TH1F *he = new TH1F("he","he",bin,mass+gOptInputKE,mass+gOptInputKE+gOptRangeKE);
  for(int istep=0; istep<gOptNstep; istep++) {
    stepKE[istep] = 0;
    for(int k=0; k<nfates; k++) {
      countfate[istep][k] = 0; if(istep==0) sigma[k] = err[k] = 0;
    }
  }

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(gOptNtpFormat, gOptRunNu);
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- get a random number generator
  RandomGen * r = RandomGen::Instance();
  time_t t; (void) time(&t);
  r->SetSeed(t);

  for(int istep=0; istep<gOptNstep; istep++) {

    if(istep>0) gOptInputKE += gOptLstep; stepKE[istep] = gOptInputKE;

    //-- generate events / print the GHEP record / add it to the ntuple
    int ievent = 0;
    while ( ievent < gOptNevents) {

      if(ievent<100) cout << " *** Generating event............ " << ievent << endl;
      
      EventRecord * evrec = new EventRecord();
      
      double energy = mass + gOptInputKE;
      //-- generate input energy (if an energy range was defined)
      if(gOptNstep==1 && gOptRangeKE>0) energy = mass + gOptInputKE + gOptRangeKE * r->RndEvg().Rndm();
      he->Fill(energy);
      double pz = TMath::Sqrt(energy*energy - mass*mass);
      
      TLorentzVector p4new(0.,0.,pz,energy);
      TLorentzVector p4tgt(0,0,0,pdglib->Find(gOptTgtPdgCode)->Mass());
      TLorentzVector x4null(0,0,0,0);
      
      evrec->AddParticle(gOptInputPdgCode,kIStInitialState,-1,-1,-1,-1,p4new,x4null);
      evrec->AddParticle(gOptTgtPdgCode  ,kIStInitialState,-1,-1,-1,-1,p4tgt,x4null);
      
      //-- Print the event record that will be input to intranuke
      //LOG("testIntranuke", pDEBUG) << "Intranuke input event: " << *evrec;
      //-- Let intranuke do its job
      intranuke->ProcessEventRecord(evrec);
      
      bool interacts = false;
      INukeFateHA_t fate = kIHAFtUndefined;
      int nn = 0, np = 0, npip = 0, npi0 = 0, npim = 0;
      TObjArrayIter piter(evrec);
      GHepParticle * p = 0;
      int icurr = -1;
      while( (p = (GHepParticle *) piter.Next()) ) {
	icurr++;
	//if(icurr<2) cout << icurr << " " << print::X4AsString(p->X4()) << endl;
	if(icurr==2) {
	  LOG("testIntranuke", pDEBUG) << " -> Event No. " << setw(4) << ievent << setw(12) << p->Pdg() << setw(5) << p->Name() << setw(3) << p->Status();
	  LOG("testIntranuke", pDEBUG) << " -> KinE = " << p->KinE() << " " << print::P4AsShortString(p->P4());
	  LOG("testIntranuke", pDEBUG) << " -> rs = " << p->X4()->Vect().Mag() << " " << print::X4AsString(p->X4()) << ", NuclearRadius=" << nuclear_radius;
	  if(p->X4()->Vect().Mag() < nuclear_radius) interacts = true;
	  if(p->Status() == kIStStableFinalState) {
	    if(interacts) {
	      if(p->Pdg() != gOptInputPdgCode) fate = kIHAFtCEx;
	      else fate = (TMath::Abs(p->P4()->Vect().Mag()-pz)<1E-5) ? kIHAFtElas : kIHAFtInelas;
	    }
	    nt2->Fill(gOptInputPdgCode,gOptInputKE,istep,energy,fate,p->Pdg(),p->KinE());
	    break;
	  }
	} else if(p->Status() == kIStStableFinalState) {
	  if(p->Pdg() == kPdgProton) np++;
	  if(p->Pdg() == kPdgNeutron) nn++;
	  if(p->Pdg() == kPdgPiP) npip++;
	  if(p->Pdg() == kPdgPiM) npim++;
	  if(p->Pdg() == kPdgPi0) npi0++;
	}
      }
      
      if(interacts) {

	int npi = npip + npi0 + npim;
	if(npi==0) { 
	  if(nn==1 && np==1) fate = kIHAFtAbsNP;
	  if(np==2 && nn==0) fate = kIHAFtAbsPP;
	  if(nn==1 && np==2) fate = kIHAFtAbsNPP;
	  if(nn==2 && np==1) fate = kIHAFtAbsNNP;
	  if(nn==2 && np==2) fate = kIHAFtAbs2N2P; // only for rescattering pion
	  if(nn==2 && np==3) fate = kIHAFtAbs2N3P; // either rescattering neutron or proton
	} else {
	  if(nn==1 && npip==1 && npi0==0) fate = kIHAFtNPip; // either rescattering neutron or proton
	  if(nn==1 && npip==1 && npi0==1) fate = kIHAFtNPipPi0;
	}

      }

      countfate[istep][fate]++;
      FateType[istep][fate] = INukeHadroFates::AsString(fate).c_str();
    
      if(ievent<100) cout << " =======> Selected fate: " << INukeHadroFates::AsString(fate) << endl;
      if(ievent<100) cout << " =======> Generated Event GHEP Record: " << *evrec << endl;
      
      //-- add event at the output ntuple
      ntpw.AddEventRecord(ievent, evrec);
      
      //-- refresh the mc job monitor
      mcjmonitor.Update(ievent,evrec);
      
      ievent++;
      delete evrec;
      
    } // end loop events

  } // end loop for KE step
  
  //-- save the generated MC events
  ntpw.Save();

  cout << endl << endl;
  cout << " A = " << NNUCLN << " nuclear radius, area are: " << nuclear_radius << " fm, " << area << " fm**2 " << endl;
  cout << " Total cross section results " << endl;
  for(int istep=0; istep<nsteps; istep++) {

    cout << " " << endl;
    cout << "  For " << pdglib->Find(gOptInputPdgCode)->GetName() << " (KE=" << stepKE[istep] << ") + " << pdglib->Find(gOptTgtPdgCode)->GetName() << " interaction " << endl;
    cout << "  -------------------------------------------------------------------------- " << endl;

    double fm2tomb = (units::fm2 / units::mb);

    int cnttot = 0;
    //int nullint = countfate[istep][0];
    double sigtot = 0, sigtoterr = 0;
    double sigtotScat = 0, sigtotAbs = 0, sigtotProd = 0;

    for(int k=1; k<nfates; k++) {
      
      cnttot += countfate[istep][k];
      double ratio = countfate[istep][k]/(double)gOptNevents;
      sigma[k] = fm2tomb * area * ratio;
      err[k]   = fm2tomb * area * TMath::Sqrt(ratio*(1-ratio)/(double)gOptNevents);
      if(err[k]==0) err[k] = fm2tomb * area * TMath::Sqrt(countfate[istep][k])/(double)gOptNevents;

      if(countfate[istep][k]>0) {
	cout << "  --> " << setw(26) << FateType[istep][k] << ": " << setw(7) << countfate[istep][k] << " events -> " << setw(7) << sigma[k] << " +- " << err[k] << " (mb)" << endl;
      }

      if(k>=1 && k<=3) sigtotScat += sigma[k];
      if(k>=4 && k<=9) sigtotAbs += sigma[k];
      if(k>=10 && k<=11) sigtotProd += sigma[k];
      
    }
    sigtot    = fm2tomb * area * cnttot/(double)gOptNevents;
    sigtoterr = fm2tomb * area * TMath::Sqrt(cnttot)/(double)gOptNevents;
    double ratioAS = (sigtotScat==0) ? 0 : sigtotAbs/(double)sigtotScat;
    
    cout << "  -------------------------------------------------------------------------- " << endl;
    cout << "  ==> " << setw(28) << " Total: " << setw(7) << cnttot << " events -> " << setw(7) << sigtot << " +- " << sigtoterr << " (mb)" << endl;
    cout << "  ==> " << setw(28) << " Ratio (abs/scat) = " << setw(7) << ratioAS << endl;
    cout << "  ==> " << setw(28) << " avg. num of int. = " << setw(7) << cnttot/(double)gOptNevents << endl;
    cout << "  ==> " << setw(28) << " no interaction   = " << setw(7) << (gOptNevents-cnttot)/(double)gOptNevents << endl; 
    cout << endl;

    double sigcex = sigma[1], sigelas = sigma[2], siginelas = sigma[3], sigabs = sigma[4]+sigma[5]+sigma[6]+sigma[7]+sigma[8]+sigma[9], sigprod = sigma[10]+sigma[11];
    double sigcexerr = err[1], sigelaserr = err[2], siginelaserr = err[3], sigabserr = err[4]+err[5]+err[6]+err[7]+err[8]+err[9], sigproderr = err[10]+err[11];
    nt->Fill(gOptInputPdgCode,stepKE[istep],istep,sigtot,sigcex,sigelas,siginelas,sigabs,sigprod,sigtoterr,sigcexerr,sigelaserr,siginelaserr,sigabserr,sigproderr);

  }
  nt->Write(); nt2->Write(); he->Write();
  f.Write();
  f.Close();

  //-- Print the algorithm id and configuration parameters
  LOG("testIntranuke", pNOTICE) << *intranuke;

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("testIntranuke", pNOTICE) << "Parsing command line arguments";

  //-- Optional arguments

  //rescattering particle PDG code:
  try {
    LOG("testIntranuke", pINFO) << "Reading rescattering particle PDG code";
    gOptInputPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'i');
    if(gOptInputPdgCode != 111 && TMath::Abs(gOptInputPdgCode) != 211 && gOptInputPdgCode != 2112 && gOptInputPdgCode != 2212) {
      LOG("testIntranuke", pFATAL) << "Incorrect PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  } catch(exceptions::CmdLineArgParserException e) {
      if(!e.ArgumentFound()) {
      LOG("testIntranuke", pFATAL) << "Unspecified PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //rescattering particle kinetic energy
  try {
    LOG("testIntranuke", pINFO) << "Reading rescattering particle KE energy";
    ////gOptInputKE = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'k');
    string ke = genie::utils::clap::CmdLineArgAsString(argc,argv,'k');

    // is it just a value or a range (comma separated set of values)
    if(ke.find(",") != string::npos) {
       // split the comma separated list
       vector<string> kerange = utils::str::Split(ke, ",");
       assert(kerange.size() == 2);
       double emin = atof(kerange[0].c_str());
       double emax = atof(kerange[1].c_str());
       assert(emax>emin && emin>0);
       gOptInputKE = emin;
       gOptRangeKE = emax-emin;
    } else {
       gOptInputKE = atof(ke.c_str());
       gOptRangeKE = -1;
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pFATAL) << "Unspecified KE - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //length of step:
  try {
    LOG("testIntranuke", pINFO) << "Reading length of step corresponding to some amount of KE";
    gOptLstep = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'l');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO) << "Unspecified length of step - Using default";
      gOptLstep = kDefOptLstep;
    }
  }

  //number of steps:
  try {
    LOG("testIntranuke", pINFO) << "Reading number of steps with some amount of KE";
    gOptNstep = genie::utils::clap::CmdLineArgAsInt(argc,argv,'m');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("testIntranuke", pINFO) << "Unspecified number of steps - Using default";
      gOptNstep = kDefOptNstep;
    }
  }

  try {
    LOG("testIntranuke", pINFO) << "Reading RO";
    gOptR0 = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'a');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptR0 = kDefOptR0;
    }
  }

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

  //output ntuple format
  int format = 1; format = 0;
  try {
    LOG("gevgen", pINFO) << "Reading requested output ntuple format";
    format = genie::utils::clap::CmdLineArgAsInt(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pINFO) << "Unspecified tree format - Using default";
    }
  }
  if(format == 0 || format == 1) gOptNtpFormat = (NtpMCFormat_t)format;

  //target PDG code:
  try {
    LOG("testIntranuke", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      gOptTgtPdgCode = kDefOptTgtPdgCode;
      //LOG("testIntranuke", pFATAL) << "Unspecified target PDG code - Exiting";
      //PrintSyntax();
      //exit(1);
    }
  }

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("testIntranuke", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gtestIntranuke [-n nev] [-t tgtpdg] [-f format] [-r run] [-a R0] -i inputpdg -k KE\n"
    << "    for inputpdg : piplus  =  211\n"
    << "                   pizero  =  111\n"
    << "                   piminus = -211\n"
    << "                   neutron = 2112\n"
    << "                   proton  = 2212\n\n"
    << "    to generate 1k events with intranuclear rescattering in the default target (pi+,Fe56)\n"
    << "    gtestIntranuke -n 1000 -i 211 -k .165\n\n"
    << "    to generate 1k events with intranuclear rescattering in the carbon target (pi+,C12)\n"
    << "    gtestIntranuke -n 1000 -i 211 -k .165 -t 1016008000 -a 1.4\n\n";
}
//____________________________________________________________________________
