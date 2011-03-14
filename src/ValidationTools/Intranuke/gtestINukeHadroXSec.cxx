//____________________________________________________________________________
/*!

\program gtestINukeHadroXSec

\brief   INTRANUKE test program. Reads-in a hadron-nucleus event file in GHEP
         format and prints-out hadron-nucleus cross sections.

         Syntax :
           gtestINukeHadroXSec -f input_ghep_file [-w]

         -f : Input file
         -w : If set, writes computed hadron cross sections to a text file

\authors Aaron Meyer and Steve Dytman

\created July 26, 2010

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/INukeHadroFates.h"
#include "HadronTransport/INukeUtils.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/CmdLnArgParser.h"

using std::endl;
using std::setw;
using std::setfill;

using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

// command line options
string gOptInpFilename = "";    ///< input event file
bool   gOptWriteOutput = false; ///< write out hadron cross sections
string gOptOutputFilename = "gevgen_hadron_xsection.txt";

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // parse command line arguments
  GetCommandLineArgs(argc,argv);

  //
  // initialize
  //

  const int nfates = 9;             // total number of possible fates
  int           countfate [nfates]; // total no. of events with given fate
  double        sigma     [nfates]; // cross sections
  double        sigma_err [nfates]; // cross section errors
  string        fatestr   [nfates] = "  ";
  INukeFateHA_t fatetype  [nfates];

  fatetype[0] = kIHAFtUndefined;   
  fatetype[1] = kIHAFtNoInteraction;
  fatetype[2] = kIHAFtCEx;
  fatetype[3] = kIHAFtElas;
  fatetype[4] = kIHAFtInelas;
  fatetype[5] = kIHAFtAbs;
  fatetype[6] = kIHAFtKo;
  fatetype[7] = kIHAFtPiProd;
  fatetype[8] = kIHAFtDCEx;
 
  for (int k=0; k<nfates; k++) {
    countfate[k] = 0; 
    sigma    [k] = 0.; 
    sigma_err[k] = 0.;
    fatestr  [k] = INukeHadroFates::AsString(fatetype[k]);
  }

  // event sample info (to be extracted from 1st event)
  int    probe_pdg  = 0;
  int    target_pdg = 0;
  double kin_energy = 0.;

  //
  // open the input ROOT file and get the event tree
  //

  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) {
    LOG("gtestINukeHadroXSec", pERROR) 
         << "No event tree found in input file!";
    return 1;
  }

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (int) tree->GetEntries();
  LOG("gtestINukeHadroXSec", pNOTICE) 
     << "Processing " << nev << " events";

  //
  // event loop
  //

  for(int ievent = 0; ievent < nev; ievent++) {

    // get next tree entry
    tree->GetEntry(ievent);

    // get the corresponding GENIE event
    EventRecord & event = *(mcrec->event);

    // extract info for the event sample
    if(ievent==0) {
       kin_energy = event.Particle(0)->KinE();
       probe_pdg  = event.Particle(0)->Pdg();
       target_pdg = event.Particle(1)->Pdg();
    }

    // analyze
    const GHepRecord * grec = dynamic_cast<const GHepRecord *> (&event);
    INukeFateHA_t fate = utils::intranuke::FindhAFate(grec);
    if(ievent<100) {
       LOG("gtestINukeHadroXSec", pNOTICE) 
          << "fate = " << INukeHadroFates::AsString(fate);
    }

    // We don't want the specific fate data, just the main (9) fate types
    switch (fate){
      case 0:   countfate[0]++; break;
      case 1:   countfate[1]++; break;
      case 2:   countfate[2]++; break;
      case 3:   countfate[3]++; break;
      case 4:   countfate[4]++; break;
      case 5:   countfate[5]++; break;
      case 6:   countfate[6]++; break;
      case 13:  countfate[8]++; break;
      default:  
	if (7<=fate && fate<=12) countfate[7]++;
	else {
         LOG("gtestINukeHadroXSec", pWARN) 
            << "Undefined fate from FindhAFate() : " << fate;
        }
	break;
    }

    // clear current mc event record
    mcrec->Clear();

  } // end event loop 
  
  //
  // output section
  //
 
  const double fm2tomb  = units::fm2 / units::mb;
  const double dnev     = (double) nev;
  const int    NR       = 3;
  const double R0       = 1.4;

  int    A              = pdg::IonPdgCodeToA(target_pdg);
  int    Z              = pdg::IonPdgCodeToZ(target_pdg);
  double nuclear_radius = NR * R0 * TMath::Power(A, 1./3.); // fm
  double area           = TMath::Pi() * TMath::Power(nuclear_radius,2);

  PDGLibrary * pdglib = PDGLibrary::Instance();
  string probe_name  = pdglib->Find(probe_pdg)->GetName();
  string target_name = pdglib->Find(target_pdg)->GetName();

  LOG("gtestINukeHadroXSec", pINFO) 
     << " Probe = " << probe_name 
     << ", KE = " << kin_energy << " GeV";
  LOG("gtestINukeHadroXSec", pINFO) 
    << " Target = " << target_name 
    << " (Z,A) = (" << Z << ", " << A 
    << "), nuclear radius = " << nuclear_radius 
    << " fm, area = " << area << " fm**2 " << '\n';

  int    cnttot     = 0;
  int    nullint    = countfate[1]; // no interactions
  double sigtot     = 0;
  double sigtoterr  = 0;
  double sigtotScat = 0;
  double sigtotAbs  = 0;

  for(int k=0; k<nfates; k++) {
    if(k!=1) {
      cnttot += countfate[k];
      double ratio = countfate[k]/dnev;
      sigma[k]     = fm2tomb * area * ratio;
      sigma_err[k] = fm2tomb * area * TMath::Sqrt(ratio*(1-ratio)/dnev);
      if(sigma_err[k]==0) {
         sigma_err[k] = fm2tomb * area * TMath::Sqrt(countfate[k])/dnev;
      }
      if(countfate[k]>0) {
        LOG("gtestINukeHadroXSec", pINFO) 
            << " --> " << setw(26) << fatestr[k] 
            << ": " << setw(7) << countfate[k] << " events -> " 
            << setw(7) << sigma[k] << " +- " << sigma_err[k] << " (mb)" << '\n';
      }
      if(k==5) {
        sigtotAbs += sigma[k];
      }
      else 
      if (k!=1) { 
        sigtotScat += sigma[k];
      }
    }//k!=1
  }//k koop

  sigtot    = fm2tomb * area * cnttot/dnev;
  sigtoterr = fm2tomb * area * TMath::Sqrt(cnttot)/dnev;

  double sigtot_noelas    = fm2tomb * area * (cnttot-countfate[3])/dnev;
  double sigtoterr_noelas = fm2tomb * area * TMath::Sqrt(cnttot-countfate[3])/dnev;

  double ratio_as = (sigtotScat==0) ? 0 : sigtotAbs/(double)sigtotScat;
  
  LOG("gtestINukeHadroXSec", pNOTICE) 
    << "\n\n --------------------------------------------------- " 
    << "\n ==> " << setw(28) << " Total: " << setw(7) << cnttot 
    << " events -> " << setw(7) << sigtot << " +- " << sigtoterr << " (mb)"
    << "\n (-> " << setw(28) << " Hadrons escaped nucleus: " 
    << setw(7) << nullint << " events)"
    << "\n ==> " << setw(28) << " Ratio (abs/scat) = " 
    << setw(7) << ratio_as
    << "\n ==> " << setw(28) << " avg. num of int. = " 
    << setw(7) << cnttot/dnev
    << "\n ==> " << setw(28) << " no interaction   = " 
    << setw(7) << (dnev-cnttot)/dnev
    << "\n ------------------------------------------------------- \n";

  if(gOptWriteOutput) 
  {
    ofstream xsec_file; 
    xsec_file.open(gOptOutputFilename.c_str(), std::ios::app);
    xsec_file << kin_energy;
    for(int k=0; k<nfates; k++) {
       xsec_file << "\t" << sigma[k] << "\t" << sigma_err[k];
    }
    xsec_file << "\t" << sigtot_noelas << "\t" << sigtoterr_noelas;
    xsec_file << "\t" << sigtot        << "\t" << sigtoterr << endl;
    xsec_file.close();
  }

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gtestINukeHadroXSec", pNOTICE) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get input ROOT file (containing a GENIE GHEP event tree)
  if( parser.OptionExists('f') ) {
    LOG("gtestINukeHadroXSec", pINFO) << "Reading input filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("gtestINukeHadroXSec", pFATAL)
       << "Unspecified input filename - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // write-out events?
  gOptWriteOutput =  parser.OptionExists('w');
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gtestINukeHadroXSec", pNOTICE)
    << "\n\n" 
    << "Syntax:" << "\n"
    << "  gtestINukeHadroXSec -f event_file [-w]"
    << "\n";
}
//____________________________________________________________________________
