//____________________________________________________________________________
/*!

\program gnc2eN

\brief   Implements S.Wood's idea to obtain eN scattering events by reweighting
         vN NC events.

         Syntax :
           gnc2eN -f filename [-n events] 

         Options:
           [] Denotes an optional argument
           -f Specifies the GENIE/ROOT file with the generated event sample
           -n Specifies how many events to analyze [default: all]

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 14, 2007

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "Conventions/Constants.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpWriter.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::constants;

// function prototypes
void   GetCommandLineArgs (int argc, char ** argv);
void   PrintSyntax        (void);
double eNXSec             (const EventRecord & eN_event);

// command-line arguments
Long64_t gOptN;                // (-n)  process so many events, all if -1
string   gOptInpFile;          // (-f) input GENIE event sample file

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments 
  GetCommandLineArgs(argc,argv);

  LOG("gnc2eN", pNOTICE) << "*** Opening GHEP data file: " << gOptInpFile;

  TFile * inp_file = new TFile(gOptInpFile.c_str(),"READ");
  TTree * event_tree = dynamic_cast <TTree *> (inp_file->Get("gtree"));
  NtpMCTreeHeader * thdr = 
            dynamic_cast <NtpMCTreeHeader *> (inp_file->Get("header"));

  LOG("gnc2eN", pNOTICE) << "*** Input tree header: " << *thdr;

  //-- figure out the TTree format (GENIE supports multiple formats)
  NtpMCFormat_t format = thdr->format;
  assert(format == kNFGHEP);

  //-- set the event record branch ptr
  NtpMCEventRecord * mcrec = 0;
  event_tree->SetBranchAddress("gmcrec", &mcrec);

  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
    event_tree->GetEntries() : TMath::Min( event_tree->GetEntries(), gOptN );

  //-- initialize an Ntuple Writer
  int gOptRunNu = 0;
  NtpWriter ntpw(kNFGHEP, gOptRunNu);
  ntpw.Initialize("gntp-eN");
  
  int ieN=0;
  for(Long64_t ivN = 0; ivN < nmax; ivN++) {

    //-- access current event tree entry & get the event record
    event_tree->GetEntry(ivN);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  vN_event   = *(mcrec->event);

    //-- use only vN NC events
    bool is_nc = vN_event.Summary()->ProcInfo().IsWeakNC();
    if(!is_nc) continue;

    LOG("gnc2eN", pNOTICE) << "Input vN NC event: " << vN_event;

    //-- copy event - switch incoming/outgoing neutrinos to electrons
    EventRecord eN_event(vN_event);
    eN_event.Particle(0)->SetPdgCode(kPdgElectron);
    eN_event.Particle(eN_event.Particle(0)->FirstDaughter())->SetPdgCode(kPdgElectron);

    //-- reweight event
    double vN_xsec = vN_event.DiffXSec();
    double eN_xsec = eNXSec(eN_event);
    double wght    = eN_xsec/vN_xsec;
    eN_event.SetWeight(wght * eN_event.Weight());

    //-- save eN event to the new event tree
    LOG("gnc2eN", pNOTICE) << "Output eN event: "<< eN_event;
    ntpw.AddEventRecord(ieN++, &eN_event);
  }

  //-- save event tree
  ntpw.Save();

  LOG("gnc2eN", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
double eNXSec(const EventRecord & /*eN_event*/)
{
  return 1.;
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gnc2eN", pNOTICE) << "*** Parsing commad line arguments";

   //number of events:
  try {
    LOG("gnc2eN", pINFO) << "Reading number of events to analyze";
    gOptN = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gnc2eN", pINFO)
              << "Unspecified number of events to analyze - Use all";
      gOptN = -1;
    }
  }

  //get GENIE event sample ROOT file (ER-format)
  try {
    LOG("gnc2eN", pINFO) << "Reading event sample filename";
    gOptInpFile = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctst", pFATAL) << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gnc2eN", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gnc2eN -f file [-n nev]\n";
}
//_________________________________________________________________________________
