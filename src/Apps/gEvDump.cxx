//____________________________________________________________________________
/*!

\program gevdump

\brief   A GENIE utility printing-out GHEP event trees.

         *** Synopsis :

         gevdump -f filename 
                [-n n1[,n2]] 
                [--event-record-print-level]

         [] denotes an optional argument

         -f 
            Specifies a GENIE GHEP/ROOT event file.
         -n 
            Specifies range of events to print-out (default: all)
         --event-record-print-level
            Allows users to set the level of information shown when the event
            record is printed in the screen. See GHepRecord::Print().

         Examples:

         1. Print out all events from /data/sample.ghep.root 
            % gevdump -f /data/sample.ghep.root

         2. Print out the first 500 events from /data/sample.ghep.root 
            % gevdump -f /data/sample.ghep.root -n 0,499

         3. Print out the event 178 from /data/sample.ghep.root 
            shell$ gevdump -f /data/sample.ghep.root -n 178

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created September 02, 2005

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#include "Tools/Flux/GJPARCNuFlux.h"
#include "Tools/Flux/GNuMIFlux.h"
#endif 

using std::string;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
void GetEventRange      (Long64_t nev, Long64_t & n1, Long64_t & n2);

Long64_t gOptNEvtL;
Long64_t gOptNEvtH;
string   gOptInpFilename;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // set print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
  
  // Add a dummy value to Dark Matter to allow reading DarkMatter files
  PDGLibrary::Instance()->AddDarkMatter( 1.0, 0.5 );
  
  //
  // open the ROOT file and get the TTree & its header
  //

  TFile file(gOptInpFilename.c_str(),"READ");

  TTree * ghep_tree = 
     dynamic_cast <TTree *> (file.Get("gtree"));
  if(!ghep_tree) {
    LOG("gevdump", pFATAL) 
        << "No GHEP event tree in input file: " << gOptInpFilename;
    gAbortingInErr=true;
    exit(1);
  }
  Long64_t nev = ghep_tree->GetEntries();
  LOG("gevdump", pFATAL) 
     << "Input GHEP event tree has " << nev 
     << ((nev==1) ? " entry." : " entries.");

  NtpMCTreeHeader * thdr = 
     dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
  LOG("gevdump", pNOTICE) 
     << "Input tree header: " << *thdr;

  //  
  // set branch addresses
  //

  // main event record branch, always present
  NtpMCEventRecord * mcrec = 0;
  ghep_tree->SetBranchAddress("gmcrec", &mcrec);

  // if the event file was created by GENIE's gevpick `cherry-picking' app 
  // (see $GENIE/src/stdapp/gEvPick.cxx) then there will be additional branches
  // holding the original event filename and event number (in that file) 
  // for each `cherry-picked' event.
  bool have_gevpick_branches = false;
  TObjString* orig_filename = 0;
  Long64_t    orig_evtnum;
  TBranch * brOrigFilename = ghep_tree->GetBranch("orig_filename");
  TBranch * brOrigEvtNum   = ghep_tree->GetBranch("orig_evtnum");
  if(brOrigFilename!=0 && brOrigEvtNum!=0) {
    have_gevpick_branches = true;
    brOrigFilename->SetAddress(&orig_filename);
    brOrigEvtNum  ->SetAddress(&orig_evtnum);
  }

  // if the event file was created by one of GENIE's specialized event generation 
  // then there may be additional branches holding flux pass-through 
  // info (flux neutrino parent info for each generated event).
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
  flux::GJPARCNuFluxPassThroughInfo * jparc_flux_info = 0;
  flux::GNuMIFluxPassThroughInfo *    gnumi_flux_info = 0;
  TBranch * brFluxInfo = ghep_tree->GetBranch("flux");
  if(brFluxInfo) {
    TObjArray * leafarr = brFluxInfo->GetListOfLeaves();
    TIter iter(leafarr);
    TLeaf * leaf = 0;
    while((leaf = (TLeaf*)iter.Next())) {
      string ltypename = leaf->GetTypeName();
      if(ltypename == "genie::flux::GJPARCNuFluxPassThroughInfo") {
        brFluxInfo->SetAddress(&jparc_flux_info);
        LOG("gevdump", pNOTICE) 
          << "Found JPARC neutrino flux pass-though info";
      }
      else
      if(ltypename == "genie::flux::GNuMIFluxPassThroughInfo") {
        brFluxInfo->SetAddress(&gnumi_flux_info);
        LOG("gevdump", pNOTICE) 
          << "Found NuMI neutrino flux pass-though info";
      }
    }//leaf
  }//flux branch
#endif


  //
  // event loop
  //

  Long64_t n1,n2;
  GetEventRange(nev,n1,n2);
  for(Long64_t i = n1; i <= n2; i++) {
    ghep_tree->GetEntry(i);

    // retrieve GHEP event record abd print it out.
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event = *(mcrec->event);
    LOG("gevdump", pNOTICE) 
       << " ** Event: " << rec_header.ievent 
       << event;

    // print info from additional tree branches that might be present
    // if the event file was created by GENIE's gevpick app.
    if(have_gevpick_branches) {
     LOG("gevdump", pNOTICE) 
        << "\n Above event was originally event: " << orig_evtnum
        << "\n in event file: " << orig_filename->GetString().Data()
        << "\n\n";
    }

    // print info from additional JPARC or NuMI flux pass-through branches
    // that might be present of the event file was created by GENIE's
    // specialized event generation applications for T2K or NuMI-expts.
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
    if(jparc_flux_info) {
      LOG("gevdump", pNOTICE) 
        << "Associated JPARC flux pass-through info for above event:"
        << *jparc_flux_info;
    }
    if(gnumi_flux_info) {
      LOG("gevdump", pNOTICE) 
        << "Associated NuMI flux pass-through info for above event:"
        << *gnumi_flux_info;
    }
#endif

    mcrec->Clear();
  }

  // clean-up

  file.Close();

  LOG("gevdump", pNOTICE)  << "Done!";
  return 0;
}
//___________________________________________________________________
void GetEventRange(Long64_t nev, Long64_t & n1, Long64_t & n2)
{
  if(gOptNEvtL == -1 && gOptNEvtH == -1) {
    // read all events
    n1=0;
    n2=nev-1;
  }
  else {
    // read a range of events
    n1 = TMath::Max((Long64_t)0,  gOptNEvtL);
    n2 = TMath::Min(nev-1,        gOptNEvtH);
    if(n2-n1 <0) {
      LOG("gevdump", pFATAL) << "Invalid event range";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    } 
  }
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevdump", pINFO) << "*** Parsing command line arguments";

  // Common run options.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if ( parser.OptionExists('f') ) {
    LOG("gevdump", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("gevdump", pFATAL) 
       << "Unspecified input filename - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // number of events:
  if ( parser.OptionExists('n') ) {
    LOG("gevdump", pINFO) << "Reading number of events to analyze";
    string nev =  parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n',",");
      if(vecn.size()!=2) {
         LOG("gevdump", pFATAL) << "Invalid syntax";
         PrintSyntax();
         gAbortingInErr = true;
         exit(1);
      }
      // read a range of events
      gOptNEvtL = vecn[0];
      gOptNEvtH = vecn[1];       
    } else {
      // read single event
      gOptNEvtL = parser.ArgAsLong('n');
      gOptNEvtH = gOptNEvtL;
    }
  } else {
    LOG("gevdump", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvtL = -1;
    gOptNEvtH = -1;
  }

  
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevdump", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gevdump -f sample.root [-n n1[,n2]] [--event-record-print-level]\n";
}
//_________________________________________________________________________________
