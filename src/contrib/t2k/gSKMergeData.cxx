//____________________________________________________________________________
/*!

\program gSKmergedata

\brief   The program merges SuperK DETSIM/reconstruction info (for those events 
         passing SuperK nu_e or nu_mu selection cuts) with the full GHEP event 
         info used for simulating those events at the first place.

         This program serves as a workaround to errors in SuperK simulation
         chain that doesn't allow the full/correct GENIE info to be propagated
         to the analysis DSTs.

         Used by the RAL analysis for calculating SuperK systematics.

         Syntax :
           gSKmergedata -l selected_events_file -t genie_file_template

         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.

         Options :
           -l  A ROOT file with a TTree containing SKDETSIM and cut info
               for SK events in addition to pointers (run nu, event nu)
               to the original GENIE GHEP event file
           -t  A template for forming GENIE GHEP filenames and looking-up
               up the complete GENIE MC truth info.
               eg `-t /opt/t2k/data/sk/2010a/numu/genie_sk.2010a.numu.%d.ghep.root'

\author  Jim Dobson
         Imperial College London

         Costas Andreopoulos, Jelena Ilic, Nick Grant
         STFC, Rutherford Appleton Laboratory

\created June 01, 2010

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#include <cstdlib>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpWriter.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using namespace genie;

void          GetCommandLineArgs (int argc, char ** argv);
EventRecord * GetEvent           (int irunnu, int ievnu);

string gOptEventList;       ///< file with events passing cuts + SKDETSIM/reconstruction info
string gOptFilenmTemplate;  ///< template for forming filename of original GENIE/GHEP files fed-in to SKDETSIM

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // create an ntuple writer
  NtpWriter ntpw(kNFGHEP, 0);
  ntpw.CustomizeFilenamePrefix("allsk");
  ntpw.Initialize();

  // add new branches in the output tree to hold (in addition to the
  // GHEP event) SK reconstructed info
  double skEnuReco = 0;
  TBranch * brEnuReco = 
     ntpw.EventTree()->Branch("skEnuReco", &skEnuReco, "skEnuReco/D");
  brEnuReco->SetAutoDelete(kFALSE);
  int skNCutsPassedCumul = 0;
  TBranch * brNCutsPassedCumul = 
     ntpw.EventTree()->Branch("skNCutsPassedCumul", &skNCutsPassedCumul, "skNCutsPassedCumul/I");
  brNCutsPassedCumul->SetAutoDelete(kFALSE);
  int skMode = 0;
  TBranch * brMode = 
     ntpw.EventTree()->Branch("skMode", &skMode, "skMode/I");
  brMode->SetAutoDelete(kFALSE);
  double skEnu = 0;
  TBranch * brEnu = 
     ntpw.EventTree()->Branch("skEnu", &skEnu, "skEnu/D");
  brEnu->SetAutoDelete(kFALSE);
  double skEvis = 0;
  TBranch * brEvis =
    ntpw.EventTree()->Branch("skEvis", &skEvis, "skEvis/D");
  brEvis->SetAutoDelete(kFALSE);
  int skCutLevel = 0;
  TBranch * brCutLevel = 
     ntpw.EventTree()->Branch("skCutLevel", &skCutLevel, "skCutLevel/I");
  brCutLevel->SetAutoDelete(kFALSE);

  // open event list file
  TFile fevlist(gOptEventList.c_str(),"read");

  // get the event list
  TTree * evlisttree = (TTree*) fevlist.Get("selection_tree");
  if(!evlisttree){
      LOG("gSKMergeData", pFATAL)
        << "Cannot find event list tree - Exiting";
      gAbortingInErr = true;
      exit(1);    
  }

  // set the branches for the runnumber,ghep event number, 
  // and other SK reco info
  int   _RunNum;
  int   _EvtNumGHEP; 
  int   _NCutsPassedCumulative;
  int   _CutLevel;
  int   _Mode;
  float _EnuReco;
  float _Enu;
  float _Evis;

  evlisttree -> SetBranchAddress ("RunNum",                &_RunNum               );
  evlisttree -> SetBranchAddress ("EvtNumGHEP",            &_EvtNumGHEP           );
  evlisttree -> SetBranchAddress ("EnuReco",               &_EnuReco              );
  evlisttree -> SetBranchAddress ("NCutsPassedCumulative", &_NCutsPassedCumulative);
  evlisttree -> SetBranchAddress ("CutLevel",              &_CutLevel             );
  evlisttree -> SetBranchAddress ("Mode",                  &_Mode                 );
  evlisttree -> SetBranchAddress ("Enu",                   &_Enu                  );
  evlisttree -> SetBranchAddress ("Evis",                  &_Evis                 );

  // loop over the list and extract each GENIE event from its 
  // original GHEP files
  int nev = evlisttree->GetEntries();
  for(int iev=0; iev<nev; iev++) {
    _RunNum = -1;
    evlisttree->GetEntry(iev);

    int irunnu = _RunNum; // run number
    int ievnu  = _EvtNumGHEP-1; // event number in given run

    EventRecord * event = GetEvent(irunnu,ievnu);
    assert(event);

    // copy info in the extra branches
    skNCutsPassedCumul = _NCutsPassedCumulative;
    skCutLevel         = _CutLevel;
    skMode             = _Mode;
    skEnu              = _Enu     * units::MeV; // SK E's in MeV. Convert to GeV.
    skEnuReco          = _EnuReco * units::MeV;
    skEvis             = _Evis    * units::MeV;

    // Check that the Enu from sk ntp matches that of GHEP record. Could
    // also do the same for Mode.
    double EnuGHEP = event->Probe()->Energy();
    int ModeGHEP = utils::ghep::NeutReactionCode(event); 

    bool same_enu = TMath::Abs(EnuGHEP-skEnu) < 0.01; 
    bool same_mode = skMode == ModeGHEP;
 
    LOG("gSKMergeData", pDEBUG) 
       << "SKDETSIM: E = " << skEnu   << ", NeutMode = " << skMode;
    LOG("gSKMergeData", pDEBUG) 
       << "GHEP:     E = " << EnuGHEP << ", NeutMode = " << ModeGHEP;

    if(!same_enu || !same_mode){
      LOG("gSKMergeData", pFATAL) 
         << " **** Mismatch between SK event and GHEP record. ****";
      LOG("gSKMergeData", pFATAL) 
         << "SKDETSIM: E = " << skEnu   << ", NeutMode = " << skMode;
      LOG("gSKMergeData", pFATAL) 
         << "GHEP:     E = " << EnuGHEP << ", NeutMode = " << ModeGHEP;
      gAbortingInErr = true;
      exit(1);
    } 

    // copy to output file
    ntpw.AddEventRecord(iev,event);

    delete event;
  }
               
  ntpw.Save();

  if(evlisttree){
    delete evlisttree; 
    evlisttree = 0;
  } 

  LOG("gSKMergeData", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
EventRecord * GetEvent(int irunnu, int ievnu)
{
//  const char filename[1000] = {'\0'};
//  Form(filename,gOptFilenmTemplate.c_str(),irunnu);

  TString filename;
  filename.Form(gOptFilenmTemplate.c_str(), irunnu);

  LOG("gSKMergeData", pNOTICE) 
      << "File template : "<< gOptFilenmTemplate 
      << ", Run nu. : "<< irunnu;

  LOG("gSKMergeData", pNOTICE) << "Opening file: "<< filename;
  TFile file(filename.Data(),"READ");

  TTree * tree = dynamic_cast <TTree *> ( file.Get("gtree")  );
  if(!tree) return 0;

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  tree->GetEntry(ievnu);
  EventRecord & event = *(mcrec->event);

  EventRecord * event_copy = new EventRecord(event);

  mcrec->Clear();
  file.Close();

  return event_copy;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gSKMergeData", pINFO) << "Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get file with events passing cuts + SKDETSIM info
  if( parser.OptionExists('l') ) {
    LOG("gSKMergeData", pINFO) 
       << "Reading event list";
    gOptEventList = parser.ArgAsString('l');
  } else {
    LOG("gSKMergeData", pFATAL) 
       << "Unspecified event list - Exiting";
    gAbortingInErr = true;
    exit(1);
  }

  // GHEP event file template
  if( parser.OptionExists('t') ) {
    LOG("gSKMergeData", pINFO) 
       << "Reading GHEP event filename template";
    gOptFilenmTemplate = parser.ArgAsString('t');
  } else {
    LOG("gSKMergeData", pFATAL) 
        << "Unspecified event filename template - Exiting";
    gAbortingInErr = true;
    exit(1);
  }

}
//_________________________________________________________________________________
