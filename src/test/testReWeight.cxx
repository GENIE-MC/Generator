//____________________________________________________________________________
/*!

\program testReWeight

\brief   A simple test program to illustrate how to use the GENIE event
         reweighting package

         To run: testReWeight -f filename
         where the filename points to a ROOT file containing a GENIE output
         TTree (ER)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 25, 2005
*/
//____________________________________________________________________________

#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>

#include "Algorithm/AlgId.h"
#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "ReWeight/MCModel.h"
#include "ReWeight/WeightCalculator.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using namespace genie;

//func prototypes
void    GetCommandLineArgs (int argc, char ** argv);
void    PrintSyntax        (void);
MCModel CreateOldXSecModel (void);
MCModel CreateNewXSecModel (void);

//input options (from command line arguments):
string gOptInpFile;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);
  LOG("test", pNOTICE) << "Input  GENIE/ROOT file: " << gOptInpFile;

  //-- define the 'old' & 'new' MCModel
  MCModel old_model = CreateOldXSecModel();
  MCModel new_model = CreateOldXSecModel();

  //-- create a weight calculator
  WeightCalculator wcalc;
  wcalc.OldCrossSectionModel(old_model);
  wcalc.NewCrossSectionModel(new_model);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFile.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("test", pINFO) << "Input tree header: " << *thdr;

  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord); // only ER trees in this test

  //-- The ER ntuple contains a single TBranch with NtpMCEventRecord
  //   objects in its leaves. Set the branch address for reading it.
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- loop over TTree NtpMC records
  for(int i = 0; i< tree->GetEntries(); i++) {

    tree->GetEntry(i);

    // get/print the next event and mc record header
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("test", pINFO) << rec_header;
    LOG("test", pINFO) << event;

    // reweight the event

    LOG("test", pINFO) << "Got the event. Computing new weight";

    double wght = wcalc.ReWeight(event);

    LOG("test", pINFO)  << "Old weight = " << event.GetWeight();
    LOG("test", pINFO)  << "New weight = " << wght;
  }

  LOG("test", pINFO)  << "Done!";
  return 0;
}
//___________________________________________________________________
MCModel CreateOldXSecModel(void)
{
  MCModel     model;
  ProcessInfo proc;
  AlgId       alg;

  proc.Set(kScQuasiElastic,  kIntWeakCC);
  alg.SetId("genie::NuclearTargetXSec","QEL-CC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScQuasiElastic,  kIntWeakNC);
  alg.SetId("genie::NuclearTargetXSec","QEL-NC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScDeepInelastic, kIntWeakCC);
  alg.SetId("genie::NuclearTargetXSec","DIS-CC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScDeepInelastic, kIntWeakNC);
  alg.SetId("genie::NuclearTargetXSec","DIS-NC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScResonant, kIntWeakCC);
  alg.SetId("genie::NuclearTargetXSec","RES-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScResonant, kIntWeakNC);
  alg.SetId("genie::NuclearTargetXSec","RES-Default");
  model.UseXSecAlg(proc,alg);

  return model;
}
//___________________________________________________________________
MCModel CreateNewXSecModel(void)
{
  MCModel     model;
  ProcessInfo proc;
  AlgId       alg;

  proc.Set(kScQuasiElastic,  kIntWeakCC);
  alg.SetId("genie::NuclearTargetXSec","QEL-CC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScQuasiElastic,  kIntWeakNC);
  alg.SetId("genie::NuclearTargetXSec","QEL-NC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScDeepInelastic, kIntWeakCC);
  alg.SetId("genie::NuclearTargetXSec","DIS-CC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScDeepInelastic, kIntWeakNC);
  alg.SetId("genie::NuclearTargetXSec","DIS-NC-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScResonant, kIntWeakCC);
  alg.SetId("genie::NuclearTargetXSec","RES-Default");
  model.UseXSecAlg(proc,alg);

  proc.Set(kScResonant, kIntWeakNC);
  alg.SetId("genie::NuclearTargetXSec","RES-Default");
  model.UseXSecAlg(proc,alg);

  return model;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  //get input ROOT file (containing a GENIE ER ntuple)
  try {
    LOG("test", pINFO) << "Reading input filename";
    gOptInpFile = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("test", pFATAL)
             << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // check input GENIE ROOT file
  bool ok = !(gSystem->AccessPathName(gOptInpFile.c_str()));
  if (!ok) {
    LOG("test", pFATAL)
      << "Input ROOT file [" << gOptInpFile << "] is not accessible";
    exit(2);
  }
}
//___________________________________________________________________
void PrintSyntax(void)
{
  LOG("test", pNOTICE)
    << "\n\n" << "Syntax:" << "\n   test -f input_filename\n";
}
//___________________________________________________________________
