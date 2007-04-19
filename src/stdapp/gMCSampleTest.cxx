//____________________________________________________________________________
/*!

\program gmctest

\brief   A GENIE utility that reads a generated event sample and generates a 
         ROOT file with loads of histograms with basic MC truth quantities.
         Typically used in validating new generator releases.

         Syntax :
           gmctest -f filename -t tgt [-n events] [-c another_filename]

         Options:
           [] Denotes an optional argument
           -f Specifies the GENIE/ROOT file with the generated event sample
           -t Target pdg code (analyzing interactions for one target at a time)
           -n Specifies how many events to analyze [default: all]
	   -c Specifies another GENIE/ROOT event sample file for comparison 

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 02, 2005

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TPostScript.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>

#include "Conventions/Constants.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpMCSummary.h"
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
void   GetCommandLineArgs(int argc, char ** argv);
void   PrintSyntax          (void);
bool   CheckRootFilename    (string filename);
string OutputFileName       (string input_file_name);
void   AnalyzeSample        (string filename);
void   Initialize           (void);
void   EventLoop            (void);
void   SaveResults          (void);
void   CleanUp              (void);
void   AddKineDir           (TFile & file);
void   AddKineCcNumuDir     (TFile & file);
void   AddKineCcNumuQelDir  (TFile & file);
void   AddKineCcNumuResDir  (TFile & file);
void   AddKineCcNumuDisDir  (TFile & file);
void   AddHMultDir          (TFile & file);
void   AddVtxDir            (TFile & file);
void   AddNtpDir            (TFile & file);
void   Plot                 (void);
void   PlotH1F              (string name, string title);

// command-line arguments
Long64_t gOptN;                // (-n)  process so many events, all if -1
Long64_t gOptTgtPdgC;          // (-t) process events only for this nucl. target
string   gOptInpFile;          // (-f) input GENIE event sample file
string   gOptInpTemplFile;     // (-c) input GENIE event sample file (template)

//
bool               gSampleComp = false;  // run sample comparisons ?
TCanvas *          gC = 0;
TPostScript *      gPS = 0;
TDirectory *       gTestedSampleDir = 0; // plot dir in test sample file
TDirectory *       gTempltSampleDir = 0; // plot dir in template sample file
NtpMCEventRecord * gMCRec = 0;
TTree *            gEventRecTree = 0;
Long64_t           gNEvt = 0;
string             gCurrInpFilename = "";
TFile *            gCurrInpFile = 0;

// summary tree
TTree * tEvtTree;
int    br_iev       = 0;
int    br_neutrino  = 0;
int    br_target    = 0;
int    br_hitnuc    = 0;
int    br_hitqrk    = 0;
bool   br_qel       = false;
bool   br_res       = false;
bool   br_dis       = false;
bool   br_coh       = false;
bool   br_em        = false;
bool   br_weakcc    = false;
bool   br_weaknc    = false;
double br_kine_xs   = 0;
double br_kine_ys   = 0; 
double br_kine_ts   = 0; 
double br_kine_Q2s  = 0; 
double br_kine_Ws   = 0;
double br_kine_x    = 0;
double br_kine_y    = 0; 
double br_kine_t    = 0; 
double br_kine_Q2   = 0; 
double br_kine_W    = 0;
double br_kine_v    = 0;
double br_Ev        = 0;
double br_El        = 0;
double br_vtxx      = 0;
double br_vtxy      = 0;
double br_vtxz      = 0;
double br_weight    = 0;
int    br_np        = 0;
int    br_nn        = 0;
int    br_npip      = 0;
int    br_npim      = 0;
int    br_npi0      = 0;
int    br_ngamma    = 0;
int    br_nKp       = 0;
int    br_nKm       = 0;
int    br_nK0       = 0;

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments 
  GetCommandLineArgs(argc,argv);

  //-- analyze tested event generation sample
  assert(CheckRootFilename(gOptInpFile));
  AnalyzeSample(gOptInpFile);

  //-- analyze tested event generation sample template - if available
  if (CheckRootFilename(gOptInpTemplFile) ) {
    AnalyzeSample(gOptInpTemplFile);
    gSampleComp = true;
  }

  //-- plot results
  //-- if a template file was available plot comparisons too
  Plot();

  LOG("gmctest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void AnalyzeSample(string filename)
{
  gCurrInpFilename = filename;

  Initialize  ();
  EventLoop   ();
  SaveResults ();
  CleanUp     ();

  LOG("gmctest", pINFO)  << "Done analyzing : " << filename;
}
//_________________________________________________________________________________
void Initialize(void)
{
  //-- open the ROOT file and get the TTree & its header

  LOG("gmctest", pNOTICE) << "*** Opening GHEP data file: " << gCurrInpFilename;

  TFile * gCurrInpFile = new TFile(gCurrInpFilename.c_str(),"READ");
  TTree * tree = dynamic_cast <TTree *> (gCurrInpFile->Get("gtree"));
  NtpMCTreeHeader * thdr = 
       dynamic_cast <NtpMCTreeHeader *> (gCurrInpFile->Get("header"));

  LOG("gmctest", pNOTICE) << "*** Input tree header: " << *thdr;

  //-- figure out the TTree format (GENIE supports multiple formats)
  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord);

  //-- set the event record branch ptr
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       tree->GetEntries() : TMath::Min( tree->GetEntries(), gOptN );

  //-- keep what is needed for the EventLoop()
  gEventRecTree = tree;
  gMCRec        = mcrec;
  gNEvt         = nmax;

  //-- create summary tree
  tEvtTree = new TTree("tEvtTree","event tree summary");
  tEvtTree->Branch("iev",    &br_iev,       "iev/I"      );
  tEvtTree->Branch("neu",    &br_neutrino,  "neu/I"      );
  tEvtTree->Branch("tgt" ,   &br_target,    "tgt/I"      );
  tEvtTree->Branch("hitnuc", &br_hitnuc,    "hitnuc/I"   );
  tEvtTree->Branch("hitqrk", &br_hitqrk,    "hitqrk/I"   );
  tEvtTree->Branch("qel",    &br_qel,       "br_qel/O"   );
  tEvtTree->Branch("res",    &br_res,       "br_res/O"   );
  tEvtTree->Branch("dis",    &br_dis,       "br_dis/O"   );
  tEvtTree->Branch("coh",    &br_coh,       "br_coh/O"   );
  tEvtTree->Branch("em",     &br_em,        "br_em/O"    );
  tEvtTree->Branch("weakcc", &br_weakcc,    "br_weakcc/O");
  tEvtTree->Branch("weaknc", &br_weaknc,    "br_weaknc/O");
  tEvtTree->Branch("xs",     &br_kine_xs,   "xs/D"       );
  tEvtTree->Branch("ys",     &br_kine_ys,   "ys/D"       );
  tEvtTree->Branch("ts",     &br_kine_ts,   "ts/D"       );
  tEvtTree->Branch("Q2s",    &br_kine_Q2s,  "Q2s/D"      );
  tEvtTree->Branch("Ws",     &br_kine_Ws,   "Ws/D"       );
  tEvtTree->Branch("x",      &br_kine_x,    "x/D"        );
  tEvtTree->Branch("y",      &br_kine_y,    "y/D"        );
  tEvtTree->Branch("t",      &br_kine_t,    "t/D"        );
  tEvtTree->Branch("Q2",     &br_kine_Q2,   "Q2/D"       );
  tEvtTree->Branch("W",      &br_kine_W,    "W/D"        );
  tEvtTree->Branch("v",      &br_kine_v,    "v/D"        );
  tEvtTree->Branch("Ev",     &br_Ev,        "Ev/D"       );
  tEvtTree->Branch("El",     &br_El,        "El/D"       );
  tEvtTree->Branch("vtxx",   &br_vtxx,      "vtxx/D"     );
  tEvtTree->Branch("vtxy",   &br_vtxy,      "vtxy/D"     );
  tEvtTree->Branch("vtxz",   &br_vtxz,      "vtxz/D"     );
  tEvtTree->Branch("wgt",    &br_weight,    "wgt/D"      );
  tEvtTree->Branch("np",     &br_np,        "np/I"       );
  tEvtTree->Branch("nn",     &br_nn,        "nn/I"       );
  tEvtTree->Branch("npip",   &br_npip,      "npip/I"     );
  tEvtTree->Branch("npim",   &br_npim,      "npim/I"     );
  tEvtTree->Branch("npi0",   &br_npi0,      "npi0/I"     );
  tEvtTree->Branch("ngamma", &br_ngamma,    "ngamma/I"   );
  tEvtTree->Branch("nKp",    &br_nKp,       "nKp/I"      );
  tEvtTree->Branch("nKm",    &br_nKm,       "nKm/I"      );
  tEvtTree->Branch("nK0",    &br_nK0,       "nK0/I"      );
}
//_________________________________________________________________________________
void EventLoop(void)
{
  LOG("gmctest", pNOTICE) << "*** Analyzing: " << gNEvt << " events";

  if ( gNEvt<0 )       return;
  if ( !gEventRecTree) return;
  if ( !gMCRec )       return;

  for(Long64_t i = 0; i < gNEvt; i++) {
    gEventRecTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    // LOG("gmctest", pINFO) << rec_header;
    // LOG("gmctest", pINFO) << event;

    // go further only if the event is physical
    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) continue;

    // nuclear target or free nucleon target
    GHepParticle * target = event.Particle(1);
    assert(target);

    // only further only if it matches the requested target
    if(target->Pdg() != gOptTgtPdgC) continue;

    // neutrino
    GHepParticle * neutrino = event.Probe();
    assert(neutrino);
    // final state primary lepton
    GHepParticle * fsl = event.FinalStatePrimaryLepton();
    assert(fsl);
    // hit nucleon
    GHepParticle * hitnucl = event.HitNucleon();
    if(!hitnucl) continue;

    const Interaction * interaction = event.Summary();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    const Kinematics &   kine       = interaction->Kine();

    double weight = event.Weight();
    double Ev     = neutrino->Energy();
    double El     = fsl->Energy();

    const TLorentzVector & k1 = *(neutrino->P4()); // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());      // l 4-p (k2)
    const TLorentzVector & p1 = *(hitnucl->P4());  // N 4-p (p1)

    double M  = kNucleonMass;

    // compute kinematical params
    // (see, for example, particle data booklet, page 210-211)

    TLorentzVector q  = k1-k2;    // q=k1-k2, 4-p transfer
    double v  = (q*p1)/M;         // v (E transfer in hit nucleon rest frame)
    double Q2 = -1 * q.M2();      // momemtum transfer
    double x  = 0.5*Q2/(M*v);     // Bjorken x
    double y  = v*M/(k1*p1);      // Inelasticity, y = q*P1/k1*P1
    double W2 = M*M + 2*M*v - Q2; // Hadronic Invariant mass ^ 2
    double W  = TMath::Sqrt(W2);

    // also, access kinematical params _exactly_ as they were selected
    // by the event generation code (can be different because of the
    // off-shell kinematics)

    bool get_selected = true;
    double xs  = kine.x (get_selected);
    double ys  = kine.y (get_selected);
    double ts  = 0; //kine.t (get_selected);
    double Q2s = kine.Q2(get_selected);
    double Ws  = kine.W (get_selected);

    // get interaction vertex (in detector coord system)
    //
    TLorentzVector * vtx = event.Vertex(); 

    // fill the summary ntuple
    //
    br_iev       = i;
    br_neutrino  = neutrino->Pdg();
    br_target    = 0;
    br_hitnuc    = 0;
    br_hitqrk    = 0;
    br_qel       = proc_info.IsQuasiElastic();
    br_res       = proc_info.IsResonant();
    br_dis       = proc_info.IsDeepInelastic();
    br_coh       = proc_info.IsCoherent();
    br_em        = proc_info.IsEM();
    br_weakcc    = proc_info.IsWeakCC();
    br_weaknc    = proc_info.IsWeakNC();
    br_kine_xs   = xs;
    br_kine_ys   = ys; 
    br_kine_ts   = ts; 
    br_kine_Q2s  = Q2s; 
    br_kine_Ws   = Ws;
    br_kine_x    = x;
    br_kine_y    = y; 
    br_kine_t    = 0; 
    br_kine_Q2   = Q2; 
    br_kine_W    = W;
    br_kine_v    = v;
    br_Ev        = Ev;
    br_El        = El;
    br_vtxx      = vtx->X();
    br_vtxy      = vtx->Y();
    br_vtxz      = vtx->Z();
    br_weight    = weight; 
    br_np        = event.NEntries(kPdgProton,  kIStStableFinalState);
    br_nn        = event.NEntries(kPdgNeutron, kIStStableFinalState);
    br_npip      = event.NEntries(kPdgPiP,     kIStStableFinalState);
    br_npim      = event.NEntries(kPdgPiM,     kIStStableFinalState);
    br_npi0      = event.NEntries(kPdgPi0,     kIStStableFinalState);
    br_ngamma    = event.NEntries(kPdgGamma,   kIStStableFinalState);
    br_nKp       = event.NEntries(kPdgKP,      kIStStableFinalState);
    br_nKm       = event.NEntries(kPdgKM,      kIStStableFinalState);
    br_nK0       = event.NEntries(kPdgK0,      kIStStableFinalState);

    tEvtTree->Fill();

    gMCRec->Clear();

  } // event loop
}
//_________________________________________________________________________________
void CleanUp(void)
{
  //-- close the input GENIE event file
  LOG("gmctest", pNOTICE) << "*** Closing input GHEP data stream";
//  gCurrInpFile->Close();
  delete gCurrInpFile;
  gCurrInpFile=0;
  gMCRec=0;
  gEventRecTree=0;
}
//_________________________________________________________________________________
void SaveResults(void)
{
  //-- build output filename based on input filename
  string outpfilename = OutputFileName(gCurrInpFilename);
  LOG("gmctest", pNOTICE) 
          << "*** Saving output histograms to: " << outpfilename;

  LOG("gmctest", pNOTICE) 
     << "number of events processed: " << tEvtTree->GetEntries();

  TFile outfile(outpfilename.c_str(),"recreate");

  AddKineDir           (outfile);
  AddKineCcNumuDir     (outfile);
  AddKineCcNumuQelDir  (outfile);
  AddKineCcNumuResDir  (outfile);
  AddKineCcNumuDisDir  (outfile);
  AddHMultDir          (outfile);
  AddVtxDir            (outfile);
  AddNtpDir            (outfile);

  outfile.Close();
}
//_________________________________________________________________________________
void AddKineDir(TFile & file)
{
  TDirectory * KineDir = file.mkdir("KineDir", "Kinematics plots - All events");
  KineDir->cd();

  tEvtTree->Draw( "x>>hx",     "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "y>>hy",     "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "t>>ht",     "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "Q2>>hQ2",   "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "W>>hW",     "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "v>>hv",     "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "xs>>hxs",   "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "ys>>hys",   "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "ts>>hts",   "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "Q2s>>hQ2s", "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "Ws>>hWs",   "wgt*(1==1)", "GOFF");

  KineDir->Write();
}
//_________________________________________________________________________________
void AddKineCcNumuDir(TFile & file)
{
  TDirectory * KineCcNumuDir = file.mkdir(
         "KineCcNumuDir", "Kinematics plots - All nu_mu CC events");
  KineCcNumuDir->cd();

  tEvtTree->Draw( "x>>hx",     "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "y>>hy",     "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "t>>ht",     "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "Q2>>hQ2",   "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "W>>hW",     "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "v>>hv",     "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "xs>>hxs",   "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "ys>>hys",   "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "ts>>hts",   "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "Q2s>>hQ2s", "wgt*(neu==14&&weakcc)", "GOFF");
  tEvtTree->Draw( "Ws>>hWs",   "wgt*(neu==14&&weakcc)", "GOFF");

  KineCcNumuDir->Write();
}
//_________________________________________________________________________________
void AddKineCcNumuQelDir(TFile & file)
{
  TDirectory * KineCcNumuQelDir = file.mkdir(
         "KineCcNumuQelDir", "Kinematics plots - nu_mu CC QEL");
  KineCcNumuQelDir->cd();

  tEvtTree->Draw( "x>>hx",     "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "y>>hy",     "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "t>>ht",     "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "Q2>>hQ2",   "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "W>>hW",     "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "v>>hv",     "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "xs>>hxs",   "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "ys>>hys",   "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "ts>>hts",   "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "Q2s>>hQ2s", "wgt*(neu==14&&weakcc&&qel)", "GOFF");
  tEvtTree->Draw( "Ws>>hWs",   "wgt*(neu==14&&weakcc&&qel)", "GOFF");

  KineCcNumuQelDir->Write();
}
//_________________________________________________________________________________
void AddKineCcNumuResDir(TFile & file)
{
  TDirectory * KineCcNumuResDir = file.mkdir(
         "KineCcNumuResDir", "Kinematics plots - nu_mu CC RES");
  KineCcNumuResDir->cd();

  tEvtTree->Draw( "x>>hx",     "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "y>>hy",     "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "t>>ht",     "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "Q2>>hQ2",   "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "W>>hW",     "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "v>>hv",     "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "xs>>hxs",   "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "ys>>hys",   "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "ts>>hts",   "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "Q2s>>hQ2s", "wgt*(neu==14&&weakcc&&res)", "GOFF");
  tEvtTree->Draw( "Ws>>hWs",   "wgt*(neu==14&&weakcc&&res)", "GOFF");

  KineCcNumuResDir->Write();
}
//_________________________________________________________________________________
void AddKineCcNumuDisDir(TFile & file)
{
  TDirectory * KineCcNumuDisDir = file.mkdir(
         "KineCcNumuDisDir", "Kinematics plots - nu_mu CC DIS");
  KineCcNumuDisDir->cd();

  tEvtTree->Draw( "x>>hx",     "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "y>>hy",     "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "t>>ht",     "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "Q2>>hQ2",   "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "W>>hW",     "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "v>>hv",     "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "xs>>hxs",   "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "ys>>hys",   "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "ts>>hts",   "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "Q2s>>hQ2s", "wgt*(neu==14&&weakcc&&dis)", "GOFF");
  tEvtTree->Draw( "Ws>>hWs",   "wgt*(neu==14&&weakcc&&dis)", "GOFF");

  KineCcNumuDisDir->Write();
}
//_________________________________________________________________________________
void AddHMultDir(TFile & file)
{
  TDirectory * HMultDir = file.mkdir(
              "HMultDir", "Hadronic multiplicities - All processes");
  HMultDir->cd();

  tEvtTree->Draw( "np>>hnp",         "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "nn>>hnn",         "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "npip>>hnpip",     "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "npim>>hnpim",     "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "ngamma/2>>hnpi0", "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "nKp>>hnKp",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "nKm>>hnKm",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw( "nK0>>hnK0",       "wgt*(1==1)", "GOFF");

  HMultDir->Write();
}
//_________________________________________________________________________________
void AddVtxDir(TFile & file)
{
  TDirectory * VtxDir = file.mkdir("VtxDir", "Event vertex plots");
  VtxDir->cd();

  tEvtTree->Draw("vtxx>>hvtxx",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw("vtxy>>hvtxy",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw("vtxz>>hvtxz",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw("vtxx:vtxy>>hvtxxy", "wgt*(1==1)", "GOFF");

  VtxDir->Write();
}
//_________________________________________________________________________________
void AddNtpDir(TFile & file)
{
  TDirectory * NtpDir = file.mkdir("NtpDir", "Summary Ntuples");
  NtpDir->cd();
  NtpDir->Add(tEvtTree);
  NtpDir->Write();
}
//_________________________________________________________________________________
void Plot(void)
{
  TFile * file1 = new TFile(OutputFileName(gOptInpFile).c_str(), "read");
  TFile * file2 = (gSampleComp) ?
                  new TFile(OutputFileName(gOptInpTemplFile).c_str(), "read") : 0;

  gC = new TCanvas("gC");
  gC->SetFillColor(0);
  gC->SetBorderMode(0);

  gPS = new TPostScript("out.ps", 112);

  // --- front page

  gPS->NewPage();
  gC->cd();

  TPavesText title(0.1, 0.6, 0.9, 0.9, 0, "tr");
  title.SetTextSize(0.04);
  title.AddText("GENIE MC Sample Test");
  title.AddText(Form("Testing all events for target pdg = %d", gOptTgtPdgC));
  title.SetFillColor(46);  
  title.SetTextColor(10);  
  title.Draw();

  TPavesText stitle(0.1, 0.2, 0.9, 0.5, 0, "tr");
  stitle.SetTextSize(0.027);
  stitle.AddText(Form("Tested event sample : %s", gOptInpFile.c_str()));
  if(gSampleComp) 
    stitle.AddText(Form("Event sample tempate : %s", gOptInpTemplFile.c_str()));
  else
    stitle.AddText("Event sample tempate : not available");
  stitle.SetFillColor(0);  
  stitle.SetTextColor(1);  
  stitle.Draw();

  gC->Update();

  // --- kinematics - all events in sample (directory: KineDir)

  gTestedSampleDir = (TDirectory*) file1->Get("KineDir");
  gTempltSampleDir = (gSampleComp) ? (TDirectory*) file2->Get("KineDir") : 0;

  PlotH1F ( "hx",   "x_{comp} - all events" );
  PlotH1F ( "hy",   "y_{comp} - all events" );
  PlotH1F ( "ht",   "t_{comp} - all events" );
  PlotH1F ( "hW",   "W_{comp} - all events" );
  PlotH1F ( "hQ2",  "Q^{2}_{comp} (GeV^{2}) - all events" );
  PlotH1F ( "hv",   "v_{comp} - all events" );
  PlotH1F ( "hxs",  "x_{sel} - all events" );
  PlotH1F ( "hys",  "y_{sel} - all events" );
  PlotH1F ( "hts",  "t_{sel} - all events" );
  PlotH1F ( "hWs",  "W_{sel} - all events" );
  PlotH1F ( "hQ2s", "Q^{2}_{sel} (GeV^{2}) - all events" );

  // --- kinematics - numu CC (directory: KineCcNumuDir)

  gTestedSampleDir = (TDirectory*) file1->Get("KineCcNumuDir");
  gTempltSampleDir = (gSampleComp) ? (TDirectory*) file2->Get("KineCcNumuDir") : 0;

  PlotH1F ( "hx",   "x_{comp} - #nu_{#mu} CC" );
  PlotH1F ( "hy",   "y_{comp} - #nu_{#mu} CC" );
  PlotH1F ( "ht",   "t_{comp} - #nu_{#mu} CC" );
  PlotH1F ( "hW",   "W_{comp} - #nu_{#mu} CC" );
  PlotH1F ( "hQ2",  "Q^{2}_{comp} (GeV^{2}) - #nu_{#mu} CC" );
  PlotH1F ( "hv",   "v_{comp} - #nu_{#mu} CC" );
  PlotH1F ( "hxs",  "x_{sel} - #nu_{#mu} CC" );
  PlotH1F ( "hys",  "y_{sel} - #nu_{#mu} CC" );
  PlotH1F ( "hts",  "t_{sel} - #nu_{#mu} CC" );
  PlotH1F ( "hWs",  "W_{sel} - #nu_{#mu} CC" );
  PlotH1F ( "hQ2s", "Q^{2}_{sel} (GeV^{2}) - #nu_{#mu} CC" );

  // --- kinematics - numu CC QEL (directory: KineCcNumuQelDir)

  gTestedSampleDir = (TDirectory*) file1->Get("KineCcNumuQelDir");
  gTempltSampleDir = (gSampleComp) ? (TDirectory*) file2->Get("KineCcNumuQelDir") : 0;

  PlotH1F ( "hx",   "x_{comp} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hy",   "y_{comp} - #nu_{#mu} CC QEL" );
  PlotH1F ( "ht",   "t_{comp} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hW",   "W_{comp} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hQ2",  "Q^{2}_{comp} (GeV^{2}) - #nu_{#mu} CC QEL" );
  PlotH1F ( "hv",   "v_{comp} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hxs",  "x_{sel} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hys",  "y_{sel} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hts",  "t_{sel} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hWs",  "W_{sel} - #nu_{#mu} CC QEL" );
  PlotH1F ( "hQ2s", "Q^{2}_{sel} (GeV^{2}) - #nu_{#mu} CC QEL" );

  // --- kinematics - numu CC RES (directory: KineCcNumuResDir)

  gTestedSampleDir = (TDirectory*) file1->Get("KineCcNumuResDir");
  gTempltSampleDir = (gSampleComp) ? (TDirectory*) file2->Get("KineCcNumuResDir") : 0;

  PlotH1F ( "hx",   "x_{comp} - #nu_{#mu} CC RES" );
  PlotH1F ( "hy",   "y_{comp} - #nu_{#mu} CC RES" );
  PlotH1F ( "ht",   "t_{comp} - #nu_{#mu} CC RES" );
  PlotH1F ( "hW",   "W_{comp} - #nu_{#mu} CC RES" );
  PlotH1F ( "hQ2",  "Q^{2}_{comp} (GeV^{2}) - #nu_{#mu} CC RES" );
  PlotH1F ( "hv",   "v_{comp} - #nu_{#mu} CC RES" );
  PlotH1F ( "hxs",  "x_{sel} - #nu_{#mu} CC RES" );
  PlotH1F ( "hys",  "y_{sel} - #nu_{#mu} CC RES" );
  PlotH1F ( "hts",  "t_{sel} - #nu_{#mu} CC RES" );
  PlotH1F ( "hWs",  "W_{sel} - #nu_{#mu} CC RES" );
  PlotH1F ( "hQ2s", "Q^{2}_{sel} (GeV^{2}) - #nu_{#mu} CC RES" );

  // --- kinematics - numu CC DIS (directory: KineCcNumuDisDir)

  gTestedSampleDir = (TDirectory*) file1->Get("KineCcNumuDisDir");
  gTempltSampleDir = (gSampleComp) ? (TDirectory*) file2->Get("KineCcNumuDisDir") : 0;

  PlotH1F ( "hx",   "x_{comp} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hy",   "y_{comp} - #nu_{#mu} CC DIS" );
  PlotH1F ( "ht",   "t_{comp} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hW",   "W_{comp} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hQ2",  "Q^{2}_{comp} (GeV^{2}) - #nu_{#mu} CC DIS" );
  PlotH1F ( "hv",   "v_{comp} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hxs",  "x_{sel} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hys",  "y_{sel} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hts",  "t_{sel} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hWs",  "W_{sel} - #nu_{#mu} CC DIS" );
  PlotH1F ( "hQ2s", "Q^{2}_{sel} (GeV^{2}) - #nu_{#mu} CC DIS" );

  gPS->Close();

  file1->Close();
  delete file1;
  if(file2) {
    file2->Close();
    delete file2;
  }
  delete gC;
  delete gPS;
}
//_________________________________________________________________________________
void PlotH1F(string name, string title)
{
  gPS->NewPage();
  gC->cd();

  TH1F * tested_hst = dynamic_cast<TH1F *> (gTestedSampleDir->Get(name.c_str()));
  TH1F * templt_hst = (gSampleComp) ?
                      dynamic_cast<TH1F *> (gTempltSampleDir->Get(name.c_str())) : 0;

  if(!tested_hst) return;

  tested_hst->SetLineColor(2);
  tested_hst->SetLineWidth(2);
  tested_hst->Draw();

  if(templt_hst) {
    templt_hst->SetLineWidth(2);
    templt_hst->SetMarkerSize(1.3);
    templt_hst->SetMarkerStyle(8);
    templt_hst->Draw("PERRSAME");
  }

  tested_hst->GetXaxis()->SetTitle(title.c_str());
  gC->Update();
}
//_________________________________________________________________________________
string OutputFileName(string inpname)
{
// Builds the output filename based on the name of the input filename
// Perfors the following conversion: name.root -> name.gmctest.root

  unsigned int L = inpname.length();

  // if the last 4 characters are "root" (ROOT file extension) then
  // remove them
  if(inpname.substr(L-4, L).find("root") != string::npos) {
    inpname.erase(L-4, L);
  }
  ostringstream name;
  name << inpname << "gmctest-" << gOptTgtPdgC << ".root";

  return gSystem->BaseName(name.str().c_str());
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gmctest", pNOTICE) << "*** Parsing commad line arguments";

   //number of events:
  try {
    LOG("gmctest", pINFO) << "Reading number of events to analyze";
    gOptN = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctest", pINFO)
              << "Unspecified number of events to analyze - Use all";
      gOptN = -1;
    }
  }

  //target PDG code:
  try {
    LOG("gmctest", pINFO) << "Reading target PDG code";
    gOptTgtPdgC = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctest", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //get GENIE event sample ROOT file (ER-format)
  try {
    LOG("gmctest", pINFO) << "Reading event sample filename";
    gOptInpFile = utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctst", pFATAL) << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //get another (template) GENIE event sample ROOT file (ER-format)
  try {
    LOG("gmctest", pINFO) << "Reading filename of event sample template";
    gOptInpTemplFile = utils::clap::CmdLineArgAsString(argc,argv,'c');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctst", pNOTICE) 
         << "Unspecified event sample template - WIll not run comparisons";
    }
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmctest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gmctest -f file -t tgtpdg [-n nev] [-c template_file]\n";
}
//_________________________________________________________________________________
bool CheckRootFilename(string filename)
{
  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gmctest", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//_________________________________________________________________________________

