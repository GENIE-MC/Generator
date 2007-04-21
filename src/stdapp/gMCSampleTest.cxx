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
#include <TText.h>
#include <TStyle.h>

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
void   InitInputEventStream (void);
void   InitOutput           (void);
void   EventLoop            (void);
void   SaveResults          (void);
void   CleanUp              (void);
void   AddEvtFracDir        (string dir, string title, int nupdg);
void   AddKineDir           (string dir, string title, int nupdg, string proc);
void   AddHMultDir          (string dir, string title, int nupdg, string proc);
void   AddVtxDir            (void);
void   AddNtpDir            (void);
void   Plot                 (void);
void   AddFrontPage         (void);
void   PlotEvtFrac          (string dir, string title);
void   PlotKine             (string dir, string title);
void   PlotHMult            (string dir, string title);
void   PlotVtx              (void);
void   PlotH1F              (string name, string title, bool keep_page = false);
void   PlotH2F              (string name, string title);

// command-line arguments
Long64_t gOptN;                // (-n)  process so many events, all if -1
Long64_t gOptTgtPdgC;          // (-t) process events only for this nucl. target
string   gOptInpFile;          // (-f) input GENIE event sample file
string   gOptInpTemplFile;     // (-c) input GENIE event sample file (template)

// internal vars
bool               gSampleComp = false;  // run sample comparisons ?
TCanvas *          gC = 0;
TPostScript *      gPS = 0;
NtpMCEventRecord * gMCRec = 0;
TTree *            gEventRecTree = 0;
Long64_t           gNEvt = 0;
string             gCurrInpFilename = "";
TFile *            gCurrInpFile = 0;
TFile *            gCurrOutFile = 0;
TFile *            gTestedSamplePlotFile = 0;
TFile *            gTempltSamplePlotFile = 0;
TDirectory *       gTestedSampleDir = 0; // plot dir in test sample file
TDirectory *       gTempltSampleDir = 0; // plot dir in template sample file

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

  InitInputEventStream();
  InitOutput();
  EventLoop();
  SaveResults();
  CleanUp();

  LOG("gmctest", pINFO)  << "Done analyzing : " << filename;
}
//_________________________________________________________________________________
void InitInputEventStream(void)
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
}
//_________________________________________________________________________________
void InitOutput(void)
{
  //-- build output filename based on input filename
  string outpfilename = OutputFileName(gCurrInpFilename);
  LOG("gmctest", pNOTICE) 
          << "*** Saving output histograms to: " << outpfilename;

  // open output file
  gCurrOutFile = new TFile(outpfilename.c_str(),"recreate");

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

    LOG("gmctest", pDEBUG) << "....... getting evt " << i;
    gEventRecTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    //LOG("gmctest", pINFO) << rec_header;
    //LOG("gmctest", pINFO) << event;

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
  LOG("gmctest", pNOTICE) << "*** Cleaning up";

  delete gCurrInpFile;
  gCurrInpFile=0;

  gCurrOutFile->Close();
  delete gCurrOutFile;
  gCurrOutFile=0;

  gMCRec=0;
  gEventRecTree=0;
}
//_________________________________________________________________________________
void SaveResults(void)
{
  LOG("gmctest", pNOTICE) 
          << "number of events processed: " << tEvtTree->GetEntries();

  AddEvtFracDir("EvtFracNumuDir", "Event fractions", 14);

  // --- add directories containing kinematics plots
  //
  AddKineDir("KineDir",          "Kinematics - All events",        0, "");
  AddKineDir("KineCcNumuDir",    "Kinematics - All nu_mu CC",     14, "weakcc");
  AddKineDir("KineCcNumuQelDir", "Kinematics - All nu_mu CC QEL", 14, "weakcc&&qel");
  AddKineDir("KineCcNumuResDir", "Kinematics - All nu_mu CC RES", 14, "weakcc&&res");
  AddKineDir("KineCcNumuDisDir", "Kinematics - All nu_mu CC DIS", 14, "weakcc&&dis");
  AddKineDir("KineNcNumuDir",    "Kinematics - All nu_mu NC",     14, "weaknc");
  AddKineDir("KineNcNumuQelDir", "Kinematics - All nu_mu NC QEL", 14, "weaknc&&qel");
  AddKineDir("KineNcNumuResDir", "Kinematics - All nu_mu NC RES", 14, "weaknc&&res");
  AddKineDir("KineNcNumuDisDir", "Kinematics - All nu_mu NC DIS", 14, "weaknc&&dis");

  // --- add directories containing hadronic multiplicity plots
  //
  AddHMultDir("HadMultDir",          "Hadronic multiplicities - All events",        0, "");
  AddHMultDir("HadMultCcNumuDir",    "Hadronic multiplicities - All nu_mu CC",     14, "weakcc");
  AddHMultDir("HadMultCcNumuQelDir", "Hadronic multiplicities - All nu_mu CC QEL", 14, "weakcc&&qel");
  AddHMultDir("HadMultCcNumuResDir", "Hadronic multiplicities - All nu_mu CC RES", 14, "weakcc&&res");
  AddHMultDir("HadMultCcNumuDisDir", "Hadronic multiplicities - All nu_mu CC DIS", 14, "weakcc&&dis");
  AddHMultDir("HadMultNcNumuDir",    "Hadronic multiplicities - All nu_mu NC",     14, "weaknc");
  AddHMultDir("HadMultNcNumuQelDir", "Hadronic multiplicities - All nu_mu NC QEL", 14, "weaknc&&qel");
  AddHMultDir("HadMultNcNumuResDir", "Hadronic multiplicities - All nu_mu NC RES", 14, "weaknc&&res");
  AddHMultDir("HadMultNcNumuDisDir", "Hadronic multiplicities - All nu_mu NC DIS", 14, "weaknc&&dis");

  // --- add directory containing event vertex plots
  //
  AddVtxDir();

  // --- add directory containing summary ntuuples
  //
  AddNtpDir();
}
//_________________________________________________________________________________
void AddKineDir(string dir, string title, int nupdg, string proc)
{
  TDirectory * tdir = gCurrOutFile->mkdir(dir.c_str(), title.c_str());
  tdir->cd();

  ostringstream condition;
  condition << "wgt*(1";
  if(nupdg)          condition << "&&neu==" << nupdg;
  if (proc.size()>0) condition << "&&" << proc;
  condition << ")";

  tEvtTree->Draw( "x>>hx",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "y>>hy",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "t>>ht",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "Q2>>hQ2",   condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "W>>hW",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "v>>hv",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "xs>>hxs",   condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "ys>>hys",   condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "ts>>hts",   condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "Q2s>>hQ2s", condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "Ws>>hWs",   condition.str().c_str(), "GOFF");

  tdir->Write(dir.c_str());
}
//_________________________________________________________________________________
void AddHMultDir(string dir, string title, int nupdg, string proc)
{
  TDirectory * tdir = gCurrOutFile->mkdir(dir.c_str(), title.c_str());
  tdir->cd();

  ostringstream condition;
  condition << "wgt*(1";
  if(nupdg)          condition << "&&neu==" << nupdg;
  if (proc.size()>0) condition << "&&" << proc;
  condition << ")";

  tEvtTree->Draw( "np>>hnp",         condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nn>>hnn",         condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "npip>>hnpip",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "npim>>hnpim",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "ngamma/2>>hnpi0", condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nKp>>hnKp",       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nKm>>hnKm",       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nK0>>hnK0",       condition.str().c_str(), "GOFF");

  tEvtTree->Draw( "np:Ws>>hnpW",         condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nn:Ws>>hnnW",         condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "npip:Ws>>hnpipW",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "npim:Ws>>hnpimW",     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "ngamma/2:Ws>>hnpi0W", condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nKp:Ws>>hnKpW",       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nKm:Ws>>hnKmW",       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "nK0:Ws>>hnK0W",       condition.str().c_str(), "GOFF");

  tdir->Write(dir.c_str());
}
//_________________________________________________________________________________
void AddEvtFracDir(string dir, string title, int nupdg)
{
  TH1F * hE       = new TH1F("hE",      "",80,0,20);
  TH1F * hEcc     = new TH1F("hEcc",    "",80,0,20);
  TH1F * hEccqel  = new TH1F("hEccqel", "",80,0,20);
  TH1F * hEccres  = new TH1F("hEccres", "",80,0,20);
  TH1F * hEccdis  = new TH1F("hEccdis", "",80,0,20);
  TH1F * hEnc     = new TH1F("hEnc",    "",80,0,20);
  TH1F * hEncqel  = new TH1F("hEncqel", "",80,0,20);
  TH1F * hEncres  = new TH1F("hEncres", "",80,0,20);
  TH1F * hEncdis  = new TH1F("hEncdis", "",80,0,20);

  TH1F * hW       = new TH1F("hW",      "",40,0,20);
  TH1F * hWcc     = new TH1F("hWcc",    "",40,0,20);
  TH1F * hWccqel  = new TH1F("hWccqel", "",40,0,20);
  TH1F * hWccres  = new TH1F("hWccres", "",40,0,20);
  TH1F * hWccdis  = new TH1F("hWccdis", "",40,0,20);
  TH1F * hWnc     = new TH1F("hWnc",    "",40,0,20);
  TH1F * hWncqel  = new TH1F("hWncqel", "",40,0,20);
  TH1F * hWncres  = new TH1F("hWncres", "",40,0,20);
  TH1F * hWncdis  = new TH1F("hWncdis", "",40,0,20);

  TH1F * hQ2      = new TH1F("hQ2",     "",80,0,40);
  TH1F * hQ2cc    = new TH1F("hQ2cc",   "",80,0,40);
  TH1F * hQ2ccqel = new TH1F("hQ2ccqel","",80,0,40);
  TH1F * hQ2ccres = new TH1F("hQ2ccres","",80,0,40);
  TH1F * hQ2ccdis = new TH1F("hQ2ccdis","",80,0,40);
  TH1F * hQ2nc    = new TH1F("hQ2nc",   "",80,0,40);
  TH1F * hQ2ncqel = new TH1F("hQ2ncqel","",80,0,40);
  TH1F * hQ2ncres = new TH1F("hQ2ncres","",80,0,40);
  TH1F * hQ2ncdis = new TH1F("hQ2ncdis","",80,0,40);

  tEvtTree->Draw("Ev>>hE",       Form("wgt*(neu==%d)",nupdg),              "GOFF");
  tEvtTree->Draw("Ev>>hEcc",     Form("wgt*(neu==%d&&weakcc)",nupdg),      "GOFF");
  tEvtTree->Draw("Ev>>hEccqel",  Form("wgt*(neu==%d&&weakcc&&qel)",nupdg), "GOFF");
  tEvtTree->Draw("Ev>>hEccres",  Form("wgt*(neu==%d&&weakcc&&res)",nupdg), "GOFF");
  tEvtTree->Draw("Ev>>hEccdis",  Form("wgt*(neu==%d&&weakcc&&dis)",nupdg), "GOFF");
  tEvtTree->Draw("Ev>>hEnc",     Form("wgt*(neu==%d&&weaknc)",nupdg),      "GOFF");
  tEvtTree->Draw("Ev>>hEncqel",  Form("wgt*(neu==%d&&weaknc&&qel)",nupdg), "GOFF");
  tEvtTree->Draw("Ev>>hEncres",  Form("wgt*(neu==%d&&weaknc&&res)",nupdg), "GOFF");
  tEvtTree->Draw("Ev>>hEncdis",  Form("wgt*(neu==%d&&weaknc&&dis)",nupdg), "GOFF");

  tEvtTree->Draw("Ws>>hW",       Form("wgt*(neu==%d)",nupdg),              "GOFF");
  tEvtTree->Draw("Ws>>hWcc",     Form("wgt*(neu==%d&&weakcc)",nupdg),      "GOFF");
  tEvtTree->Draw("Ws>>hWccqel",  Form("wgt*(neu==%d&&weakcc&&qel)",nupdg), "GOFF");
  tEvtTree->Draw("Ws>>hWccres",  Form("wgt*(neu==%d&&weakcc&&res)",nupdg), "GOFF");
  tEvtTree->Draw("Ws>>hWccdis",  Form("wgt*(neu==%d&&weakcc&&dis)",nupdg), "GOFF");
  tEvtTree->Draw("Ws>>hWnc",     Form("wgt*(neu==%d&&weaknc)",nupdg),      "GOFF");
  tEvtTree->Draw("Ws>>hWncqel",  Form("wgt*(neu==%d&&weaknc&&qel)",nupdg), "GOFF");
  tEvtTree->Draw("Ws>>hWncres",  Form("wgt*(neu==%d&&weaknc&&res)",nupdg), "GOFF");
  tEvtTree->Draw("Ws>>hWncdis",  Form("wgt*(neu==%d&&weaknc&&dis)",nupdg), "GOFF");

  tEvtTree->Draw("Q2s>>hQ2",     Form("wgt*(neu==%d)",nupdg),              "GOFF");
  tEvtTree->Draw("Q2s>>hQ2cc",   Form("wgt*(neu==%d&&weakcc)",nupdg),      "GOFF");
  tEvtTree->Draw("Q2s>>hQ2ccqel",Form("wgt*(neu==%d&&weakcc&&qel)",nupdg), "GOFF");
  tEvtTree->Draw("Q2s>>hQ2ccres",Form("wgt*(neu==%d&&weakcc&&res)",nupdg), "GOFF");
  tEvtTree->Draw("Q2s>>hQ2ccdis",Form("wgt*(neu==%d&&weakcc&&dis)",nupdg), "GOFF");
  tEvtTree->Draw("Q2s>>hQ2nc",   Form("wgt*(neu==%d&&weaknc)",nupdg),      "GOFF");
  tEvtTree->Draw("Q2s>>hQ2ncqel",Form("wgt*(neu==%d&&weaknc&&qel)",nupdg), "GOFF");
  tEvtTree->Draw("Q2s>>hQ2ncres",Form("wgt*(neu==%d&&weaknc&&res)",nupdg), "GOFF");
  tEvtTree->Draw("Q2s>>hQ2ncdis",Form("wgt*(neu==%d&&weaknc&&dis)",nupdg), "GOFF");

  hE       -> Sumw2();
  hEcc     -> Sumw2();
  hEccqel  -> Sumw2();
  hEccres  -> Sumw2();
  hEccdis  -> Sumw2();
  hEnc     -> Sumw2();
  hEncqel  -> Sumw2();
  hEncres  -> Sumw2();
  hEncdis  -> Sumw2();
  hW       -> Sumw2();
  hWcc     -> Sumw2();
  hWccqel  -> Sumw2();
  hWccres  -> Sumw2();
  hWccdis  -> Sumw2();
  hWnc     -> Sumw2();
  hWncqel  -> Sumw2();
  hWncres  -> Sumw2();
  hWncdis  -> Sumw2();
  hQ2      -> Sumw2();
  hQ2cc    -> Sumw2();
  hQ2ccqel -> Sumw2();
  hQ2ccres -> Sumw2();
  hQ2ccdis -> Sumw2();
  hQ2nc    -> Sumw2();
  hQ2ncqel -> Sumw2();
  hQ2ncres -> Sumw2();
  hQ2ncdis -> Sumw2();

  hEcc     -> Divide(hE);
  hEccqel  -> Divide(hE);
  hEccres  -> Divide(hE);
  hEccdis  -> Divide(hE);
  hEnc     -> Divide(hE);
  hEncqel  -> Divide(hE);
  hEncres  -> Divide(hE);
  hEncdis  -> Divide(hE);
  hWcc     -> Divide(hW);
  hWccqel  -> Divide(hW);
  hWccres  -> Divide(hW);
  hWccdis  -> Divide(hW);
  hWnc     -> Divide(hW);
  hWncqel  -> Divide(hW);
  hWncres  -> Divide(hW);
  hWncdis  -> Divide(hW);
  hQ2cc    -> Divide(hQ2);
  hQ2ccqel -> Divide(hQ2);
  hQ2ccres -> Divide(hQ2);
  hQ2ccdis -> Divide(hQ2);
  hQ2nc    -> Divide(hQ2);
  hQ2ncqel -> Divide(hQ2);
  hQ2ncres -> Divide(hQ2);
  hQ2ncdis -> Divide(hQ2);

  TDirectory * tdir = gCurrOutFile->mkdir(dir.c_str(), title.c_str());
  tdir->cd();

  tdir -> Add ( hEcc     );
  tdir -> Add ( hEccqel  );
  tdir -> Add ( hEccres  );
  tdir -> Add ( hEccdis  );
  tdir -> Add ( hEnc     );
  tdir -> Add ( hEncqel  );
  tdir -> Add ( hEncres  );
  tdir -> Add ( hEncdis  );
  tdir -> Add ( hWcc     );
  tdir -> Add ( hWccqel  );
  tdir -> Add ( hWccres  );
  tdir -> Add ( hWccdis  );
  tdir -> Add ( hWnc     );
  tdir -> Add ( hWncqel  );
  tdir -> Add ( hWncres  );
  tdir -> Add ( hWncdis  );
  tdir -> Add ( hQ2cc    );
  tdir -> Add ( hQ2ccqel );
  tdir -> Add ( hQ2ccres );
  tdir -> Add ( hQ2ccdis );
  tdir -> Add ( hQ2nc    );
  tdir -> Add ( hQ2ncqel );
  tdir -> Add ( hQ2ncres );
  tdir -> Add ( hQ2ncdis );

  tdir->Write(dir.c_str());
}
//_________________________________________________________________________________
void AddVtxDir(void)
{
  TDirectory * VtxDir = gCurrOutFile->mkdir("VtxDir", "Event vertex plots");
  VtxDir->cd();

  tEvtTree->Draw("vtxx>>hvtxx",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw("vtxy>>hvtxy",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw("vtxz>>hvtxz",       "wgt*(1==1)", "GOFF");
  tEvtTree->Draw("vtxx:vtxy>>hvtxxy", "wgt*(1==1)", "GOFF");

  VtxDir->Write();
}
//_________________________________________________________________________________
void AddNtpDir(void)
{
  TDirectory * NtpDir = gCurrOutFile->mkdir("NtpDir", "Summary Ntuples");
  NtpDir->cd();
  NtpDir->Add(tEvtTree);
  NtpDir->Write();
}
//_________________________________________________________________________________
void Plot(void)
{
  gTestedSamplePlotFile = new TFile(OutputFileName(gOptInpFile).c_str(), "read");
  gTempltSamplePlotFile = (gSampleComp) ?
                 new TFile(OutputFileName(gOptInpTemplFile).c_str(), "read") : 0;

  gC = new TCanvas("gC");
  gC->SetFillColor(0);
  gC->SetBorderMode(0);
  gC->SetBorderMode(0);

  gPS = new TPostScript("out.ps", 112);

  // --- front page
  //
  AddFrontPage();

  // --- event fraction plots
  //
  PlotEvtFrac("EvtFracNumuDir", "Event fractions");

  // --- kinematics plots
  //
  PlotKine("KineDir",             "All events"       );
  PlotKine("KineCcNumuDir",       "#nu_{#mu} CC"     );
  PlotKine("KineCcNumuQelDir",    "#nu_{#mu} CC QEL" );
  PlotKine("KineCcNumuResDir",    "#nu_{#mu} CC RES" );
  PlotKine("KineCcNumuDisDir",    "#nu_{#mu} CC DIS" );
  PlotKine("KineNcNumuDir",       "#nu_{#mu} NC"     );
  PlotKine("KineNcNumuQelDir",    "#nu_{#mu} NC QEL" );
  PlotKine("KineNcNumuResDir",    "#nu_{#mu} NC RES" );
  PlotKine("KineNcNumuDisDir",    "#nu_{#mu} NC DIS" );

  // --- hadronic multiplicity plots
  //
  PlotHMult("HadMultDir",             "All events"       );
  PlotHMult("HadMultCcNumuDir",       "#nu_{#mu} CC"     );
  PlotHMult("HadMultCcNumuQelDir",    "#nu_{#mu} CC QEL" );
  PlotHMult("HadMultCcNumuResDir",    "#nu_{#mu} CC RES" );
  PlotHMult("HadMultCcNumuDisDir",    "#nu_{#mu} CC DIS" );
  PlotHMult("HadMultNcNumuDir",       "#nu_{#mu} NC"     );
  PlotHMult("HadMultNcNumuQelDir",    "#nu_{#mu} NC QEL" );
  PlotHMult("HadMultNcNumuResDir",    "#nu_{#mu} NC RES" );
  PlotHMult("HadMultNcNumuDisDir",    "#nu_{#mu} NC DIS" );

  // --- vertex position plots
  //
  PlotVtx();

  gPS->Close();

  gTestedSamplePlotFile->Close();
  delete gTestedSamplePlotFile;
  if(gTempltSamplePlotFile) {
    gTempltSamplePlotFile->Close();
    delete gTempltSamplePlotFile;
  }
  delete gC;
  delete gPS;

  gTestedSamplePlotFile=0;
  gTempltSamplePlotFile=0;
  gC=0;
  gPS=0;
}
//_________________________________________________________________________________
void AddFrontPage(void)
{
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
}
//_________________________________________________________________________________
void PlotEvtFrac(string dir, string title)
{
  gTestedSampleDir = (TDirectory*) gTestedSamplePlotFile->Get(dir.c_str());
  gTempltSampleDir = (gSampleComp) ? 
                     (TDirectory*) gTempltSamplePlotFile->Get(dir.c_str()) : 0;

  PlotH1F ( "hEcc",   (title + ", in energy bins, CC events").c_str());
  PlotH1F ( "hEccqel",(title + ", in energy bins, CC QEL events").c_str(),true);
  PlotH1F ( "hEccres",(title + ", in energy bins, CC RES events").c_str(),true);
  PlotH1F ( "hEccdis",(title + ", in energy bins, CC DIS events").c_str(),true);
  PlotH1F ( "hEnc",   (title + ", in energy bins, NC events").c_str());
  PlotH1F ( "hEncqel",(title + ", in energy bins, NC QEL events").c_str());
  PlotH1F ( "hEncres",(title + ", in energy bins, NC RES events").c_str());
  PlotH1F ( "hEncdis",(title + ", in energy bins, NC DIS events").c_str());

  PlotH1F ( "hWcc",   (title + ", in W bins, CC events").c_str());
  PlotH1F ( "hWccqel",(title + ", in W bins, CC QEL events").c_str());
  PlotH1F ( "hWccres",(title + ", in W bins, CC RES events").c_str());
  PlotH1F ( "hWccdis",(title + ", in W bins, CC DIS events").c_str());
  PlotH1F ( "hWnc",   (title + ", in W bins, NC events").c_str());
  PlotH1F ( "hWncqel",(title + ", in W bins, NC QEL events").c_str());
  PlotH1F ( "hWncres",(title + ", in W bins, NC RES events").c_str());
  PlotH1F ( "hWncdis",(title + ", in W bins, NC DIS events").c_str());

  PlotH1F ( "hQ2cc",   (title + ", in Q2 bins, CC events").c_str());
  PlotH1F ( "hQ2ccqel",(title + ", in Q2 bins, CC QEL events").c_str());
  PlotH1F ( "hQ2ccres",(title + ", in Q2 bins, CC RES events").c_str());
  PlotH1F ( "hQ2ccdis",(title + ", in Q2 bins, CC DIS events").c_str());
  PlotH1F ( "hQ2nc",   (title + ", in Q2 bins, NC events").c_str());
  PlotH1F ( "hQ2ncqel",(title + ", in Q2 bins, NC QEL events").c_str());
  PlotH1F ( "hQ2ncres",(title + ", in Q2 bins, NC RES events").c_str());
  PlotH1F ( "hQ2ncdis",(title + ", in Q2 bins, NC DIS events").c_str());
}
//_________________________________________________________________________________
void PlotKine(string dir, string title)
{
  gTestedSampleDir = (TDirectory*) gTestedSamplePlotFile->Get(dir.c_str());
  gTempltSampleDir = (gSampleComp) ? 
                     (TDirectory*) gTempltSamplePlotFile->Get(dir.c_str()) : 0;

  PlotH1F ( "hx",   (title + ", x_{comp}").c_str());
  PlotH1F ( "hy",   (title + ", y_{comp}").c_str());
  PlotH1F ( "ht",   (title + ", t_{comp}").c_str());
  PlotH1F ( "hW",   (title + ", W_{comp}").c_str());
  PlotH1F ( "hQ2",  (title + ", Q^{2}_{comp} (GeV^{2})").c_str());
  PlotH1F ( "hv",   (title + ", v_{comp}").c_str());
  PlotH1F ( "hxs",  (title + ", x_{sel}").c_str());
  PlotH1F ( "hys",  (title + ", y_{sel}").c_str());
  PlotH1F ( "hts",  (title + ", t_{sel}").c_str());
  PlotH1F ( "hWs",  (title + ", W_{sel}").c_str());
  PlotH1F ( "hQ2s", (title + ", Q^{2}_{sel} (GeV^{2})").c_str());
}
//_________________________________________________________________________________
void PlotHMult(string dir, string title)
{
  gTestedSampleDir = (TDirectory*) gTestedSamplePlotFile->Get(dir.c_str());
  gTempltSampleDir = (gSampleComp) ? 
                     (TDirectory*) gTempltSamplePlotFile->Get(dir.c_str()) : 0;

  PlotH1F ( "hnp",   (title + ", num of protons").c_str());
  PlotH1F ( "hnn",   (title + ", num of neutrons").c_str());
  PlotH1F ( "hnpip", (title + ", num of #pi^{+}").c_str());
  PlotH1F ( "hnpim", (title + ", num of #pi^{-}").c_str());
  PlotH1F ( "hnpi0", (title + ", num of #pi^{0}").c_str());
  PlotH1F ( "hnKp",  (title + ", num of K^{+}").c_str());
  PlotH1F ( "hnKm",  (title + ", num of K^{-}").c_str());
  PlotH1F ( "hnK0",  (title + ", num of K^{0}+#bar{K^{0}}").c_str());

  PlotH2F ( "hnpW",   (title + ", num of protons vs W (GeV)").c_str());
  PlotH2F ( "hnnW",   (title + ", num of neutrons vs W (GeV)").c_str());
  PlotH2F ( "hnpipW", (title + ", num of #pi^{+} vs W (GeV)").c_str());
  PlotH2F ( "hnpimW", (title + ", num of #pi^{-} vs W (GeV)").c_str());
  PlotH2F ( "hnpi0W", (title + ", num of #pi^{0} vs W (GeV)").c_str());
  PlotH2F ( "hnKpW",  (title + ", num of K^{+} vs W (GeV)").c_str());
  PlotH2F ( "hnKmW",  (title + ", num of K^{-} vs W (GeV)").c_str());
  PlotH2F ( "hnK0W",  (title + ", num of K^{0}+#bar{K^{0}} vs W (GeV)").c_str());
}
//_________________________________________________________________________________
void PlotVtx(void)
{
  gTestedSampleDir = (TDirectory*) gTestedSamplePlotFile->Get("VtxDir");
  gTempltSampleDir = (gSampleComp) ? 
                     (TDirectory*) gTempltSamplePlotFile->Get("VtxDir") : 0;

  PlotH1F ( "vtxx",  "vertex x");
  PlotH1F ( "vtxy",  "vertex y");
  PlotH1F ( "vtxz",  "vertex z");
}
//_________________________________________________________________________________
void PlotH1F(string name, string title, bool keep_page)
{
  if(!keep_page) gPS->NewPage();
  gC->cd();

  TH1F * tested_hst = dynamic_cast<TH1F *> (gTestedSampleDir->Get(name.c_str()));
  TH1F * templt_hst = (gSampleComp) ?
                      dynamic_cast<TH1F *> (gTempltSampleDir->Get(name.c_str())) : 0;

  // plot histogram from test sample
  if(!tested_hst) return;
  tested_hst->SetLineColor(2);
  tested_hst->SetLineWidth(2);
  tested_hst->Draw();

  // plot same hisogram from template sample (if any) 
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
void PlotH2F(string name, string title)
{
  gPS->NewPage();
  gC->cd();

  TH2F * tested_hst = dynamic_cast<TH2F *> (gTestedSampleDir->Get(name.c_str()));
  TH2F * templt_hst = (gSampleComp) ?
                      dynamic_cast<TH2F *> (gTempltSampleDir->Get(name.c_str())) : 0;

  // plot histogram from test sample
  if(!tested_hst) return;
  //gStyle->SetPalette(1);
  //tested_hst->Draw("COLZ");
  tested_hst->SetFillColor(3);
  tested_hst->Draw();

  TProfile * hpx = tested_hst->ProfileX();
  hpx->SetLineColor(2);
  hpx->SetLineWidth(2);
  hpx->SetMarkerStyle(20);
  hpx->SetMarkerSize(1.3);
  hpx->Draw("SAME");

  // plot same hisogram from template sample (if any)
  if(templt_hst) {
    TProfile * hpxt = templt_hst->ProfileX();
    hpxt->SetLineColor(1);
    hpxt->SetLineWidth(2);
    hpxt->SetMarkerStyle(8);
    hpxt->SetMarkerSize(1.3);
    hpxt->Draw("SAME");
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

