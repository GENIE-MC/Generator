//____________________________________________________________________________
/*!

\program gmctest

\brief   A GENIE utility that reads a generated GENIE event tree and generates 
         a ROOT file with loads of histograms with basic MC truth quantities.
         It also generates the 'definite' summary ROOT ntuple.
         Typically used in validating generator changes.

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
void   Init                 (void);
void   EventLoop            (void);
void   SaveResults          (void);
void   CleanUp              (void);
void   AddEvtFracDir        (string dir, string title, int nupdg);  
void   AddKineDir           (string dir, string title, int nupdg, string proc); 
void   AddHMultDir          (string dir, string title, int nupdg, string proc); 
void   AddHP4Dir            (string dir, string title, int nupdg, string proc); 
void   AddVtxDir            (void); 
void   AddNtpDir            (void);
void   Plot                 (void);
void   AddFrontPage         (void);
void   PlotEvtFrac          (string dir, string title); 
void   PlotKine             (string dir, string title); 
void   PlotHMult            (string dir, string title); 
void   PlotHP4              (string dir, string title); 
void   PlotVtx              (void);                     
void   PlotH1F              (string name, string title, bool keep_page = false);
void   PlotH1F_2            (string name1, string name2, string title, bool keep_page = false);
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
double             gEmin  = 999, gEmax  = -999; // Ev range of events in input sample
double             gWmin  = 999, gWmax  = -999; // W  range of events in input sample
double             gQ2min = 999, gQ2max = -999; // Q2 range of events in input sample

// various control constants
const int    kNPmax      = 100;   //
const double kEBinSz     = 0.250; //
const double kWBinSz     = 0.200; //
const double kQ2BinSz    = 0.100; //

// summary tree
TTree * tEvtTree;
int    br_iev       = 0;      //
int    br_neutrino  = 0;      //
int    br_target    = 0;      //
int    br_hitnuc    = 0;      //
int    br_hitqrk    = 0;      //
bool   br_qel       = false;  //
bool   br_res       = false;  //
bool   br_dis       = false;  //
bool   br_coh       = false;  //
bool   br_em        = false;  //
bool   br_weakcc    = false;  //
bool   br_weaknc    = false;  //
double br_kine_xs   = 0;      //
double br_kine_ys   = 0;      //
double br_kine_ts   = 0;      //
double br_kine_Q2s  = 0;      //
double br_kine_Ws   = 0;      //
double br_kine_x    = 0;      //
double br_kine_y    = 0;      //
double br_kine_t    = 0;      //
double br_kine_Q2   = 0;      //
double br_kine_W    = 0;      //
double br_kine_v    = 0;      //
double br_Ev        = 0;      //
double br_El        = 0;      //
double br_vtxx      = 0;      // vertex position x 
double br_vtxy      = 0;      // vertex position y
double br_vtxz      = 0;      // vertex position z
double br_weight    = 0;      // event weight
int    br_np        = 0;      // number of f/s p
int    br_nn        = 0;      // number of f/s n
int    br_npip      = 0;      // number of f/s pi+
int    br_npim      = 0;      // number of f/s pi-
int    br_npi0      = 0;      // number of f/s pi0
int    br_ngamma    = 0;      // number of f/s gamma
int    br_nKp       = 0;      // number of f/s K+
int    br_nKm       = 0;      // number of f/s K-
int    br_nK0       = 0;      // number of f/s K0
int    br_hmod      = 0;      // hadronization model (-1: nodis, 0:kno, 1:string, 2:cluster, 3:indep)
int    br_nhep      = 0;      // number of GHEP record entries
int    br_pdg[kNPmax];        //
int    br_ist[kNPmax];        //
double br_px [kNPmax];        //
double br_py [kNPmax];        //
double br_pz [kNPmax];        //
double br_p  [kNPmax];        //
double br_E  [kNPmax];        //
double br_KE [kNPmax];        //
double br_x  [kNPmax];        //
double br_y  [kNPmax];        //
double br_z  [kNPmax];        //
int    br_da1[kNPmax];        //
int    br_da2[kNPmax];        //
int    br_mom[kNPmax];        //
int    br_n_hprim = 0;        // number of particles at the prim hadronic system (before transport)
int    br_pdg_hprim[kNPmax];  //
double br_px_hprim [kNPmax];  //
double br_py_hprim [kNPmax];  //
double br_pz_hprim [kNPmax];  //
double br_p_hprim  [kNPmax];  //
double br_E_hprim  [kNPmax];  //
double br_KE_hprim [kNPmax];  //
double br_had_E;              // hadronic energy 
double br_had_px;             // -//-     px      
double br_had_py;             // -//-     py      
double br_had_pz;             // -//-     pz      
double br_vishad_E;           // 'visible' hadronic energy 
double br_vishad_px;          // -//-     px      
double br_vishad_py;          // -//-     py      
double br_vishad_pz;          // -//-     pz      
double br_mishad_E;           // missing hadronic energy (because of intranuclear rescattering)
double br_mishad_px;          //      -//-        px      -//-
double br_mishad_py;          //      -//-        py      -//-
double br_mishad_pz;          //      -//-        pz      -//-

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
    gSampleComp = true;
    AnalyzeSample(gOptInpTemplFile);
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

  Init();
  EventLoop();
  SaveResults();
  CleanUp();

  LOG("gmctest", pINFO)  << "Done analyzing : " << filename;
}
//_________________________________________________________________________________
void Init(void)
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

  //-- build output filename based on input filename
  string outpfilename = OutputFileName(gCurrInpFilename);
  LOG("gmctest", pNOTICE) 
          << "*** Saving output histograms to: " << outpfilename;

  // open output file
  gCurrOutFile = new TFile(outpfilename.c_str(),"recreate");

  //-- create summary tree
  tEvtTree = new TTree("tEvtTree","event tree summary");
  tEvtTree->Branch("iev",       &br_iev,        "iev/I"       );
  tEvtTree->Branch("neu",       &br_neutrino ,  "neu/I"       );
  tEvtTree->Branch("tgt" ,      &br_target,     "tgt/I"       );
  tEvtTree->Branch("hitnuc",    &br_hitnuc,     "hitnuc/I"    );
  tEvtTree->Branch("hitqrk",    &br_hitqrk,     "hitqrk/I"    );
  tEvtTree->Branch("qel",       &br_qel,        "br_qel/O"    );
  tEvtTree->Branch("res",       &br_res,        "br_res/O"    );
  tEvtTree->Branch("dis",       &br_dis,        "br_dis/O"    );
  tEvtTree->Branch("coh",       &br_coh,        "br_coh/O"    );
  tEvtTree->Branch("em",        &br_em,         "br_em/O"     );
  tEvtTree->Branch("weakcc",    &br_weakcc,     "br_weakcc/O" );
  tEvtTree->Branch("weaknc",    &br_weaknc,     "br_weaknc/O" );
  tEvtTree->Branch("xs",        &br_kine_xs,    "xs/D"        );
  tEvtTree->Branch("ys",        &br_kine_ys,    "ys/D"        );
  tEvtTree->Branch("ts",        &br_kine_ts,    "ts/D"        );
  tEvtTree->Branch("Q2s",       &br_kine_Q2s ,  "Q2s/D"       );
  tEvtTree->Branch("Ws",        &br_kine_Ws,    "Ws/D"        );
  tEvtTree->Branch("x",         &br_kine_x,     "x/D"         );
  tEvtTree->Branch("y",         &br_kine_y,     "y/D"         );
  tEvtTree->Branch("t",         &br_kine_t,     "t/D"         );
  tEvtTree->Branch("Q2",        &br_kine_Q2,    "Q2/D"        );
  tEvtTree->Branch("W",         &br_kine_W,     "W/D"         );
  tEvtTree->Branch("v",         &br_kine_v,     "v/D"         );
  tEvtTree->Branch("Ev",        &br_Ev,         "Ev/D"        );
  tEvtTree->Branch("El",        &br_El,         "El/D"        );
  tEvtTree->Branch("vtxx",      &br_vtxx,       "vtxx/D"      );
  tEvtTree->Branch("vtxy",      &br_vtxy,       "vtxy/D"      );
  tEvtTree->Branch("vtxz",      &br_vtxz,       "vtxz/D"      );
  tEvtTree->Branch("wgt",       &br_weight,     "wgt/D"       );
  tEvtTree->Branch("np",        &br_np,         "np/I"        );
  tEvtTree->Branch("nn",        &br_nn,         "nn/I"        );
  tEvtTree->Branch("npip",      &br_npip,       "npip/I"      );
  tEvtTree->Branch("npim",      &br_npim,       "npim/I"      );
  tEvtTree->Branch("npi0",      &br_npi0,       "npi0/I"      );
  tEvtTree->Branch("ngamma",    &br_ngamma,     "ngamma/I"    );
  tEvtTree->Branch("nKp",       &br_nKp,        "nKp/I"       );
  tEvtTree->Branch("nKm",       &br_nKm,        "nKm/I"       );
  tEvtTree->Branch("nK0",       &br_nK0,        "nK0/I"       );
  tEvtTree->Branch("hmod",      &br_hmod,       "hmod/I"      );
  tEvtTree->Branch("nhep",      &br_nhep,       "nhep/I"      );
  tEvtTree->Branch("pdg",        br_pdg,        "pdg[nhep]/I" );
  tEvtTree->Branch("ist",        br_ist,        "ist[nhep]/I" );
  tEvtTree->Branch("px",         br_px,         "px[nhep]/D"  );
  tEvtTree->Branch("py",         br_py,         "py[nhep]/D"  );
  tEvtTree->Branch("pz",         br_pz,         "pz[nhep]/D"  );
  tEvtTree->Branch("p" ,         br_p,          "p[nhep]/D"   );
  tEvtTree->Branch("E",          br_E,          "E[nhep]/D"   );
  tEvtTree->Branch("KE",         br_KE,         "KE[nhep]/D"  );
  tEvtTree->Branch("x",          br_x,          "x[nhep]/D"   );
  tEvtTree->Branch("y",          br_y,          "y[nhep]/D"   );
  tEvtTree->Branch("z",          br_z,          "z[nhep]/D"   );
  tEvtTree->Branch("fdaughter",  br_da1,        "fdaughter[nhep]/I" );
  tEvtTree->Branch("ldaughter",  br_da2,        "ldaughter[nhep]/I" );
  tEvtTree->Branch("mom",        br_mom,        "mom[nhep]/I" );
  tEvtTree->Branch("n_hprim",   &br_n_hprim,    "n_hprim/I"            );
  tEvtTree->Branch("pdg_hprim",  br_pdg_hprim,  "pdg_hprim[n_hprim]/I" );
  tEvtTree->Branch("px_hprim",   br_px_hprim,   "px_hprim[n_hprim]/D"  );
  tEvtTree->Branch("py_hprim",   br_py_hprim,   "py_hprim[n_hprim]/D"  );
  tEvtTree->Branch("pz_hprim",   br_pz_hprim,   "pz_hprim[n_hprim]/D"  );
  tEvtTree->Branch("p_hprim" ,   br_p_hprim,    "p_hprim[n_hprim]/D"   );
  tEvtTree->Branch("E_hprim",    br_E_hprim,    "E_hprim[n_hprim]/D"   );
  tEvtTree->Branch("KE_hprim",   br_KE_hprim,   "KE_hprim[n_hprim]/D"  );
  tEvtTree->Branch("had_E",     &br_had_E,      "had_E/D"              );
  tEvtTree->Branch("had_px",    &br_had_px,     "had_px/D"             );
  tEvtTree->Branch("had_py",    &br_had_py,     "had_py/D"             );
  tEvtTree->Branch("had_pz",    &br_had_pz,     "had_pz/D"             );
  tEvtTree->Branch("vishad_E",  &br_vishad_E,   "vishad_E/D"           );
  tEvtTree->Branch("vishad_px", &br_vishad_px,  "vishad_px/D"          );
  tEvtTree->Branch("vishad_py", &br_vishad_py,  "vishad_py/D"          );
  tEvtTree->Branch("vishad_pz", &br_vishad_pz,  "vishad_pz/D"          );
  tEvtTree->Branch("mishad_E",  &br_mishad_E,   "mishad_E/D"           );
  tEvtTree->Branch("mishad_px", &br_mishad_px,  "mishad_px/D"          );
  tEvtTree->Branch("mishad_py", &br_mishad_py,  "mishad_py/D"          );
  tEvtTree->Branch("mishad_pz", &br_mishad_pz,  "mishad_pz/D"          );

  br_nhep    = 0; 
  br_n_hprim = 0; 
  for(int k=0; k<kNPmax; k++) {
     br_pdg      [k] =  0;        
     br_ist      [k] = -1;        
     br_px       [k] =  0;        
     br_py       [k] =  0;        
     br_pz       [k] =  0;        
     br_p        [k] =  0;        
     br_E        [k] =  0;        
     br_KE       [k] =  0;        
     br_x        [k] =  0;        
     br_y        [k] =  0;        
     br_z        [k] =  0;        
     br_da1      [k] = -1;        
     br_da2      [k] = -1;        
     br_mom      [k] = -1;        
     br_pdg_hprim[k] =  0;  
     br_px_hprim [k] =  0;  
     br_py_hprim [k] =  0;  
     br_pz_hprim [k] =  0;  
     br_p_hprim  [k] =  0;  
     br_E_hprim  [k] =  0;  
     br_KE_hprim [k] =  0;  
  }
}
//_________________________________________________________________________________
void EventLoop(void)
{
  LOG("gmctest", pNOTICE) << "*** Analyzing: " << gNEvt << " events";

  if (gNEvt<0)        return;
  if (!gEventRecTree) return;
  if (!gMCRec)        return;

  for(Long64_t i = 0; i < gNEvt; i++) {

    gEventRecTree->GetEntry(i);

    NtpMCRecHeader rec_header = gMCRec->hdr;
    EventRecord &  event      = *(gMCRec->event);

    //LOG("gmctest", pINFO) << rec_header;
    //LOG("gmctest", pINFO) << event;

    // go further only if the event is physical
    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) {
      gMCRec->Clear();
      continue;
    }

    // nuclear target or free nucleon target
    GHepParticle * target = event.Particle(1);
    assert(target);

    // only further only if it matches the requested target
    if(target->Pdg() != gOptTgtPdgC) {
      gMCRec->Clear();
      continue;
    }

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

    bool is_qel    = proc_info.IsQuasiElastic();
    bool is_res    = proc_info.IsResonant();
    bool is_dis    = proc_info.IsDeepInelastic();
    bool is_coh    = proc_info.IsCoherent();
    bool is_em     = proc_info.IsEM();
    bool is_weakcc = proc_info.IsWeakCC();
    bool is_weaknc = proc_info.IsWeakNC();

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

    // f/s multiplicities
    //
    int np      = event.NEntries(kPdgProton,  kIStStableFinalState);
    int nn      = event.NEntries(kPdgNeutron, kIStStableFinalState);
    int npip    = event.NEntries(kPdgPiP,     kIStStableFinalState);
    int npim    = event.NEntries(kPdgPiM,     kIStStableFinalState);
    int npi0    = event.NEntries(kPdgPi0,     kIStStableFinalState);
    int ngamma  = event.NEntries(kPdgGamma,   kIStStableFinalState);
    int nKp     = event.NEntries(kPdgKP,      kIStStableFinalState);
    int nKm     = event.NEntries(kPdgKM,      kIStStableFinalState);
    int nK0     = event.NEntries(kPdgK0,      kIStStableFinalState);

    // update ranges
    //
    if(!gSampleComp) {
      gEmin  = TMath::Min( gEmin,  Ev    );
      gEmax  = TMath::Max( gEmax,  Ev    );
      gWmin  = TMath::Min( gWmin,  W     );
      gWmax  = TMath::Max( gWmax,  W     );
      gQ2min = TMath::Min( gQ2min, Q2    );
      gQ2max = TMath::Max( gQ2max, Q2    );
    }

    // init hadronization model flag
    if(is_dis) br_hmod =  0;  
    else       br_hmod = -1;

    // init root position for hadronic system (before intranuclear transport)
    int ihad_root = -1;
    if(!is_dis) ihad_root = event.HitNucleonPosition();

    // copy GHEP record into a flat array
    //
    TObjArrayIter piter(&event);
    GHepParticle * p = 0;
    unsigned int ip=0;
    while( (p = (GHepParticle *) piter.Next()) )
    {
      br_pdg[ip] = p->Pdg(); 
      br_ist[ip] = p->Status(); 
      br_px [ip] = p->Px(); 
      br_py [ip] = p->Py(); 
      br_pz [ip] = p->Pz(); 
      br_p  [ip] = p->P4()->Vect().Mag(); 
      br_E  [ip] = p->Energy(); 
      br_KE [ip] = p->KinE(); 
      br_x  [ip] = p->Vx(); 
      br_y  [ip] = p->Vy(); 
      br_z  [ip] = p->Vz(); 
      br_da1[ip] = p->FirstDaughter();
      br_da2[ip] = p->LastDaughter();
      br_mom[ip] = p->FirstMother();

      if(is_dis) {
         if(p->Pdg() == kPdgHadronicSyst) { ihad_root = ip; }

         if      (p->Pdg() == kPdgString ) { br_hmod=1; ihad_root = ip; }
         else if (p->Pdg() == kPdgCluster) { br_hmod=2; ihad_root = ip; }
         else if (p->Pdg() == kPdgIndep  ) { br_hmod=3; ihad_root = ip; }
      }

      ip++;
    }
    br_nhep = event.GetEntries();

    if(ihad_root>0) {
      int ih=0;
      GHepParticle * hadroot = event.Particle(ihad_root);    
      for(int j=hadroot->FirstDaughter(); j<hadroot->LastDaughter(); j++) {
        GHepParticle * hj = event.Particle(j);    
        if(hj->Status() == kIStHadronInTheNucleus) {
           br_pdg_hprim[ih] = hj->Pdg();  
           br_px_hprim [ih] = hj->Px();  
           br_py_hprim [ih] = hj->Py();  
           br_pz_hprim [ih] = hj->Pz();  
           br_p_hprim  [ih] = hj->P4()->Vect().Mag();  
           br_E_hprim  [ih] = hj->E();  
           br_KE_hprim [ih] = hj->KinE();  
           ih++;
        } else if(hj->Status() == kIStDecayedState) {
           for(int k=hj->FirstDaughter(); k<hadroot->LastDaughter(); k++) {
              GHepParticle * hk = event.Particle(k);    
              if(hk->Status() == kIStHadronInTheNucleus) {
                 br_pdg_hprim[ih] = hk->Pdg();  
                 br_px_hprim [ih] = hk->Px();  
                 br_py_hprim [ih] = hk->Py();  
                 br_pz_hprim [ih] = hk->Pz();  
                 br_p_hprim  [ih] = hk->P4()->Vect().Mag();  
                 br_E_hprim  [ih] = hk->E();  
                 br_KE_hprim [ih] = hk->KinE();  
                 ih++;
              } // hadronic decay products
          }//decay products
        } // has decayed before transport
      } // j
      br_n_hprim = ih; 
    } else {
      br_n_hprim = 0; 
    }

    br_had_E     = 0.;
    br_had_px    = 0.;
    br_had_py    = 0.;
    br_had_pz    = 0.;
    br_vishad_E  = 0.;
    br_vishad_px = 0.;
    br_vishad_py = 0.;
    br_vishad_pz = 0.;
    br_mishad_E  = 0.;
    br_mishad_px = 0.;
    br_mishad_py = 0.;
    br_mishad_pz = 0.;
    if(ihad_root) {
	br_had_E  = event.Particle(ihad_root)->Energy();
	br_had_px = event.Particle(ihad_root)->Px();
	br_had_py = event.Particle(ihad_root)->Py();
	br_had_pz = event.Particle(ihad_root)->Pz();
    }


    br_iev       = i;
    br_neutrino  = neutrino->Pdg();
    br_target    = 0;
    br_hitnuc    = 0;
    br_hitqrk    = 0;
    br_qel       = is_qel;
    br_res       = is_res;
    br_dis       = is_dis;
    br_coh       = is_coh;
    br_em        = is_em;
    br_weakcc    = is_weakcc;
    br_weaknc    = is_weaknc;
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
    br_np        = np;
    br_nn        = nn;
    br_npip      = npip;
    br_npim      = npim;
    br_npi0      = npi0;
    br_ngamma    = ngamma;
    br_nKp       = nKp;
    br_nKm       = nKm;
    br_nK0       = nK0;

    tEvtTree->Fill();

    gMCRec->Clear();

  } // event loop

//  assert(gEmin < gEmax);
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

  // -- add directory containing event fraction plots as a function of kinematical 
  //    quantities
  //
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

  // --- add directories containing hadronic 4-momentum plots
  //
  AddHP4Dir("HadP4Dir",          "Hadronic 4-momentum - All events",        0, "");
  AddHP4Dir("HadP4CcNumuDir",    "Hadronic 4-momentum - All nu_mu CC",     14, "weakcc");
  AddHP4Dir("HadP4CcNumuQelDir", "Hadronic 4-momentum - All nu_mu CC QEL", 14, "weakcc&&qel");
  AddHP4Dir("HadP4CcNumuResDir", "Hadronic 4-momentum - All nu_mu CC RES", 14, "weakcc&&res");
  AddHP4Dir("HadP4CcNumuDisDir", "Hadronic 4-momentum - All nu_mu CC DIS", 14, "weakcc&&dis");
  AddHP4Dir("HadP4NcNumuDir",    "Hadronic 4-momentum - All nu_mu NC",     14, "weaknc");
  AddHP4Dir("HadP4NcNumuQelDir", "Hadronic 4-momentum - All nu_mu NC QEL", 14, "weaknc&&qel");
  AddHP4Dir("HadP4NcNumuResDir", "Hadronic 4-momentum - All nu_mu NC RES", 14, "weaknc&&res");
  AddHP4Dir("HadP4NcNumuDisDir", "Hadronic 4-momentum - All nu_mu NC DIS", 14, "weaknc&&dis");

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

  tEvtTree->Draw( "x>>hx",                       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "y>>hy",                       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "t>>ht",                       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "log10(Q2)>>hlog10Q2",         condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "W>>hW",                       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "v>>hv",                       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "xs>>hxs",                     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "ys>>hys",                     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "ts>>hts",                     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "log10(Q2s)>>hlog10Q2s",       condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "Ws>>hWs",                     condition.str().c_str(), "GOFF");
  tEvtTree->Draw( "log10(Q2s):Ws>>hlog10Q2Ws",   condition.str().c_str(), "GOFF");

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
void AddHP4Dir(string dir, string title, int nupdg, string proc)
{
  TDirectory * tdir = gCurrOutFile->mkdir(dir.c_str(), title.c_str());
  tdir->cd();

  ostringstream condition;
  condition << "wgt*(1";
  if(nupdg)          condition << "&&neu==" << nupdg;
  if (proc.size()>0) condition << "&&" << proc;

  tEvtTree->Draw( "px>>hpx_pip",  (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "py>>hpy_pip",  (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "pz>>hpz_pip",  (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "p>>hp_pip",    (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E>>hE_pip",    (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "KE>>hKE_pip",  (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");

  tEvtTree->Draw( "px>>hpx_pim",  (condition.str() + "&&ist==1&&pdg==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "py>>hpy_pim",  (condition.str() + "&&ist==1&&pdg==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "pz>>hpz_pim",  (condition.str() + "&&ist==1&&pdg==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "p>>hp_pim",    (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E>>hE_pim",    (condition.str() + "&&ist==1&&pdg==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "KE>>hKE_pim",  (condition.str() + "&&ist==1&&pdg==-211)").c_str(), "GOFF");

  tEvtTree->Draw( "px>>hpx_p",    (condition.str() + "&&ist==1&&pdg==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "py>>hpy_p",    (condition.str() + "&&ist==1&&pdg==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "pz>>hpz_p",    (condition.str() + "&&ist==1&&pdg==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "p>>hp_p",      (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E>>hE_p",      (condition.str() + "&&ist==1&&pdg==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "KE>>hKE_p",    (condition.str() + "&&ist==1&&pdg==2212)").c_str(), "GOFF");

  tEvtTree->Draw( "px>>hpx_n",    (condition.str() + "&&ist==1&&pdg==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "py>>hpy_n",    (condition.str() + "&&ist==1&&pdg==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "pz>>hpz_n",    (condition.str() + "&&ist==1&&pdg==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "p>>hp_n",      (condition.str() + "&&ist==1&&pdg==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E>>hE_n",      (condition.str() + "&&ist==1&&pdg==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "KE>>hKE_n",    (condition.str() + "&&ist==1&&pdg==2112)").c_str(), "GOFF");


  // add plots before hadron transport (primary hadronic system)

  tEvtTree->Draw( "px_hprim>>hpx_pip_hprim",  (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "py_hprim>>hpy_pip_hprim",  (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "pz_hprim>>hpz_pip_hprim",  (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "p_hprim>>hp_pip_hprim",    (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E_hprim>>hE_pip_hprim",    (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "KE_hprim>>hKE_pip_hprim",  (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");

  tEvtTree->Draw( "px_hprim>>hpx_pim_hprim",  (condition.str() + "&&pdg_hprim==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "py_hprim>>hpy_pim_hprim",  (condition.str() + "&&pdg_hprim==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "pz_hprim>>hpz_pim_hprim",  (condition.str() + "&&pdg_hprim==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "p_hprim>>hp_pim_hprim",    (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E_hprim>>hE_pim_hprim",    (condition.str() + "&&pdg_hprim==-211)").c_str(), "GOFF");
  tEvtTree->Draw( "KE_hprim>>hKE_pim_hprim",  (condition.str() + "&&pdg_hprim==-211)").c_str(), "GOFF");

  tEvtTree->Draw( "px_hprim>>hpx_hp_prim",    (condition.str() + "&&pdg_hprim==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "py_hprim>>hpy_hp_prim",    (condition.str() + "&&pdg_hprim==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "pz_hprim>>hpz_hp_prim",    (condition.str() + "&&pdg_hprim==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "p_hprim>>hp_p_hprim",      (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E_hprim>>hE_p_hprim",      (condition.str() + "&&pdg_hprim==2212)").c_str(), "GOFF");
  tEvtTree->Draw( "KE_hprim>>hKE_p_hprim",    (condition.str() + "&&pdg_hprim==2212)").c_str(), "GOFF");

  tEvtTree->Draw( "px_hprim>>hpx_n_hprim",    (condition.str() + "&&pdg_hprim==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "py_hprim>>hpy_n_hprim",    (condition.str() + "&&pdg_hprim==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "pz_hprim>>hpz_n_hprim",    (condition.str() + "&&pdg_hprim==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "p_hprim>>hp_n_hprim",      (condition.str() + "&&pdg_hprim==211)").c_str(),  "GOFF");
  tEvtTree->Draw( "E_hprim>>hE_n_hprim",      (condition.str() + "&&pdg_hprim==2112)").c_str(), "GOFF");
  tEvtTree->Draw( "KE_hprim>>hKE_n_hprim",    (condition.str() + "&&pdg_hprim==2112)").c_str(), "GOFF");

  tdir->Write(dir.c_str());
}
//_________________________________________________________________________________
void AddEvtFracDir(string dir, string title, int nupdg)
{
  double dE      = gEmax  - gEmin;
  double dW      = gWmax  - gWmin;
  double dQ2     = gQ2max - gQ2min;
  int    nbe     = dE     / kEBinSz;
  int    nbw     = dW     / kWBinSz;
  int    nbq2    = dQ2   / kQ2BinSz;

  TH1F * hE       = new TH1F("hE",      "", nbe, gEmin, gEmax);
  TH1F * hEcc     = new TH1F("hEcc",    "", nbe, gEmin, gEmax);
  TH1F * hEccqel  = new TH1F("hEccqel", "", nbe, gEmin, gEmax);
  TH1F * hEccres  = new TH1F("hEccres", "", nbe, gEmin, gEmax);
  TH1F * hEccdis  = new TH1F("hEccdis", "", nbe, gEmin, gEmax);
  TH1F * hEnc     = new TH1F("hEnc",    "", nbe, gEmin, gEmax);
  TH1F * hEncqel  = new TH1F("hEncqel", "", nbe, gEmin, gEmax);
  TH1F * hEncres  = new TH1F("hEncres", "", nbe, gEmin, gEmax);
  TH1F * hEncdis  = new TH1F("hEncdis", "", nbe, gEmin, gEmax);

  TH1F * hW       = new TH1F("hW",      "", nbw, gWmin, gWmax);
  TH1F * hWcc     = new TH1F("hWcc",    "", nbw, gWmin, gWmax);
  TH1F * hWccqel  = new TH1F("hWccqel", "", nbw, gWmin, gWmax);
  TH1F * hWccres  = new TH1F("hWccres", "", nbw, gWmin, gWmax);
  TH1F * hWccdis  = new TH1F("hWccdis", "", nbw, gWmin, gWmax);
  TH1F * hWnc     = new TH1F("hWnc",    "", nbw, gWmin, gWmax);
  TH1F * hWncqel  = new TH1F("hWncqel", "", nbw, gWmin, gWmax);
  TH1F * hWncres  = new TH1F("hWncres", "", nbw, gWmin, gWmax);
  TH1F * hWncdis  = new TH1F("hWncdis", "", nbw, gWmin, gWmax);

  TH1F * hQ2      = new TH1F("hQ2",     "", nbq2, gQ2min, gQ2max);
  TH1F * hQ2cc    = new TH1F("hQ2cc",   "", nbq2, gQ2min, gQ2max);
  TH1F * hQ2ccqel = new TH1F("hQ2ccqel","", nbq2, gQ2min, gQ2max);
  TH1F * hQ2ccres = new TH1F("hQ2ccres","", nbq2, gQ2min, gQ2max);
  TH1F * hQ2ccdis = new TH1F("hQ2ccdis","", nbq2, gQ2min, gQ2max);
  TH1F * hQ2nc    = new TH1F("hQ2nc",   "", nbq2, gQ2min, gQ2max);
  TH1F * hQ2ncqel = new TH1F("hQ2ncqel","", nbq2, gQ2min, gQ2max);
  TH1F * hQ2ncres = new TH1F("hQ2ncres","", nbq2, gQ2min, gQ2max);
  TH1F * hQ2ncdis = new TH1F("hQ2ncdis","", nbq2, gQ2min, gQ2max);

  hE       -> SetMaximum(1.2);
  hEcc     -> SetMaximum(1.2);
  hEccqel  -> SetMaximum(1.2);
  hEccres  -> SetMaximum(1.2);
  hEccdis  -> SetMaximum(1.2);
  hEnc     -> SetMaximum(1.2);
  hEncqel  -> SetMaximum(1.2);
  hEncres  -> SetMaximum(1.2);
  hEncdis  -> SetMaximum(1.2);

  hW       -> SetMaximum(1.2);
  hWcc     -> SetMaximum(1.2);
  hWccqel  -> SetMaximum(1.2);
  hWccres  -> SetMaximum(1.2);
  hWccdis  -> SetMaximum(1.2);
  hWnc     -> SetMaximum(1.2);
  hWncqel  -> SetMaximum(1.2);
  hWncres  -> SetMaximum(1.2);
  hWncdis  -> SetMaximum(1.2);

  hQ2      -> SetMaximum(1.2);
  hQ2cc    -> SetMaximum(1.2);
  hQ2ccqel -> SetMaximum(1.2);
  hQ2ccres -> SetMaximum(1.2);
  hQ2ccdis -> SetMaximum(1.2);
  hQ2nc    -> SetMaximum(1.2);
  hQ2ncqel -> SetMaximum(1.2);
  hQ2ncres -> SetMaximum(1.2);
  hQ2ncdis -> SetMaximum(1.2);

  hE       -> SetMinimum(0.0);
  hEcc     -> SetMinimum(0.0);
  hEccqel  -> SetMinimum(0.0);
  hEccres  -> SetMinimum(0.0);
  hEccdis  -> SetMinimum(0.0);
  hEnc     -> SetMinimum(0.0);
  hEncqel  -> SetMinimum(0.0);
  hEncres  -> SetMinimum(0.0);
  hEncdis  -> SetMinimum(0.0);

  hW       -> SetMinimum(0.0);
  hWcc     -> SetMinimum(0.0);
  hWccqel  -> SetMinimum(0.0);
  hWccres  -> SetMinimum(0.0);
  hWccdis  -> SetMinimum(0.0);
  hWnc     -> SetMinimum(0.0);
  hWncqel  -> SetMinimum(0.0);
  hWncres  -> SetMinimum(0.0);
  hWncdis  -> SetMinimum(0.0);

  hQ2      -> SetMinimum(0.0);
  hQ2cc    -> SetMinimum(0.0);
  hQ2ccqel -> SetMinimum(0.0);
  hQ2ccres -> SetMinimum(0.0);
  hQ2ccdis -> SetMinimum(0.0);
  hQ2nc    -> SetMinimum(0.0);
  hQ2ncqel -> SetMinimum(0.0);
  hQ2ncres -> SetMinimum(0.0);
  hQ2ncdis -> SetMinimum(0.0);

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

  ostringstream psfilename;
  psfilename << "mc-sample-test-" << gOptTgtPdgC << ".ps";
  gPS = new TPostScript(psfilename.str().c_str(), 112);

  // --- front page
  //
  AddFrontPage();

  // --- event fraction plots
  //
  //PlotEvtFrac("EvtFracNumuDir", "Event fractions");

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

  // --- hadronic 4-momentum plots
  //
  PlotHP4("HadP4Dir",             "All events"       );
  PlotHP4("HadP4CcNumuDir",       "#nu_{#mu} CC"     );
  PlotHP4("HadP4CcNumuQelDir",    "#nu_{#mu} CC QEL" );
  PlotHP4("HadP4CcNumuResDir",    "#nu_{#mu} CC RES" );
  PlotHP4("HadP4CcNumuDisDir",    "#nu_{#mu} CC DIS" );
  PlotHP4("HadP4NcNumuDir",       "#nu_{#mu} NC"     );
  PlotHP4("HadP4NcNumuQelDir",    "#nu_{#mu} NC QEL" );
  PlotHP4("HadP4NcNumuResDir",    "#nu_{#mu} NC RES" );
  PlotHP4("HadP4NcNumuDisDir",    "#nu_{#mu} NC DIS" );

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

  PlotH1F ( "hx",        (title + ", x_{comp}").c_str());
  PlotH1F ( "hy",        (title + ", y_{comp}").c_str());
  PlotH1F ( "ht",        (title + ", t_{comp}").c_str());
  PlotH1F ( "hW",        (title + ", W_{comp}").c_str());
  PlotH1F ( "hlog10Q2",  (title + ", log_{10}Q^{2}_{comp} (GeV^{2})").c_str());
  PlotH1F ( "hv",        (title + ", v_{comp}").c_str());
  PlotH1F ( "hxs",       (title + ", x_{sel}").c_str());
  PlotH1F ( "hys",       (title + ", y_{sel}").c_str());
  PlotH1F ( "hts",       (title + ", t_{sel}").c_str());
  PlotH1F ( "hWs",       (title + ", W_{sel}").c_str());
  PlotH1F ( "hlog10Q2s", (title + ", log_{10}Q^{2}_{sel} (GeV^{2})").c_str());
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
void PlotHP4(string dir, string title)
{
  gTestedSampleDir = (TDirectory*) gTestedSamplePlotFile->Get(dir.c_str());
  gTempltSampleDir = (gSampleComp) ? 
                     (TDirectory*) gTempltSamplePlotFile->Get(dir.c_str()) : 0;

  PlotH1F_2 ( "hpx_pip", "hpx_pip_hprim", (title + ", f/s #pi^{+} Px (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hpy_pip", "hpy_pip_hprim", (title + ", f/s #pi^{+} Py (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hpz_pip", "hpz_pip_hprim", (title + ", f/s #pi^{+} Pz (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hp_pip",  "hp_pip_hprim",  (title + ", f/s #pi^{+} P  (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hE_pip",  "hE_pip_hprim",  (title + ", f/s #pi^{+} E  (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hKE_pip", "hKE_pip_hprim", (title + ", f/s #pi^{+} KE (dashed/open: before hadron transport)").c_str() );

  PlotH1F_2 ( "hpx_pim", "hpx_pim_hprim", (title + ", f/s #pi^{-} Px (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hpy_pim", "hpy_pim_hprim", (title + ", f/s #pi^{-} Py (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hpz_pim", "hpz_pim_hprim", (title + ", f/s #pi^{-} Pz (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hp_pim",  "hp_pim_hprim",  (title + ", f/s #pi^{-} P  (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hE_pim",  "hE_pim_hprim",  (title + ", f/s #pi^{-} E  (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hKE_pim", "hKE_pim_hprim", (title + ", f/s #pi^{-} KE (dashed/open: before hadron transport)").c_str() );

  PlotH1F_2 ( "hpx_p",   "hpx_p_hprim",   (title + ", f/s proton Px (dashed/open: before hadron transport)").c_str()  );
  PlotH1F_2 ( "hpy_p",   "hpy_p_hprim",   (title + ", f/s proton Py (dashed/open: before hadron transport)").c_str()  );
  PlotH1F_2 ( "hpz_p",   "hpz_p_hprim",   (title + ", f/s proton Pz (dashed/open: before hadron transport)").c_str()  );
  PlotH1F_2 ( "hp_p",    "hp_p_hprim",    (title + ", f/s proton P  (dashed/open: before hadron transport)").c_str()  );
  PlotH1F_2 ( "hE_p",    "hE_p_hprim",    (title + ", f/s proton E  (dashed/open: before hadron transport)").c_str()  );
  PlotH1F_2 ( "hKE_p",   "hKE_p_hprim",   (title + ", f/s proton KE (dashed/open: before hadron transport)").c_str()  );

  PlotH1F_2 ( "hpx_n",   "hpx_n_hprim",   (title + ", f/s neutron Px (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hpy_n",   "hpy_n_hprim",   (title + ", f/s neutron Py (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hpz_n",   "hpz_n_hprim",   (title + ", f/s neutron Pz (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hp_n",    "hp_n_hprim",    (title + ", f/s neutron P  (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hE_n",    "hE_n_hprim",    (title + ", f/s neutron E  (dashed/open: before hadron transport)").c_str() );
  PlotH1F_2 ( "hKE_n",   "hKE_n_hprim",   (title + ", f/s neutron KE (dashed/open: before hadron transport)").c_str() );
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

  if(!tested_hst) return;
/*
  if(tested_hst->GetEntries() <= 0) {
    TPavesText warn(0.1, 0.6, 0.9, 0.9, 0, "tr");
    warn.SetTextSize(0.04);
    warn.AddText("EMPTY HISTOGRAM");
    warn.AddText(title.c_str());
    warn.SetFillColor(46);  
    warn.SetTextColor(10);  
    warn.Draw();
    return;
  }
*/
  // plot histogram from test sample
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
void PlotH1F_2(string name1, string name2, string title, bool keep_page)
{
  if(!keep_page) gPS->NewPage();
  gC->cd();

  TH1F * tested_hst1 = dynamic_cast<TH1F *> (gTestedSampleDir->Get(name1.c_str()));
  TH1F * tested_hst2 = dynamic_cast<TH1F *> (gTestedSampleDir->Get(name2.c_str()));
  TH1F * templt_hst1 = 0;
  TH1F * templt_hst2 = 0;
  if(gSampleComp) {
      templt_hst1 = dynamic_cast<TH1F *> (gTempltSampleDir->Get(name1.c_str()));
      templt_hst2 = dynamic_cast<TH1F *> (gTempltSampleDir->Get(name2.c_str()));
  }
  if(!tested_hst1 || !tested_hst2) return;

  // plot histogram from test sample
  tested_hst1->SetLineColor(2);
  tested_hst1->SetLineWidth(2);
  tested_hst1->SetLineStyle(1);
  tested_hst1->Draw();

  tested_hst2->SetLineColor(2);
  tested_hst2->SetLineWidth(2);
  tested_hst2->SetLineStyle(2);
  tested_hst2->Draw("SAME");

  // plot same hisogram from template sample (if any) 
  if(templt_hst1 && templt_hst2) {
    templt_hst1->SetLineWidth(2);
    templt_hst1->SetLineStyle(1);
    templt_hst1->SetMarkerSize(1.3);
    templt_hst1->SetMarkerStyle(8);
    templt_hst1->Draw("PERRSAME");

    templt_hst2->SetLineWidth(2);
    templt_hst2->SetLineStyle(2);
    templt_hst2->SetMarkerSize(1.3);
    templt_hst2->SetMarkerStyle(4);
    templt_hst2->Draw("PERRSAME");
  }

  tested_hst1->GetXaxis()->SetTitle(title.c_str());
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
  tested_hst->SetFillColor(3);
  tested_hst->Draw("BOX");

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
    hpxt->Draw("SAMEBOX1");
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
  if(filename.size() == 0) return false;

  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gmctest", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//_________________________________________________________________________________

