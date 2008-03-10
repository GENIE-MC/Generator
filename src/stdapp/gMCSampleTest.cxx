//____________________________________________________________________________
/*!

\program gmctest

\brief   A GENIE utility that reads a generated GENIE event tree and generates 
         a ROOT file with loads of histograms with basic MC truth quantities.
         It also generates the 'definite' summary ROOT ntuple.
         Typically used in validating generator changes.

         Syntax :
           gmctest -m mode -f sample [-n events] [-r reference_sample] 

         Options:
           [] Denotes an optional argument
           -f Specifies the GENIE/ROOT file with the generated event sample
           -n Specifies how many events to analyze [default: all]
	   -r Specifies another GENIE/ROOT event sample file for comparison 
           -m mode
	      0: create summary tree
	      1: create plots (and compare with reference sample, if any)
	      
\example gmctest -m 0 -f /path/gntp.1.ghep.root -n 10000
         gmctest -m 1 -f /path/gntp.1.gst.root -r /path/gntp.2.gst.root
		      
\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created September 02, 2005

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

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
#include <TLegend.h>

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
using std::vector;

using namespace genie;
using namespace genie::constants;

// function prototypes
void   GetCommandLineArgs   (int argc, char ** argv);
void   PrintSyntax          (void);
bool   CheckRootFilename    (string filename);
string OutputFileName       (string input_file_name, int mod);
void   CreateSummaryTree    (string filename);
void   CreatePlots          (string filename, string filename_ref);

// command-line arguments
Long64_t gOptN           = -1; // (-n) process so many events, all if -1
string   gOptInpFile     = ""; // (-f) input GENIE event sample file
string   gOptInpFileRef  = ""; // (-r) input GENIE event sample file (reference)
Int_t    gOptMode        = -1; // (-m) mode

// various control constants
const int kNPmax = 100;   

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- scan the command line arguments 
  GetCommandLineArgs(argc,argv);

  //-- analyze tested event generation sample
  if(gOptMode==0) {
     CreateSummaryTree(gOptInpFile);
  }
  
  if(gOptMode==1) {
     CreatePlots(gOptInpFile,gOptInpFileRef);
  }
  
  LOG("gmctest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void CreateSummaryTree(string inp_filename)
{
  if(!CheckRootFilename(inp_filename)) {
    LOG("gmctest", pERROR) << "Input file: " << inp_filename << " doesn't exist";
    return;
  }

  //-- define branch variables
  //
  int    brIev         = 0;      // Event number 
  int    brNeutrino    = 0;      // Neutrino pdg code
  int    brTarget      = 0;      // Nuclear target pdg code (10LZZZAAAI)
  int    brHitNuc      = 0;      // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
  int    brHitQrk      = 0;      // Hit quark pdg code        (set for DIS events only)
  bool   brFromSea     = false;  // Hit quark is from sea     (set for DIS events only)
  bool   brResId       = 0;      // Produced baryon resonance (set for resonance events only)
  bool   brIsQel       = false;  // Is QEL?
  bool   brIsRes       = false;  // Is RES?
  bool   brIsDis       = false;  // Is DIS?
  bool   brIsCohPi     = false;  // Is COHPi?
  bool   brIsCohEl     = false;  // Is COHEl?
  bool   brIsImd       = false;  // Is IMD?
  bool   brIsNuEL      = false;  // Is ve elastic?
  bool   brIsCC        = false;  // Is CC?
  bool   brIsNC        = false;  // Is NC?
  bool   brIsCharmPro  = false;  // Produces charm?
  double brWeight      = 0;      // Event weight
  double brKineXs      = 0;      // Bjorken x (selected)
  double brKineYs      = 0;      // Inelasticity y (selected)
  double brKineTs      = 0;      // Energy transfer to nucleus at COHPi events (selected)
  double brKineQ2s     = 0;      // Momentum transfer Q^2 (selected)
  double brKineWs      = 0;      // Hadronic invariant mass W (selected)
  double brKineX       = 0;      // Bjorken x  (computed from the event record)
  double brKineY       = 0;      // Inelasticity y (computed from the event record)
  double brKineT       = 0;      // Energy transfer to nucleus at COHPi events (computed from the event record)
  double brKineQ2      = 0;      // Momentum transfer Q^2 (computed from the event record)
  double brKineW       = 0;      // Hadronic invariant mass W (computed from the event record)
  double brEv          = 0;      // Neutrino energy (neutrino assumed in +z direction)
  double brEn          = 0;      // Initial state hit nucleon energy
  double brPxn         = 0;      // Initial state hit nucleon px
  double brPyn         = 0;      // Initial state hit nucleon py
  double brPzn         = 0;      // Initial state hit nucleon pz
  double brEl          = 0;      // Final state primary lepton energy
  double brPxl         = 0;      // Final state primary lepton px
  double brPyl         = 0;      // Final state primary lepton py
  double brPzl         = 0;      // Final state primary lepton pz 
  int    brNfP         = 0;      // Number of final state p's + \bar{p}'s (after intranuclear rescattering)
  int    brNfN         = 0;      // Number of final state n's + \bar{n}'s
  int    brNfPip       = 0;      // Number of final state pi+'s
  int    brNfPim       = 0;      // Number of final state pi-'s
  int    brNfPi0       = 0;      // Number of 'final state' pi0's (
  int    brNfKp        = 0;      // Number of final state K+'s
  int    brNfKm        = 0;      // Number of final state K-'s
  int    brNfK0        = 0;      // Number of final state K0's + \bar{K0}'s
  int    brNfEM        = 0;      // Number of final state gammas and e-/e+ (excluding pi0 decay products)
  int    brNfOther     = 0;      // Number of heavier final state hadrons (D+,D-,D0,Ds+,Ds-,Lamda,Sigma,Lamda_c,Sigma_c,...)
  int    brNiP         = 0;      // Number of 'primary' (: before intranuclear rescattering) p's + \bar{p}'s  
  int    brNiN         = 0;      // Number of 'primary' n's + \bar{n}'s  
  int    brNiPip       = 0;      // Number of 'primary' pi+'s 
  int    brNiPim       = 0;      // Number of 'primary' pi-'s 
  int    brNiPi0       = 0;      // Number of 'primary' pi0's 
  int    brNiKp        = 0;      // Number of 'primary' K+'s  
  int    brNiKm        = 0;      // Number of 'primary' K-'s  
  int    brNiK0        = 0;      // Number of 'primary' K0's + \bar{K0}'s 
  int    brNiEM        = 0;      // Number of 'primary' gammas and e-/e+ (eg from resonance decays)
  int    brNiOther     = 0;      // Number of 'primary' hadron shower particles
  int    brNf          = 0;      // Number of final state particles in hadronic system
  int    brPdgf[kNPmax];        // Pdg code of i^th final state particle in hadronic system
  double brEf  [kNPmax];        // Energy   of i^th final state particle in hadronic system
  double brPxf [kNPmax];        // Px       of i^th final state particle in hadronic system
  double brPyf [kNPmax];        // Py       of i^th final state particle in hadronic system
  double brPzf [kNPmax];        // Pz       of i^th final state particle in hadronic system
  int    brNi          = 0;      // Number of particles in 'primary' hadronic system (before intranuclear rescattering)
  int    brPdgi[kNPmax];        // Pdg code of i^th particle in 'primary' hadronic system 
  double brEi  [kNPmax];        // Energy   of i^th particle in 'primary' hadronic system 
  double brPxi [kNPmax];        // Px       of i^th particle in 'primary' hadronic system 
  double brPyi [kNPmax];        // Py       of i^th particle in 'primary' hadronic system 
  double brPzi [kNPmax];        // Pz       of i^th particle in 'primary' hadronic system 

  //-- build output filename based on input filename
  string outp_filename = OutputFileName(inp_filename,0);
  LOG("gmctest", pNOTICE) 
       << "*** Saving summary tree to: " << outp_filename;

  TFile fout(outp_filename.c_str(),"recreate");

  //-- create output summary tree & create the tree branches
  //
  TTree * s_tree = new TTree("gst","GENIE Summary Event Tree");

  //-- create tree branches
  //
  s_tree->Branch("iev",       &brIev,       "iev/I"       );
  s_tree->Branch("neu",	      &brNeutrino,  "neu/I"	  );
  s_tree->Branch("tgt",       &brTarget,    "tgt/I"	  );
  s_tree->Branch("hitnuc",    &brHitNuc,    "hitnuc/I"    );
  s_tree->Branch("hitqrk",    &brHitQrk,    "hitqrk/I"    );
  s_tree->Branch("resid",     &brResId,	    "resid/I"	  );
  s_tree->Branch("sea",	      &brFromSea,   "sea/O"	  );
  s_tree->Branch("qel",	      &brIsQel,	    "qel/O"	  );
  s_tree->Branch("res",	      &brIsRes,	    "res/O"	  );
  s_tree->Branch("dis",	      &brIsDis,	    "dis/O"	  );
  s_tree->Branch("cohpi",     &brIsCohPi,   "cohpi/O"	  );
  s_tree->Branch("cohel",     &brIsCohEl,   "cohel/O"	  );
  s_tree->Branch("imd",	      &brIsImd,	    "imd/O"	  );
  s_tree->Branch("nuel",      &brIsNuEL,    "nuel/O"	  );
  s_tree->Branch("cc",	      &brIsCC,	    "cc/O"	  );
  s_tree->Branch("nc",	      &brIsNC,	    "nc/O"	  );
  s_tree->Branch("charm",     &brIsCharmPro,"charm/O"	  );
  s_tree->Branch("wght",      &brWeight,    "wght/D"	  );
  s_tree->Branch("xs",	      &brKineXs,    "xs/D"	  );
  s_tree->Branch("ys",	      &brKineYs,    "ys/D"	  );
  s_tree->Branch("ts",	      &brKineTs,    "ts/D"	  );
  s_tree->Branch("Q2s",	      &brKineQ2s,   "Q2s/D"	  );
  s_tree->Branch("Ws",	      &brKineWs,    "Ws/D"	  );
  s_tree->Branch("x",	      &brKineX,	    "x/D"	  );
  s_tree->Branch("y",	      &brKineY,	    "y/D"	  );
  s_tree->Branch("t",	      &brKineT,	    "t/D"	  );
  s_tree->Branch("Q2",	      &brKineQ2,    "Q2/D"	  );
  s_tree->Branch("W",	      &brKineW,	    "W/D"	  );
  s_tree->Branch("Ev",	      &brEv,	    "Ev/D"	  );
  s_tree->Branch("En",	      &brEn,	    "En/D"	  );
  s_tree->Branch("pxn",	      &brPxn,	    "pxn/D"	  );
  s_tree->Branch("pyn",	      &brPyn,	    "pyn/D"	  );
  s_tree->Branch("pzn",	      &brPzn,	    "pzn/D"	  );
  s_tree->Branch("El",	      &brEl,	    "El/D"	  );
  s_tree->Branch("pxl",	      &brPxl,	    "pxl/D"	  );
  s_tree->Branch("pyl",	      &brPyl,	    "pyl/D"	  );
  s_tree->Branch("pzl",	      &brPzl,	    "pzl/D"	  );
  s_tree->Branch("nfp",	      &brNfP,	    "nfp/I"	  );
  s_tree->Branch("nfn",	      &brNfN,	    "nfn/I"	  );
  s_tree->Branch("nfpip",     &brNfPip,	    "nfpip/I"	  );
  s_tree->Branch("nfpim",     &brNfPim,	    "nfpim/I"	  );
  s_tree->Branch("nfpi0",     &brNfPi0,	    "nfpi0/I"	  );
  s_tree->Branch("nfkp",      &brNfKp,	    "nfkp/I"	  );
  s_tree->Branch("nfkm",      &brNfKm,	    "nfkm/I"	  );
  s_tree->Branch("nfk0",      &brNfK0,	    "nfk0/I"	  );
  s_tree->Branch("nfem",      &brNfEM,	    "nfem/I"	  );
  s_tree->Branch("nfother",   &brNfOther,   "nfother/I"   );
  s_tree->Branch("nip",	      &brNiP,	    "np/I"	  );
  s_tree->Branch("nin",	      &brNiN,	    "nn/I"	  );
  s_tree->Branch("nipip",     &brNiPip,	    "npip/I"	  );
  s_tree->Branch("nipim",     &brNiPim,	    "npim/I"	  );
  s_tree->Branch("nipi0",     &brNiPi0,	    "npi0/I"	  );
  s_tree->Branch("nikp",      &brNiKp,	    "nkp/I"	  );
  s_tree->Branch("nikm",      &brNiKm,	    "nkm/I"	  );
  s_tree->Branch("nik0",      &brNiK0,	    "nk0/I"	  );
  s_tree->Branch("niem",      &brNiEM,	    "niem/I"	  );
  s_tree->Branch("niother",   &brNiOther,   "niother/I"   );
  s_tree->Branch("ni",	      &brNi,	    "ni/I"	  );
  s_tree->Branch("pdgi",       brPdgi,	    "pdgi[ni]/I " );
  s_tree->Branch("Ei",	       brEi,	    "Ei[ni]/D"    );
  s_tree->Branch("pxi",	       brPxi,	    "pxi[ni]/D"   );
  s_tree->Branch("pyi",	       brPyi,	    "pyi[ni]/D"   );
  s_tree->Branch("pzi",	       brPzi,	    "pzi[ni]/D"   );
  s_tree->Branch("nf",	      &brNf,	    "nf/I"	  );
  s_tree->Branch("pdgf",       brPdgf,	    "pdgf[nf]/I " );
  s_tree->Branch("Ef",	       brEf,	    "Ef[nf]/D"    );
  s_tree->Branch("pxf",	       brPxf,	    "pxf[nf]/D"   );
  s_tree->Branch("pyf",	       brPyf,	    "pyf[nf]/D"   );
  s_tree->Branch("pzf",	       brPzf,	    "pzf[nf]/D"   );

  //-- open the ROOT file and get the TTree & its header
  TFile fin(inp_filename.c_str(),"READ");
  TTree *           er_tree = 0;
  NtpMCTreeHeader * thdr    = 0;
  er_tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr    = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  if (!er_tree) {
    LOG("gmctest", pERROR) << "Null input ER tree";
    return;
  }

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //NtpMCFormat_t format = thdr->format;
  //assert(format == kNFGHEP);

  //-- get the mc record
  NtpMCEventRecord * mcrec = 0;
  er_tree->SetBranchAddress("gmcrec", &mcrec);

  if (!mcrec) {
    LOG("gmctest", pERROR) << "Null MC record";
    return;
  }
  
  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       er_tree->GetEntries() : TMath::Min( er_tree->GetEntries(), gOptN );

  if (nmax<0) {
    LOG("gmctest", pERROR) << "Number of events = 0";
    return;
  }
  
  LOG("gmctest", pNOTICE) << "*** Analyzing: " << nmax << " events";

  TLorentzVector pdummy(0,0,0,0);

  for(Long64_t iev = 0; iev < nmax; iev++) {

    er_tree->GetEntry(iev);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gmctest", pINFO) << rec_header;
    LOG("gmctest", pINFO) << event;

    // go further only if the event is physical
    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) {
      LOG("gmctest", pINFO) << "Skipping unphysical event";
      mcrec->Clear();
      continue;
    }

    // clean-up arrays
    //
    for(int j=0; j<kNPmax; j++) {
       brPdgi[j] = 0;     
       brEi  [j] = 0;     
       brPxi [j] = 0;     
       brPyi [j] = 0;     
       brPzi [j] = 0;     
       brPdgf[j] = 0;     
       brEf  [j] = 0;     
       brPxf [j] = 0;     
       brPyf [j] = 0;     
       brPzf [j] = 0;     
    }

    // Computing event characteristics
    //

    //input particles
    GHepParticle * neutrino = event.Probe();
    assert(neutrino);
    GHepParticle * target = event.Particle(1);
    assert(target);
    GHepParticle * fsl = event.FinalStatePrimaryLepton();
    assert(fsl);
    GHepParticle * hitnucl = event.HitNucleon();
  
    //summary info
    const Interaction * interaction = event.Summary();
    const InitialState & init_state = interaction->InitState();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    const Kinematics &   kine       = interaction->Kine();
    const XclsTag &      xcls       = interaction->ExclTag();
    const Target &       tgt        = init_state.Tgt();

    //process id
    bool is_qel    = proc_info.IsQuasiElastic();
    bool is_res    = proc_info.IsResonant();
    bool is_dis    = proc_info.IsDeepInelastic();
    bool is_cohpi  = proc_info.IsCoherentPiProd();
    bool is_cohel  = proc_info.IsCoherentElas();
    bool is_imd    = proc_info.IsInverseMuDecay();
    bool is_nuel   = proc_info.IsNuElectronElastic();
    bool is_weakcc = proc_info.IsWeakCC();
    bool is_weaknc = proc_info.IsWeakNC();
    bool is_coh    = is_cohpi | is_cohel;

    if(!hitnucl) { assert(is_coh || is_imd || is_nuel); }
  
    // hit quark 
    // set only for DIS events
    int  qrk  = (is_dis) ? tgt.HitQrkPdg() : 0;     
    bool seaq = (is_dis) ? tgt.HitSeaQrk() : false; 

    // resonance id ($GENIE/src/BaryonResonance/BaryonResonance.h)
    // set only for resonance neutrinoproduction
    int resid = (is_res) ? xcls.Resonance() : 0;

    // (qel or dis) charm production?
    bool charm = xcls.IsCharmEvent();
    
    //weight
    double weight = event.Weight();

    //input 4-momenta    
    const TLorentzVector & k1 = *(neutrino->P4());                     // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());                          // l 4-p (k2)
    const TLorentzVector & p1 = (hitnucl) ? *(hitnucl->P4()) : pdummy; // N 4-p (p1)      
     
    // compute kinematical params
    //
    double M  = kNucleonMass;
    TLorentzVector q  = k1-k2;                     // q=k1-k2, 4-p transfer
    double Q2 = -1 * q.M2();                       // momemtum transfer
    double v  = (hitnucl) ? q.Energy()       : -1; // v (E transfer in hit nucleon rest frame)
    double x  = (hitnucl) ? 0.5*Q2/(M*v)     : -1; // Bjorken x
    double y  = (hitnucl) ? v/k1.Energy()    : -1; // Inelasticity, y = q*P1/k1*P1
    double W2 = (hitnucl) ? M*M + 2*M*v - Q2 : -1; // Hadronic Invariant mass ^ 2
    double W  = (hitnucl) ? TMath::Sqrt(W2)  : -1; 
    double t  = 0;

    LOG("gmctest", pDEBUG) 
       << "[Calc] Q2 = " << Q2 << ", W = " << W 
       << ", x = " << x << ", y = " << y << ", t = " << t;

    // also, access kinematical params _exactly_ as they were selected internally
    // (possibly using off-shell kinematics)
    //
    bool get_selected = true;
    double xs  = kine.x (get_selected);
    double ys  = kine.y (get_selected);
    double ts  = (is_cohpi) ? kine.t (get_selected) : -1;
    double Q2s = kine.Q2(get_selected);
    double Ws  = kine.W (get_selected);

    LOG("gmctest", pDEBUG) 
       << "[Select] Q2 = " << Q2s << ", W = " << Ws 
       << ", x = " << xs << ", y = " << ys << ", t = " << ts;

    // Extract more info on the hadronic system
    // Only for QEL/RES/DIS/COH events
    //
    bool study_hadsyst = (is_qel || is_res || is_dis || is_coh);
    
    //
    TObjArrayIter piter(&event);
    GHepParticle * p = 0;
    int ip=-1;

    //
    // Extract info on the final state system originating from the
    // hadronic vertex (includes intranuclear rescattering mc)
    //
    // Notes:
    //  ** include f/s  p,n,\bar{p},\bar{n}
    //  ** include f/s pi+, pi-
    //  ** include **decayed** pi0 & ommit their decay products
    //  ** include f/s K+, K-, K0, \bar{K0}
    //  ** include gammas/e+/e- but not the ones coming from decaying pi0's (pi0's are counted)
    //  ** include f/s D+, D-, D0, \bar{D0}, Ds+, Ds-, Sigma's, Omega's, Lambda's, Sigma_{c}'s,...
    //  ** baryon resonances should have been decayed early on: include decay products
    //  ** eta,eta',rho0,rho+,rho-,omega,phi should have been decayed early on: include decay products
    //

    LOG("gmctest", pDEBUG) << "Extracting final state hadronic system";

    vector<int> final_had_syst;
    while( (p = (GHepParticle *) piter.Next()) && study_hadsyst)
    {
      ip++;
//      if(!is_coh && ip < TMath::Max(hitnucl->FirstDaughter(), event.FinalStatePrimaryLeptonPosition()+1)) continue;
      if(!is_coh && event.Particle(ip)->FirstMother()==0) continue;
      if(p->IsFake()) continue;
      int pdgc = p->Pdg();
      int ist  = p->Status();
      if(ist==kIStStableFinalState) {
         if (pdgc == kPdgGamma || pdgc == kPdgElectron || pdgc == kPdgPositron)  {
            int igmom = p->FirstMother();
            if(igmom!=-1) {
               if(event.Particle(igmom)->Pdg() != kPdgPi0) { final_had_syst.push_back(ip); }
            }
         } else {
            final_had_syst.push_back(ip);
         }
      }
      if(ist==kIStDecayedState && pdgc==kPdgPi0) {
         final_had_syst.push_back(ip);
      }
    }//particle-loop

    if( count(final_had_syst.begin(), final_had_syst.end(), -1) > 0) {
        mcrec->Clear();
 	continue;
    }

    //
    // Extract info on the primary hadronic system (before any intranuclear rescattering)
    // * For DIS: 
    //   Low-W events hadronized by KNO:
    //       Find the HadronicSyst special particle & get its daughters.
    //   High-W events hadronized by JETSET: 
    //       Find the HadronicSyst special particle & get its daughters. Find the JETSET
    //       special particle ('cluster','string','indep') and take its own daughters.
    //       Neglect particles decayed internally by JETSET
    // * For RES:
    //   Find the hit nucleon and lookup its 1st daughter (intermediate resonance).
    //   Get the resonance decay products.
    // * For QEL:
    //   Get the 1st daughter of the hit nucleon
    // * For other processes:
    //   Skip...
    //
    // For free nucleon targets (no intranuclear rescattering) the primary hadronic system
    // is 'identical' with the final state hadronic system
    //

    LOG("gmctest", pDEBUG) << "Extracting primary hadronic system";

    vector<int> prim_had_syst;
    if(study_hadsyst) {
      if(!target->IsNucleus() || (is_cohel||is_cohpi)) {
         vector<int>::const_iterator hiter = final_had_syst.begin();
         for( ; hiter != final_had_syst.end(); ++hiter) {
           prim_had_syst.push_back(*hiter);
         }
      } 
      else {
         int ihadbase=0;
         if(is_dis) {
           ihadbase = event.FinalStateHadronicSystemPosition();
           int idx = event.Particle(ihadbase)->LastDaughter() + 1;
           p = event.Particle(idx);
           if(p->Pdg()==kPdgCluster || p->Pdg()==kPdgString || p->Pdg()==kPdgIndep) ihadbase=idx;
         }
         if(is_qel || is_res) {
           ihadbase = hitnucl->FirstDaughter();
         }
         assert(ihadbase>0);

         int idx1 = event.Particle(ihadbase)->FirstDaughter();
         int idx2 = event.Particle(ihadbase)->LastDaughter();
         for(int i=idx1; i<=idx2; i++) {
            p = event.Particle(i);
            if(p->IsFake()) continue;
            int ist = p->Status();
            // handle decayed dis states
            if(is_dis && ist==kIStDISPreFragmHadronicState) {
               for(int j=p->FirstDaughter(); j<=p->LastDaughter(); j++) prim_had_syst.push_back(j);
            } 
            // handle decayed resonances (whose decay products may be resonances that decay further)
            else if(is_res && ist==kIStDecayedState) {
                for(int j=p->FirstDaughter(); j<=p->LastDaughter(); j++) {
                   GHepParticle * pd = event.Particle(j);
                   if(pd->Status()==kIStDecayedState) {
                      for(int k=pd->FirstDaughter(); k<=pd->LastDaughter(); k++) prim_had_syst.push_back(k);
                   } else {
                      prim_had_syst.push_back(j); 
                   }       
                }
            } else {
                 prim_had_syst.push_back(i);
            }
         }//i
         // also include gammas from nuclear de-excitations (appearing in the daughter list of the 
         // hit nucleus, earlier than the primary hadronic system extracted above)
         for(int i = target->FirstDaughter(); i <= target->LastDaughter(); i++) {
           if(i<0) continue;
           if(event.Particle(i)->Status()==kIStStableFinalState) { prim_had_syst.push_back(i); }
         }
      }//freenuc?
    }//study_hadsystem?

    if( count(prim_had_syst.begin(), prim_had_syst.end(), -1) > 0) {
        mcrec->Clear();
 	continue;
    }

    //
    // Al information has been assembled -- Start filling up the tree branches
    //
    brIev        = iev;      
    brNeutrino   = neutrino->Pdg();      
    brTarget     = target->Pdg();      
    brHitNuc     = (hitnucl) ? hitnucl->Pdg() : 0;      
    brHitQrk     = qrk;     
    brFromSea    = seaq;  
    brResId      = resid;
    brIsQel      = is_qel;
    brIsRes      = is_res;
    brIsDis      = is_dis;  
    brIsCohPi    = is_cohpi;  
    brIsCohEl    = is_cohel;  
    brIsImd      = is_imd;  
    brIsNuEL     = is_nuel;  
    brIsCC       = is_weakcc;  
    brIsNC       = is_weaknc;  
    brIsCharmPro = charm;
    brWeight     = weight;      
    brKineXs     = xs;      
    brKineYs     = ys;      
    brKineTs     = ts;      
    brKineQ2s    = Q2s;            
    brKineWs     = Ws;      
    brKineX      = x;      
    brKineY      = y;      
    brKineT      = t;      
    brKineQ2     = Q2;      
    brKineW      = W;      
    brEv         = k1.Energy();      
    brEn         = (hitnucl) ? p1.Energy() : 0;      
    brPxn        = (hitnucl) ? p1.Px()     : 0;      
    brPyn        = (hitnucl) ? p1.Py()     : 0;      
    brPzn        = (hitnucl) ? p1.Pz()     : 0;            
    brEl         = k2.Energy();      
    brPxl        = k2.Px();      
    brPyl        = k2.Py();      
    brPzl        = k2.Pz();      

    // prim had syst
    brNiP        = 0;
    brNiN        = 0;    
    brNiPip      = 0;    
    brNiPim      = 0;    
    brNiPi0      = 0;    
    brNiKp       = 0;  
    brNiKm       = 0;  
    brNiK0       = 0;  
    brNiEM       = 0;  
    brNiOther    = 0;  
    brNi = prim_had_syst.size();
    for(int j=0; j<brNi; j++) {
      p = event.Particle(prim_had_syst[j]);
      assert(p);
      brPdgi[j] = p->Pdg();     
      brEi  [j] = p->Energy();     
      brPxi [j] = p->Px();     
      brPyi [j] = p->Py();     
      brPzi [j] = p->Pz();     

      if      (p->Pdg() == kPdgProton  || p->Pdg() == kPdgAntiProton)   brNiP++;
      else if (p->Pdg() == kPdgNeutron || p->Pdg() == kPdgAntiNeutron)  brNiN++;
      else if (p->Pdg() == kPdgPiP) brNiPip++; 
      else if (p->Pdg() == kPdgPiM) brNiPim++; 
      else if (p->Pdg() == kPdgPi0) brNiPi0++; 
      else if (p->Pdg() == kPdgKP)  brNiKp++;  
      else if (p->Pdg() == kPdgKM)  brNiKm++;  
      else if (p->Pdg() == kPdgK0    || p->Pdg() == kPdgAntiK0)  brNiK0++; 
      else if (p->Pdg() == kPdgGamma || p->Pdg() == kPdgElectron || p->Pdg() == kPdgPositron) brNiEM++;
      else brNiOther++;

      LOG("gmctest", pINFO) 
        << "Counting in primary hadronic system: idx = " << prim_had_syst[j]
        << " -> " << p->Name();
    }

    LOG("gmctest", pINFO) 
     << "N(p):"             << brNiP
     << ", N(n):"           << brNiN
     << ", N(pi+):"         << brNiPip
     << ", N(pi-):"         << brNiPim
     << ", N(pi0):"         << brNiPi0
     << ", N(K+,K-,K0):"    << brNiKp+brNiKm+brNiK0
     << ", N(gamma,e-,e+):" << brNiEM
     << ", N(etc):"         << brNiOther << "\n";

    // f/s had syst
    brNfP        = 0;
    brNfN        = 0;    
    brNfPip      = 0;    
    brNfPim      = 0;    
    brNfPi0      = 0;    
    brNfKp       = 0;  
    brNfKm       = 0;  
    brNfK0       = 0;  
    brNfEM       = 0;  
    brNfOther    = 0;  

    brNf = final_had_syst.size();
    for(int j=0; j<brNf; j++) {
      p = event.Particle(final_had_syst[j]);
      assert(p);
      brPdgf[j] = p->Pdg();     
      brEf  [j] = p->Energy();     
      brPxf [j] = p->Px();     
      brPyf [j] = p->Py();     
      brPzf [j] = p->Pz();     

      if      (p->Pdg() == kPdgProton  || p->Pdg() == kPdgAntiProton)   brNfP++;
      else if (p->Pdg() == kPdgNeutron || p->Pdg() == kPdgAntiNeutron)  brNfN++;
      else if (p->Pdg() == kPdgPiP) brNfPip++; 
      else if (p->Pdg() == kPdgPiM) brNfPim++; 
      else if (p->Pdg() == kPdgPi0) brNfPi0++; 
      else if (p->Pdg() == kPdgKP)  brNfKp++;  
      else if (p->Pdg() == kPdgKM)  brNfKm++;  
      else if (p->Pdg() == kPdgK0    || p->Pdg() == kPdgAntiK0)  brNfK0++; 
      else if (p->Pdg() == kPdgGamma || p->Pdg() == kPdgElectron || p->Pdg() == kPdgPositron) brNfEM++;
      else brNfOther++;

      LOG("gmctest", pINFO) 
        << "Counting in f/s system from hadronic vtx: idx = " << final_had_syst[j]
        << " -> " << p->Name();
    }

    LOG("gmctest", pINFO) 
     << "N(p):"             << brNfP
     << ", N(n):"           << brNfN
     << ", N(pi+):"         << brNfPip
     << ", N(pi-):"         << brNfPim
     << ", N(pi0):"         << brNfPi0
     << ", N(K+,K-,K0):"    << brNfKp+brNfKm+brNfK0
     << ", N(gamma,e-,e+):" << brNfEM
     << ", N(etc):"         << brNfOther << "\n";

    s_tree->Fill();

    mcrec->Clear();

  } // event loop

  fin.Close();

  fout.Write();
  fout.Close();
}
//_________________________________________________________________________________
void CreatePlots(string inp_filename, string inp_filename_ref)
{
  TFile * fin_0 = 0;
  TFile * fin_1 = 0;
  TTree * gst_0 = 0;
  TTree * gst_1 = 0;

  if(!CheckRootFilename(inp_filename)) {
    LOG("gmctest", pERROR) << "Input file: " << inp_filename << " doesn't exist";
    return;
  }
  fin_0 = new TFile(inp_filename.c_str(),"READ");
  gst_0 = (TTree *) fin_0->Get("gst");
  assert(gst_0);
  
  if(CheckRootFilename(inp_filename_ref)) {
     fin_1 = new TFile(inp_filename_ref.c_str(),"READ");
     gst_1 = (TTree *) fin_1->Get("gst");
     assert(gst_1);
  }
  
  gst_0->SetLineColor(kBlack);
  gst_0->SetLineWidth(3);
  if(gst_1) {
    gst_1->SetLineColor(kRed);
    gst_1->SetMarkerColor(kRed);
    gst_1->SetLineWidth(2);
    gst_1->SetMarkerStyle(20);
    gst_1->SetMarkerSize(1);
  }
  
  // Set global plot style
  //
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetHistTopMargin(0.33);
  gStyle->SetHistMinimumZero(true);
  
  // Plotting options
  //
  bool monoenergetic_sample = true;
  bool show_coh_plots       = true;
  bool show_calc_kinematics = true;
  bool show_mult_per_proc   = true;
  bool show_primary_hadsyst = true;
  
  gst_0->Draw("1","tgt>1000010010","GOFF");
  show_coh_plots = (gst_0->GetSelectedRows() > 0); 

  TCanvas * c = new TCanvas("c","",20,20,500,650);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  c->SetGridx();
  c->SetGridy();

  TLegend * ls = new TLegend(0.20,0.94,0.99,0.99);
  ls->SetFillColor(0);
  ls->SetBorderSize(0);

  string ps_filename = OutputFileName(inp_filename,1);
  TPostScript * ps = new TPostScript(ps_filename.c_str(), 111);
  
  //
  // SECTION: PS File Header
  //

  ps->NewPage();
  c->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText("GENIE Event Sample Comparisons");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText("Event Sample:");
  //hdr.AddText(sample);
  hdr.AddText(" ");
  hdr.AddText("Notes:");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.Draw();
  c->Update();

  //
  // SECTION: Event Numbers
  //

  TH1F * h0num = new TH1F("h0num","", 2, 0., 2.);
  TH1F * h1num = new TH1F("h1num","", 2, 0., 2.);
  float n0=0, n1=0;
  ps->NewPage();
  c->Range(0,0,100,100);
  TPavesText evn(10,10,90,90,3,"tr");
  evn.AddText("Event Numbers:");
  evn.AddText("  ");
 
            gst_0->Draw("1>>h0num","1","goff");
  if(gst_1) gst_1->Draw("1>>h1num","1","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("ALL    : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","qel","goff");
  if(gst_1) gst_1->Draw("1>>h1num","qel","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("QEL    : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","qel&&cc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","qel&&cc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("QEL-CC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","qel&&nc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","qel&&nc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("QEL-NC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","res","goff");
  if(gst_1) gst_1->Draw("1>>h1num","res","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("RES    : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","res&&cc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","res&&cc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("RES-CC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","res&&nc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","res&&nc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("RES-NC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","dis","goff");
  if(gst_1) gst_1->Draw("1>>h1num","dis","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("DIS    : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","dis&&cc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","dis&&cc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("DIS-CC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","dis&&cc&&charm","goff");
  if(gst_1) gst_1->Draw("1>>h1num","dis&&cc&&charm","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("(charm): %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","dis&&nc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","dis&&nc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("DIS-NC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","cohpi","goff");
  if(gst_1) gst_1->Draw("1>>h1num","cohpi","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("COHPi    : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","cohpi&&cc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","cohpi&&cc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("COHPi-CC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","cohpi&&nc","goff");
  if(gst_1) gst_1->Draw("1>>h1num","cohpi&&nc","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("COHPi-NC : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

            gst_0->Draw("1>>h0num","imd","goff");
  if(gst_1) gst_1->Draw("1>>h1num","imd","goff");
  n0 =           h0num->GetEntries();
  n1 = (gst_1) ? h1num->GetEntries() : 0;
  evn.AddText( Form("IMD    : %7.0f [test sample], %7.0f [ref sample]", n0, n1) );

  evn.Draw();
  c->Update();

  if(!monoenergetic_sample) {
              gst_0->Draw("Ev","","");
    if(gst_1) gst_1->Draw("Ev","","perrsame");
    ls->Clear();
    ls->SetHeader("Neutrino Energy Spectrum");
    ls->Draw();
    c->Update();  
  }

  //
  // SECTION: Kinematics
  //
  ps->NewPage();
  c->Clear();
  c->Range(0,0,100,100);
  TPavesText hdrk(10,40,90,70,3,"tr");
  hdrk.AddText("Selected Kinematical Quantities");
  hdrk.AddText(" ");
  hdrk.Draw();
  c->Update();

  //------ selected Q2 for all events
  ps->NewPage();
  gst_0->Draw("Q2s","","");
  if(gst_1) gst_1->Draw("Q2s","","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for all events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for QEL
  ps->NewPage();
  gst_0->Draw("Q2s","qel","");
  if(gst_1) gst_1->Draw("Q2s","qel","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for QEL events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for QEL CC
  ps->NewPage();
  gst_0->Draw("Q2s","qel&&cc","");
  if(gst_1) gst_1->Draw("Q2s","qel&&cc","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for QEL CC events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for QEL NC
  ps->NewPage();
  gst_0->Draw("Q2s","qel&&nc","");
  if(gst_1) gst_1->Draw("Q2s","qel&&nc","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for QEL NC events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for RES
  ps->NewPage();
  gst_0->Draw("Q2s","res","");
  if(gst_1) gst_1->Draw("Q2s","res","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for RES events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for RES CC
  ps->NewPage();
  gst_0->Draw("Q2s","res&&cc","");
  if(gst_1) gst_1->Draw("Q2s","res&&cc","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for RES CC events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for RES NC
  ps->NewPage();
  gst_0->Draw("Q2s","res&&nc","");
  if(gst_1) gst_1->Draw("Q2s","res&&nc","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for RES NC events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for DIS
  ps->NewPage();
  gst_0->Draw("Q2s","dis","");
  if(gst_1) gst_1->Draw("Q2s","dis","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for DIS events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for DIS CC
  ps->NewPage();
  gst_0->Draw("Q2s","dis&&cc","");
  if(gst_1) gst_1->Draw("Q2s","dis&&cc","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for DIS CC events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for DIS NC
  ps->NewPage();
  gst_0->Draw("Q2s","dis&&nc","");
  if(gst_1) gst_1->Draw("Q2s","dis&&nc","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for DIS NC events");
  ls->Draw();
  c->Update();

  //------ selected Q2 for Charm/DIS
  ps->NewPage();
  gst_0->Draw("Q2s","dis&&charm","");
  if(gst_1) gst_1->Draw("Q2s","dis&&charm","perrsame");
  ls->Clear();
  ls->SetHeader("selected Q2 for Charm/DIS events");
  ls->Draw();
  c->Update();

  if(show_coh_plots) {
     //------ selected Q2 for COHPi
     ps->NewPage();
     gst_0->Draw("Q2s","cohpi","");
     if(gst_1) gst_1->Draw("Q2s","cohpi","perrsame");
     ls->Clear();
     ls->SetHeader("selected Q2 for COHPi events");
     ls->Draw();
     c->Update();

     //------ selected Q2 for COHPi CC
     ps->NewPage();
     gst_0->Draw("Q2s","cohpi&&cc","");
     if(gst_1) gst_1->Draw("Q2s","cohpi&&cc","perrsame");
     ls->Clear();
     ls->SetHeader("selected Q2 for COHPi CC events");
     ls->Draw();
     c->Update();

     //------ selected Q2 for COHPi NC
     ps->NewPage();
     gst_0->Draw("Q2s","cohpi&&nc","");
     if(gst_1) gst_1->Draw("Q2s","cohpi&&nc","perrsame");
     ls->Clear();
     ls->SetHeader("selected Q2 for COHPi NC events");
     ls->Draw();
     c->Update();
  }
  
  //------ selected W for all events
  ps->NewPage();
  gst_0->Draw("Ws","","");
  if(gst_1) gst_1->Draw("Ws","","perrsame");
  ls->Clear();
  ls->SetHeader("selected W for all events");
  ls->Draw();
  c->Update();

  //------ selected W for QEL
  ps->NewPage();
  gst_0->Draw("Ws","qel","");
  if(gst_1) gst_1->Draw("Ws","qel","perrsame");
  ls->Clear();
  ls->SetHeader("selected W for QEL events");
  ls->Draw();
  c->Update();

  //------ selected W for RES
  ps->NewPage();
  gst_0->Draw("Ws","res","");
  if(gst_1) gst_1->Draw("Ws","res","perrsame");
  ls->Clear();
  ls->SetHeader("selected W for RES events");
  ls->Draw();
  c->Update();

  //------ selected W for DIS
  ps->NewPage();
  gst_0->Draw("Ws","dis","");
  if(gst_1) gst_1->Draw("Ws","dis","perrsame");
  ls->Clear();
  ls->SetHeader("selected W for DIS events");
  ls->Draw();
  c->Update();

  //------ selected W for DIS CC
  ps->NewPage();
  gst_0->Draw("Ws","dis&&cc","");
  if(gst_1) gst_1->Draw("Ws","dis&&cc","perrsame");
  ls->Clear();
  ls->SetHeader("selected W for DIS CC events");
  ls->Draw();
  c->Update();

  //------ selected W for DIS NC
  ps->NewPage();
  gst_0->Draw("Ws","dis&&nc","");
  if(gst_1) gst_1->Draw("Ws","dis&&nc","perrsame");
  ls->Clear();
  ls->SetHeader("selected W for DIS NC events");
  ls->Draw();
  c->Update();

  //------ selected x for all events
  ps->NewPage();
  gst_0->Draw("xs","","");
  if(gst_1) gst_1->Draw("xs","","perrsame");
  ls->Clear();
  ls->SetHeader("selected x for all events");
  ls->Draw();
  c->Update();

  //------ selected x for QEL
  ps->NewPage();
  gst_0->Draw("xs","qel","");
  if(gst_1) gst_1->Draw("xs","qel","perrsame");
  ls->Clear();
  ls->SetHeader("selected x for QEL events");
  ls->Draw();
  c->Update();

  //------ selected x for RES
  ps->NewPage();
  gst_0->Draw("xs","res","");
  if(gst_1) gst_1->Draw("xs","res","perrsame");
  ls->Clear();
  ls->SetHeader("selected x for RES events");
  ls->Draw();
  c->Update();

  //------ selected x for DIS
  ps->NewPage();
  gst_0->Draw("xs","dis","");
  if(gst_1) gst_1->Draw("xs","dis","perrsame");
  ls->Clear();
  ls->SetHeader("selected x for DIS events");
  ls->Draw();
  c->Update();

  //------ selected x for DIS CC
  ps->NewPage();
  gst_0->Draw("xs","dis&&cc","");
  if(gst_1) gst_1->Draw("xs","dis&&cc","perrsame");
  ls->Clear();
  ls->SetHeader("selected x for DIS CC events");
  ls->Draw();
  c->Update();

  //------ selected x for DIS NC
  ps->NewPage();
  gst_0->Draw("xs","dis&&nc","");
  if(gst_1) gst_1->Draw("xs","dis&&nc","perrsame");
  ls->Clear();
  ls->SetHeader("selected x for DIS NC events");
  ls->Draw();
  c->Update();

  //------ selected x for Charm/DIS 
  ps->NewPage();
  gst_0->Draw("xs","dis&&charm","");
  if(gst_1) gst_1->Draw("xs","dis&&charm","perrsame");
  ls->Clear();
  ls->SetHeader("selected x for Charm/DIS events");
  ls->Draw();
  c->Update();

  if(show_coh_plots) {
     //------ selected x for COHPi
     ps->NewPage();
     gst_0->Draw("xs","cohpi","");
     if(gst_1) gst_1->Draw("xs","cohpi","perrsame");
     ls->Clear();
     ls->SetHeader("selected x for COHPi events");
     ls->Draw();
     c->Update();

     //------ selected x for COHPi CC
     ps->NewPage();
     gst_0->Draw("xs","cohpi&&cc","");
     if(gst_1) gst_1->Draw("xs","cohpi&&cc","perrsame");
     ls->Clear();
     ls->SetHeader("selected x for COHPi CC events");
     ls->Draw();
     c->Update();

     //------ selected x for COHPi NC
     ps->NewPage();
     gst_0->Draw("xs","cohpi&&nc","");
     if(gst_1) gst_1->Draw("xs","cohpi&&nc","perrsame");
     ls->Clear();
     ls->SetHeader("selected x for COHPi NC events");
     ls->Draw();
     c->Update();
  }

  //------ selected y for all events
  ps->NewPage();
  gst_0->Draw("ys","","");
  if(gst_1) gst_1->Draw("ys","","perrsame");
  ls->Clear();
  ls->SetHeader("selected y for all events");
  ls->Draw();
  c->Update();

  //------ selected y for QEL
  ps->NewPage();
  gst_0->Draw("ys","qel","");
  if(gst_1) gst_1->Draw("ys","qel","perrsame");
  ls->Clear();
  ls->SetHeader("selected y for QEL events");
  ls->Draw();
  c->Update();

  //------ selected y for RES
  ps->NewPage();
  gst_0->Draw("ys","res","");
  if(gst_1) gst_1->Draw("ys","res","perrsame");
  ls->Clear();
  ls->SetHeader("selected y for RES events");
  ls->Draw();
  c->Update();

  //------ selected y for DIS
  ps->NewPage();
  gst_0->Draw("ys","dis","");
  if(gst_1) gst_1->Draw("ys","dis","perrsame");
  ls->Clear();
  ls->SetHeader("selected y for DIS events");
  ls->Draw();
  c->Update();

  //------ selected y for DIS CC
  ps->NewPage();
  gst_0->Draw("ys","dis&&cc","");
  if(gst_1) gst_1->Draw("ys","dis&&cc","perrsame");
  ls->Clear();
  ls->SetHeader("selected y for DIS CC events");
  ls->Draw();
  c->Update();

  //------ selected y for DIS NC
  ps->NewPage();
  gst_0->Draw("ys","dis&&nc","");
  if(gst_1) gst_1->Draw("ys","dis&&nc","perrsame");
  ls->Clear();
  ls->SetHeader("selected y for DIS NC events");
  ls->Draw();
  c->Update();

  //------ selected y for Charm/DIS 
  ps->NewPage();
  gst_0->Draw("ys","dis&&charm","");
  if(gst_1) gst_1->Draw("ys","dis&&charm","perrsame");
  ls->Clear();
  ls->SetHeader("selected y for Charm/DIS events");
  ls->Draw();
  c->Update();

  if(show_coh_plots) {
     //------ selected y for COHPi
     ps->NewPage();
     gst_0->Draw("ys","cohpi","");
     if(gst_1) gst_1->Draw("ys","cohpi","perrsame");
     ls->Clear();
     ls->SetHeader("selected y for COHPi events");
     ls->Draw();
     c->Update();

     //------ selected y for COHPi CC
     ps->NewPage();
     gst_0->Draw("ys","cohpi&&cc","");
     if(gst_1) gst_1->Draw("ys","cohpi&&cc","perrsame");
     ls->Clear();
     ls->SetHeader("selected y for COHPi CC events");
     ls->Draw();
     c->Update();

     //------ selected y for COHPi NC
     ps->NewPage();
     gst_0->Draw("ys","cohpi&&nc","");
     if(gst_1) gst_1->Draw("ys","cohpi&&nc","perrsame");
     ls->Clear();
     ls->SetHeader("selected y for COHPi NC events");
     ls->Draw();
     c->Update();

     //------ selected t for COHPi
     ps->NewPage();
     gst_0->Draw("ts","cohpi","");
     if(gst_1) gst_1->Draw("ts","cohpi","perrsame");
     ls->Clear();
     ls->SetHeader("selected t for COHPi events");
     ls->Draw();
     c->Update();
  }

  if(show_calc_kinematics) {

     //
     // SECTION: Computed Kinematics 
     //
     ps->NewPage();
     c->Clear();
     c->Range(0,0,100,100);
     TPavesText hdrck(10,40,90,70,3,"tr");
     hdrck.AddText("Kinematical Quantities");
     hdrck.AddText(" ");
     hdrck.AddText(" ");
     hdrck.AddText("Similar to the previous set of plots but");
     hdrck.AddText("showing 'computed' rather than 'selected' variables");
     hdrck.Draw();
     c->Update();

     //------ Q2 for all events
     ps->NewPage();
               gst_0->Draw("Q2","","");
     if(gst_1) gst_1->Draw("Q2","","perrsame");
     ls->Clear();
     ls->SetHeader("computed Q2 for all events");
     ls->Draw();
     c->Update();

     //------ Q2 for QEL
     ps->NewPage();
               gst_0->Draw("Q2","qel","");
     if(gst_1) gst_1->Draw("Q2","qel","perrsame");
     ls->Clear();
     ls->SetHeader("computed Q2 for QEL events");
     ls->Draw();
     c->Update();

     //------ Q2 for RES
     ps->NewPage();
               gst_0->Draw("Q2","res","");
     if(gst_1) gst_1->Draw("Q2","res","perrsame");
     ls->Clear();
     ls->SetHeader("computed Q2 for RES events");
     ls->Draw();
     c->Update();

     //------ Q2 for DIS
     ps->NewPage();
               gst_0->Draw("Q2","dis","");
     if(gst_1) gst_1->Draw("Q2","dis","perrsame");
     ls->Clear();
     ls->SetHeader("computed Q2 for DIS events");
     ls->Draw();
     c->Update();

     //------ x for all events
     ps->NewPage();
     gst_0->Draw("x","","");
     if(gst_1) gst_1->Draw("x","","perrsame");
     ls->Clear();
     ls->SetHeader("computed x for all events");
     ls->Draw();
     c->Update();

     //------ x for QEL
     ps->NewPage();
     gst_0->Draw("x","qel","");
     if(gst_1) gst_1->Draw("x","qel","perrsame");
     ls->Clear();
     ls->SetHeader("computed x for QEL events");
     ls->Draw();
     c->Update();

     //------ x for RES
     ps->NewPage();
     gst_0->Draw("x","res","");
     if(gst_1) gst_1->Draw("x","res","perrsame");
     ls->Clear();
     ls->SetHeader("computed x for RES events");
     ls->Draw();
     c->Update();

     //------ x for DIS
     ps->NewPage();
     gst_0->Draw("x","dis","");
     if(gst_1) gst_1->Draw("x","dis","perrsame");
     ls->Clear();
     ls->SetHeader("computed x for DIS events");
     ls->Draw();
     c->Update();

     //------ y for all events
     ps->NewPage();
     gst_0->Draw("y","","");
     if(gst_1) gst_1->Draw("y","","perrsame");
     ls->Clear();
     ls->SetHeader("computed y for all events");
     ls->Draw();
     c->Update();

     //------ y for QEL
     ps->NewPage();
     gst_0->Draw("y","qel","");
     if(gst_1) gst_1->Draw("y","qel","perrsame");
     ls->Clear();
     ls->SetHeader("computed y for QEL events");
     ls->Draw();
     c->Update();

     //------ y for RES
     ps->NewPage();
     gst_0->Draw("y","res","");
     if(gst_1) gst_1->Draw("y","res","perrsame");
     ls->Clear();
     ls->SetHeader("computed y for RES events");
     ls->Draw();
     c->Update();

     //------ y for DIS
     ps->NewPage();
     gst_0->Draw("y","dis","");
     if(gst_1) gst_1->Draw("y","dis","perrsame");
     ls->Clear();
     ls->SetHeader("computed y for DIS events");
     ls->Draw();
     c->Update();

  }//show?


  //
  // SECTION: Initial State nucleon
  //
  ps->NewPage();
  c->Clear();
  c->Range(0,0,100,100);
  TPavesText hdrinuc(10,40,90,70,3,"tr");
  hdrinuc.AddText("Initial state nucleon 4-Momentum");
  hdrinuc.Draw();
  c->Update();

  //------ selected hit nucleon px
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxn","","");
  if(gst_1) gst_1->Draw("pxn","","perrsame");
  c->cd(2);
  gst_0->Draw("pyn","","");
  if(gst_1) gst_1->Draw("pyn","","perrsame");
  c->cd(3);
  gst_0->Draw("pzn","","");
  if(gst_1) gst_1->Draw("pzn","","perrsame");
  c->cd(4);
  gst_0->Draw("En","En>.2","");
  if(gst_1) gst_1->Draw("En","En>.2","perrsame");
  c->Update();

  //
  // SECTION: Final State Primary Lepton
  //
  ps->NewPage();
  c->Clear();
  c->Range(0,0,100,100);
  TPavesText hdrfsl(10,40,90,70,3,"tr");
  hdrfsl.AddText("Final State Primary Lepton 4-Momentum");
  hdrfsl.Draw();
  c->Update();

  //------ f/s primary lepton : all events
  ps->NewPage();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxl","","");
  if(gst_1) gst_1->Draw("pxl","","perrsame");
  c->cd(2);
  gst_0->Draw("pyl","","");
  if(gst_1) gst_1->Draw("pyl","","perrsame");
  c->cd(3);
  gst_0->Draw("pzl","","");
  if(gst_1) gst_1->Draw("pzl","","perrsame");
  c->cd(4);
  gst_0->Draw("El","","");
  if(gst_1) gst_1->Draw("El","","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state primary lepton 4-p: All events");
  ls->Draw();
  c->Update();

  //------ f/s primary lepton : all CC events
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxl","cc","");
  if(gst_1) gst_1->Draw("pxl","cc","perrsame");
  c->cd(2);
  gst_0->Draw("pyl","cc","");
  if(gst_1) gst_1->Draw("pyl","cc","perrsame");
  c->cd(3);
  gst_0->Draw("pzl","cc","");
  if(gst_1) gst_1->Draw("pzl","cc","perrsame");
  c->cd(4);
  gst_0->Draw("El","cc","");
  if(gst_1) gst_1->Draw("El","cc","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state primary lepton 4-p: All CC events");
  ls->Draw();
  c->Update();

  //------ f/s primary lepton : all NC events
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxl","nc","");
  if(gst_1) gst_1->Draw("pxl","nc","perrsame");
  c->cd(2);
  gst_0->Draw("pyl","nc","");
  if(gst_1) gst_1->Draw("pyl","nc","perrsame");
  c->cd(3);
  gst_0->Draw("pzl","nc","");
  if(gst_1) gst_1->Draw("pzl","nc","perrsame");
  c->cd(4);
  gst_0->Draw("El","nc","");
  if(gst_1) gst_1->Draw("El","nc","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state primary lepton 4-p: All NC events");
  ls->Draw();
  c->Update();

  //------ f/s primary lepton : QEL events
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxl","qel","");
  if(gst_1) gst_1->Draw("pxl","qel","perrsame");
  c->cd(2);
  gst_0->Draw("pyl","qel","");
  if(gst_1) gst_1->Draw("pyl","qel","perrsame");
  c->cd(3);
  gst_0->Draw("pzl","qel","");
  if(gst_1) gst_1->Draw("pzl","qel","perrsame");
  c->cd(4);
  gst_0->Draw("El","qel","");
  if(gst_1) gst_1->Draw("El","qel","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state primary lepton 4-p: QEL events");
  ls->Draw();
  c->Update();

  //------ f/s primary lepton : RES events
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxl","res","");
  if(gst_1) gst_1->Draw("pxl","res","perrsame");
  c->cd(2);
  gst_0->Draw("pyl","res","");
  if(gst_1) gst_1->Draw("pyl","res","perrsame");
  c->cd(3);
  gst_0->Draw("pzl","res","");
  if(gst_1) gst_1->Draw("pzl","res","perrsame");
  c->cd(4);
  gst_0->Draw("El","res","");
  if(gst_1) gst_1->Draw("El","res","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state primary lepton 4-p: RES events");
  ls->Draw();
  c->Update();

  //------ f/s primary lepton : DIS events
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxl","dis","");
  if(gst_1) gst_1->Draw("pxl","dis","perrsame");
  c->cd(2);
  gst_0->Draw("pyl","dis","");
  if(gst_1) gst_1->Draw("pyl","dis","perrsame");
  c->cd(3);
  gst_0->Draw("pzl","dis","");
  if(gst_1) gst_1->Draw("pzl","dis","perrsame");
  c->cd(4);
  gst_0->Draw("El","dis","");
  if(gst_1) gst_1->Draw("El","dis","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state primary lepton 4-p: All DIS events");
  ls->Draw();
  c->Update();

  if(show_coh_plots) {
     //------ f/s primary lepton : COHPi events
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxl","cohpi","");
     if(gst_1) gst_1->Draw("pxl","cohpi","perrsame");
     c->cd(2);
     gst_0->Draw("pyl","cohpi","");
     if(gst_1) gst_1->Draw("pyl","cohpi","perrsame");
     c->cd(3);
     gst_0->Draw("pzl","cohpi","");
     if(gst_1) gst_1->Draw("pzl","cohpi","perrsame");
     c->cd(4);
     gst_0->Draw("El","cohpi","");
     if(gst_1) gst_1->Draw("El","cohpi","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state primary lepton 4-p: COHPi events");
     ls->Draw();
     c->Update();
  }

  //
  // SECTION: Final State Hadronic System Multiplicities & 4P
  //
  ps->NewPage();
  c->Clear();
  c->Range(0,0,100,100);
  TPavesText hdrfhad(10,40,90,70,3,"tr");
  hdrfhad.AddText("Final State Hadronic System");
  hdrfhad.AddText("Multiplicities and 4-Momenta");
  hdrfhad.AddText(" ");
  hdrfhad.AddText(" ");
  hdrfhad.AddText(" ");
  hdrfhad.AddText(" ");
  hdrfhad.AddText("Note:");
  hdrfhad.AddText("For nuclear targets these plots include the effect");
  hdrfhad.AddText("of intranuclear hadron transport / rescattering");
  hdrfhad.Draw();
  c->Update();

  //------ number of final state p
  ps->NewPage();
  gst_0->Draw("nfp","","");
  if(gst_1) gst_1->Draw("nfp","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state protons");
  ls->Draw();
  c->Update();

  //------ number of final state n
  ps->NewPage();
  gst_0->Draw("nfn","","");
  if(gst_1) gst_1->Draw("nfn","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state neutrons");
  ls->Draw();
  c->Update();

  //------ number of final state pi+
  ps->NewPage();
  gst_0->Draw("nfpip","","");
  if(gst_1) gst_1->Draw("nfpip","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state pi+");
  ls->Draw();
  c->Update();

  //------ number of final state pi-
  ps->NewPage();
  gst_0->Draw("nfpim","","");
  if(gst_1) gst_1->Draw("nfpim","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state pi-");
  ls->Draw();
  c->Update();

  //------ number of final state pi0
  ps->NewPage();
  gst_0->Draw("nfpi0","","");
  if(gst_1) gst_1->Draw("nfpi0","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state pi0");
  ls->Draw();
  c->Update();

  //------ number of final state K+
  ps->NewPage();
  gst_0->Draw("nfkp","","");
  if(gst_1) gst_1->Draw("nfkp","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state K+");
  ls->Draw();
  c->Update();

  //------ number of final state K-
  ps->NewPage();
  gst_0->Draw("nfkm","","");
  if(gst_1) gst_1->Draw("nfkm","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state K-");
  ls->Draw();
  c->Update();

  //------ number of final state K0
  ps->NewPage();
  gst_0->Draw("nfk0","","");
  if(gst_1) gst_1->Draw("nfk0","","perrsame");
  ls->Clear();
  ls->SetHeader("Number of final state K0");
  ls->Draw();
  c->Update();

  //------ momentum of final state p
  ps->NewPage();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxf","pdgf==2212","");
  if(gst_1) gst_1->Draw("pxf","pdgf==2212","perrsame");
  c->cd(2);
  gst_0->Draw("pyf","pdgf==2212","");
  if(gst_1) gst_1->Draw("pyf","pdgf==2212","perrsame");
  c->cd(3);
  gst_0->Draw("pzf","pdgf==2212","");
  if(gst_1) gst_1->Draw("pzf","pdgf==2212","perrsame"); 
  c->cd(4);
  gst_0->Draw("Ef","pdgf==2212","");
  if(gst_1) gst_1->Draw("Ef","pdgf==2212","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state protons 4-momentum");
  ls->Draw();
  c->Update();

  //------ momentum of final state n
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxf","pdgf==2112","");
  if(gst_1) gst_1->Draw("pxf","pdgf==2112","perrsame");
  c->cd(2);
  gst_0->Draw("pyf","pdgf==2112","");
  if(gst_1) gst_1->Draw("pyf","pdgf==2112","perrsame");
  c->cd(3);
  gst_0->Draw("pzf","pdgf==2112","");
  if(gst_1) gst_1->Draw("pzf","pdgf==2112","perrsame");
  c->cd(4);
  gst_0->Draw("Ef","pdgf==2112","");
  if(gst_1) gst_1->Draw("Ef","pdgf==2112","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state neutrons 4-momentum");
  ls->Draw();
  c->Update();

  //------ momentum of final state pi0
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxf","pdgf==111","");
  if(gst_1) gst_1->Draw("pxf","pdgf==111","perrsame");
  c->cd(2);
  gst_0->Draw("pyf","pdgf==111","");
  if(gst_1) gst_1->Draw("pyf","pdgf==111","perrsame");
  c->cd(3);
  gst_0->Draw("pzf","pdgf==111","");
  if(gst_1) gst_1->Draw("pzf","pdgf==111","perrsame");
  c->cd(4);
  gst_0->Draw("Ef","pdgf==111","");
  if(gst_1) gst_1->Draw("Ef","pdgf==111","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state pi0's 4-momentum");
  ls->Draw();
  c->Update();

  //------ momentum of final state pi+
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxf","pdgf==211","");
  if(gst_1) gst_1->Draw("pxf","pdgf==211","perrsame");
  c->cd(2);
  gst_0->Draw("pyf","pdgf==211","");
  if(gst_1) gst_1->Draw("pyf","pdgf==211","perrsame");
  c->cd(3);
  gst_0->Draw("pzf","pdgf==211","");
  if(gst_1) gst_1->Draw("pzf","pdgf==211","perrsame");
  c->cd(4);
  gst_0->Draw("Ef","pdgf==211","");
  if(gst_1) gst_1->Draw("Ef","pdgf==211","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state pi+'s 4-momentum");
  ls->Draw();
  c->Update();

  //------ momentum of final state pi+
  ps->NewPage();
  c->Clear();
  c->Divide(2,2);
  c->cd(1);
  gst_0->Draw("pxf","pdgf==-211","");
  if(gst_1) gst_1->Draw("pxf","pdgf==-211","perrsame");
  c->cd(2);
  gst_0->Draw("pyf","pdgf==-211","");
  if(gst_1) gst_1->Draw("pyf","pdgf==-211","perrsame");
  c->cd(3);
  gst_0->Draw("pzf","pdgf==-211","");
  if(gst_1) gst_1->Draw("pzf","pdgf==-211","perrsame");
  c->cd(4);
  gst_0->Draw("Ef","pdgf==-211","");
  if(gst_1) gst_1->Draw("Ef","pdgf==-211","perrsame");
  c->cd();
  ls->Clear();
  ls->SetHeader("Final state pi-'s 4-momentum");
  ls->Draw();
  c->Update();

  if(show_mult_per_proc) {

     //
     // similarly but for QEL events only
     //

     //------ number of final state p /QEL
     ps->NewPage();
     if(gst_1) gst_1->Draw("nfp","qel","");
     gst_0->Draw("nfp","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state protons / QEL only");
     ls->Draw();
     c->Update();

     //------ number of final state n /QEL
     ps->NewPage();
     if(gst_1) gst_1->Draw("nfn","qel","");
     gst_0->Draw("nfn","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state neutrons / QEL only");
     ls->Draw();
     c->Update();

     //------ number of final state pi+ /QEL
     ps->NewPage();
     gst_0->Draw("nfpip","qel","");
     if(gst_1) gst_1->Draw("nfpip","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi+ / QEL only");
     ls->Draw();
     c->Update();

     //------ number of final state pi- /QEL
     ps->NewPage();
     gst_0->Draw("nfpim","qel","");
     if(gst_1) gst_1->Draw("nfpim","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi- / QEL only");
     ls->Draw();
     c->Update();

     //------ number of final state pi0 /QEL
     ps->NewPage();
     gst_0->Draw("nfpi0","qel","");
     if(gst_1) gst_1->Draw("nfpi0","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi0 / QEL only");
     ls->Draw();
     c->Update();

     //------ number of final state K+ /QEL
     ps->NewPage();
     gst_0->Draw("nfkp","qel","");
     if(gst_1) gst_1->Draw("nfkp","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K+ / QEL only");
     ls->Draw();
     c->Update();

     //------ number of final state K- /QEL
     ps->NewPage();
     gst_0->Draw("nfkm","qel","");
     if(gst_1) gst_1->Draw("nfkm","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K- / QEL only");
     ls->Draw();
     c->Update();

     //------ number of final state K0 /QEL
     ps->NewPage();
     gst_0->Draw("nfk0","qel","");
     if(gst_1) gst_1->Draw("nfk0","qel","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K0 / QEL only");
     ls->Draw();
     c->Update();

     //------ momentum of final state p /QEL
     ps->NewPage();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","qel&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pxf","qel&&pdgf==2212","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","qel&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pyf","qel&&pdgf==2212","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","qel&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pzf","qel&&pdgf==2212","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","qel&&pdgf==2212","");
     if(gst_1) gst_1->Draw("Ef","qel&&pdgf==2212","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state protons 4-momentum / QEL only");
     ls->Draw();
     c->Update();

     //------ momentum of final state n /QEL
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","qel&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pxf","qel&&pdgf==2112","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","qel&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pyf","qel&&pdgf==2112","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","qel&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pzf","qel&&pdgf==2112","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","qel&&pdgf==2112","");
     if(gst_1) gst_1->Draw("Ef","qel&&pdgf==2112","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state neutrons 4-momentum / QEL only");
     ls->Draw();
     c->Update();

     //------ momentum of final state pi0 /QEL
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","qel&&pdgf==111","");
     if(gst_1) gst_1->Draw("pxf","qel&&pdgf==111","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","qel&&pdgf==111","");
     if(gst_1) gst_1->Draw("pyf","qel&&pdgf==111","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","qel&&pdgf==111","");
     if(gst_1) gst_1->Draw("pzf","qel&&pdgf==111","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","qel&&pdgf==111","");
     if(gst_1) gst_1->Draw("Ef","qel&&pdgf==111","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi0's 4-momentum / QEL only");
     ls->Draw();
     c->Update();

     //------ momentum of final state pi+ /QEL
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","qel&&pdgf==211","");
     if(gst_1) gst_1->Draw("pxf","qel&&pdgf==211","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","qel&&pdgf==211","");
     if(gst_1) gst_1->Draw("pyf","qel&&pdgf==211","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","qel&&pdgf==211","");
     if(gst_1) gst_1->Draw("pzf","qel&&pdgf==211","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","qel&&pdgf==211","");
     if(gst_1) gst_1->Draw("Ef","qel&&pdgf==211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi+'s 4-momentum / QEL only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state pi+ /QEL
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","qel&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pxf","qel&&pdgf==-211","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","qel&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pyf","qel&&pdgf==-211","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","qel&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pzf","qel&&pdgf==-211","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","qel&&pdgf==-211","");
     if(gst_1) gst_1->Draw("Ef","qel&&pdgf==-211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi-'s 4-momentum/ QEL only");
     ls->Draw();
     c->Update();
     
     
     //
     // similarly but for RES events only
     //
     
     //------ number of final state p /RES
     ps->NewPage();
               gst_0->Draw("nfp","res","");
     if(gst_1) gst_1->Draw("nfp","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state protons / RES only");
     ls->Draw();
     c->Update();
     
     //------ number of final state n /RES
     ps->NewPage();
               gst_0->Draw("nfn","res","");
     if(gst_1) gst_1->Draw("nfn","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state neutrons / RES only");
     ls->Draw();
     c->Update();
     
     //------ number of final state pi+ /RES
     ps->NewPage();
               gst_0->Draw("nfpip","res","");
     if(gst_1) gst_1->Draw("nfpip","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi+ / RES only");
     ls->Draw();
     c->Update();
     
     //------ number of final state pi- /RES
     ps->NewPage();
               gst_0->Draw("nfpim","res","");
     if(gst_1) gst_1->Draw("nfpim","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi- / RES only");
     ls->Draw();
     c->Update();

     //------ number of final state pi0 /RES
     ps->NewPage();
               gst_0->Draw("nfpi0","res","");
     if(gst_1) gst_1->Draw("nfpi0","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi0 / RES only");
     ls->Draw();
     c->Update();
     
     //------ number of final state K+ /RES
     ps->NewPage();
     gst_0->Draw("nfkp","res","");
     if(gst_1) gst_1->Draw("nfkp","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K+ / RES only");
     ls->Draw();
     c->Update();
     
     //------ number of final state K- /RES
     ps->NewPage();
     gst_0->Draw("nfkm","res","");
     if(gst_1) gst_1->Draw("nfkm","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K- / RES only");
     ls->Draw();
     c->Update();
     
     //------ number of final state K0 /RES
     ps->NewPage();
     gst_0->Draw("nfk0","res","");
     if(gst_1) gst_1->Draw("nfk0","res","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K0 / RES only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state p /QEL
     ps->NewPage();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","res&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pxf","res&&pdgf==2212","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","res&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pyf","res&&pdgf==2212","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","res&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pzf","res&&pdgf==2212","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","res&&pdgf==2212","");
     if(gst_1) gst_1->Draw("Ef","res&&pdgf==2212","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state protons 4-momentum / RES only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state n /RES
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","res&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pxf","res&&pdgf==2112","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","res&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pyf","res&&pdgf==2112","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","res&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pzf","res&&pdgf==2112","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","res&&pdgf==2112","");
     if(gst_1) gst_1->Draw("Ef","res&&pdgf==2112","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state neutrons 4-momentum / RES only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state pi0 /RES
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","res&&pdgf==111","");
     if(gst_1) gst_1->Draw("pxf","res&&pdgf==111","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","res&&pdgf==111","");
     if(gst_1) gst_1->Draw("pyf","res&&pdgf==111","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","res&&pdgf==111","");
     if(gst_1) gst_1->Draw("pzf","res&&pdgf==111","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","res&&pdgf==111","");
     if(gst_1) gst_1->Draw("Ef","res&&pdgf==111","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi0's 4-momentum / RES only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state pi+ /RES
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","res&&pdgf==211","");
     if(gst_1) gst_1->Draw("pxf","res&&pdgf==211","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","res&&pdgf==211","");
     if(gst_1) gst_1->Draw("pyf","res&&pdgf==211","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","res&&pdgf==211","");
     if(gst_1) gst_1->Draw("pzf","res&&pdgf==211","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","res&&pdgf==211","");
     if(gst_1) gst_1->Draw("Ef","res&&pdgf==211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi+'s 4-momentum / RES only");
     ls->Draw();
     c->Update();

     //------ momentum of final state pi+ /RES
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","res&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pxf","res&&pdgf==-211","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","res&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pyf","res&&pdgf==-211","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","res&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pzf","res&&pdgf==-211","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","res&&pdgf==-211","");
     if(gst_1) gst_1->Draw("Ef","res&&pdgf==-211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi-'s 4-momentum/ RES only");
     ls->Draw();
     c->Update();
     
     //
     // similarly but for DIS events only
     //
     
     //------ number of final state p /DIS
     ps->NewPage();
     gst_0->Draw("nfp","dis","");
     if(gst_1) gst_1->Draw("nfp","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state protons / DIS only");
     ls->Draw();
     c->Update();
     
     //------ number of final state n /DIS
     ps->NewPage();
     gst_0->Draw("nfn","dis","");
     if(gst_1) gst_1->Draw("nfn","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state neutrons / DIS only");
     ls->Draw();
     c->Update();
     
     //------ number of final state pi+ /DIS
     ps->NewPage();
     gst_0->Draw("nfpip","dis","");
     if(gst_1) gst_1->Draw("nfpip","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi+ / DIS only");
     ls->Draw();
     c->Update();
     
     //------ number of final state pi- /DIS
     ps->NewPage();
     gst_0->Draw("nfpim","dis","");
     if(gst_1) gst_1->Draw("nfpim","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi- / DIS only");
     ls->Draw();
     c->Update();
     
     //------ number of final state pi0 /DIS
     ps->NewPage();
     gst_0->Draw("nfpi0","dis","");
     if(gst_1) gst_1->Draw("nfpi0","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state pi0 / DIS only");
     ls->Draw();
     c->Update();
     
     //------ number of final state K+ /DIS
     ps->NewPage();
     gst_0->Draw("nfkp","dis","");
     if(gst_1) gst_1->Draw("nfkp","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K+ / DIS only");
     ls->Draw();
     c->Update();
     
     //------ number of final state K- /DIS
     ps->NewPage();
     gst_0->Draw("nfkm","dis","");
     if(gst_1) gst_1->Draw("nfkm","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K- / DIS only");
     ls->Draw();
     c->Update();
     
     //------ number of final state K0 /DIS
     ps->NewPage();
     gst_0->Draw("nfk0","dis","");
     if(gst_1) gst_1->Draw("nfk0","dis","perrsame");
     ls->Clear();
     ls->SetHeader("Number of final state K0 / DIS only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state p /DIS
     ps->NewPage();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","dis&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pxf","dis&&pdgf==2212","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","dis&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pyf","dis&&pdgf==2212","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","dis&&pdgf==2212","");
     if(gst_1) gst_1->Draw("pzf","dis&&pdgf==2212","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","dis&&pdgf==2212","");
     if(gst_1) gst_1->Draw("Ef","dis&&pdgf==2212","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state protons 4-momentum / DIS only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state n /DIS
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","dis&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pxf","dis&&pdgf==2112","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","dis&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pyf","dis&&pdgf==2112","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","dis&&pdgf==2112","");
     if(gst_1) gst_1->Draw("pzf","dis&&pdgf==2112","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","dis&&pdgf==2112","");
     if(gst_1) gst_1->Draw("Ef","dis&&pdgf==2112","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state neutrons 4-momentum / DIS only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state pi0 /DIS
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","dis&&pdgf==111","");
     if(gst_1) gst_1->Draw("pxf","dis&&pdgf==111","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","dis&&pdgf==111","");
     if(gst_1) gst_1->Draw("pyf","dis&&pdgf==111","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","dis&&pdgf==111","");
     if(gst_1) gst_1->Draw("pzf","dis&&pdgf==111","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","dis&&pdgf==111","");
     if(gst_1) gst_1->Draw("Ef","dis&&pdgf==111","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi0's 4-momentum / DIS only");
     ls->Draw();
     c->Update();

     //------ momentum of final state pi+ /DIS
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","dis&&pdgf==211","");
     if(gst_1) gst_1->Draw("pxf","dis&&pdgf==211","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","dis&&pdgf==211","");
     if(gst_1) gst_1->Draw("pyf","dis&&pdgf==211","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","dis&&pdgf==211","");
     if(gst_1) gst_1->Draw("pzf","dis&&pdgf==211","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","dis&&pdgf==211","");
     if(gst_1) gst_1->Draw("Ef","dis&&pdgf==211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi+'s 4-momentum / DIS only");
     ls->Draw();
     c->Update();
     
     //------ momentum of final state pi+ /DIS
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxf","dis&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pxf","dis&&pdgf==-211","perrsame");
     c->cd(2);
     gst_0->Draw("pyf","dis&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pyf","dis&&pdgf==-211","perrsame");
     c->cd(3);
     gst_0->Draw("pzf","dis&&pdgf==-211","");
     if(gst_1) gst_1->Draw("pzf","dis&&pdgf==-211","perrsame");
     c->cd(4);
     gst_0->Draw("Ef","dis&&pdgf==-211","");
     if(gst_1) gst_1->Draw("Ef","dis&&pdgf==-211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Final state pi-'s 4-momentum/ DIS only");
     ls->Draw();
     c->Update();
     
  } // per-proc

  //
  // SECTION: Primary Hadronic System Multiplicities & 4P
  //
  if(show_primary_hadsyst) {
     
     ps->NewPage();
     c->Clear();
     c->Range(0,0,100,100);
     TPavesText hdrihad(10,40,90,70,3,"tr");
     hdrihad.AddText("Parimary Hadronic System");
     hdrihad.AddText("Multiplicities and 4-Momenta");
     hdrihad.AddText(" ");
     hdrihad.AddText(" ");
     hdrihad.AddText(" ");
     hdrihad.AddText(" ");
     hdrihad.AddText("Note:");
     hdrihad.AddText("For nuclear targets these plots show the hadronic system");
     hdrihad.AddText("BEFORE any intranuclear hadron transport / rescattering");
     hdrihad.Draw();
     c->Update();
     
     //------ number of prim p
     ps->NewPage();
     gst_0->Draw("nip","","");
     if(gst_1) gst_1->Draw("nfp","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of protons");
     ls->Draw();
     c->Update();
     
     //------ number of prim n
     ps->NewPage();
     gst_0->Draw("nin","","");
     if(gst_1) gst_1->Draw("nfn","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of neutrons");
     ls->Draw();
     c->Update();
     
     //------ number of prim pi+
     ps->NewPage();
     gst_0->Draw("nipip","","");
     if(gst_1) gst_1->Draw("nfpip","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of pi+");
     ls->Draw();
     c->Update();
     
     //------ number of prim pi-
     ps->NewPage();
     gst_0->Draw("nipim","","");
     if(gst_1) gst_1->Draw("nfpim","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of pi-");
     ls->Draw();
     c->Update();
     
     //------ number of prim pi0
     ps->NewPage();
     gst_0->Draw("nipi0","","");
     if(gst_1) gst_1->Draw("nfpi0","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of pi0");
     ls->Draw();
     c->Update();
     
     //------ number of prim K+
     ps->NewPage();
     gst_0->Draw("nikp","","");
     if(gst_1) gst_1->Draw("nikp","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of K+");
     ls->Draw();
     c->Update();
     
     //------ number of prim K-
     ps->NewPage();
     gst_0->Draw("nikm","","");
     if(gst_1) gst_1->Draw("nikm","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of K-");
     ls->Draw();
     c->Update();
     
     //------ number of prim K0
     ps->NewPage();
     gst_0->Draw("nik0","","");
     if(gst_1) gst_1->Draw("nik0","","perrsame");
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: Number of K0");
     ls->Draw();
     c->Update();
     
     //------ momentum of prim, p
     ps->NewPage();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxi","pdgi==2212","");
     if(gst_1) gst_1->Draw("pxi","pdgi==2212","perrsame");
     c->cd(2);
     gst_0->Draw("pyi","pdgi==2212","");
     if(gst_1) gst_1->Draw("pyi","pdgi==2212","perrsame");
     c->cd(3);
     gst_0->Draw("pzi","pdgi==2212","");
     if(gst_1) gst_1->Draw("pzi","pdgi==2212","perrsame");
     c->cd(4);
     gst_0->Draw("Ei","pdgi==2212","");
     if(gst_1) gst_1->Draw("Ei","pdgi==2212","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: proton 4-momentum");
     ls->Draw();
     c->Update();
     
     //------ momentum of prim. n
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxi","pdgi==2112","");
     if(gst_1) gst_1->Draw("pxi","pdgi==2112","perrsame");
     c->cd(2);
     gst_0->Draw("pyi","pdgi==2112","");
     if(gst_1) gst_1->Draw("pyi","pdgi==2112","perrsame");
     c->cd(3);
     gst_0->Draw("pzi","pdgi==2112","");
     if(gst_1) gst_1->Draw("pzi","pdgi==2112","perrsame");
     c->cd(4);
     gst_0->Draw("Ei","pdgi==2112","");
     if(gst_1) gst_1->Draw("Ei","pdgi==2112","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: neutron 4-momentum");
     ls->Draw();
     c->Update();
     
     //------ momentum of prim. pi0
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxi","pdgi==111","");
     if(gst_1) gst_1->Draw("pxi","pdgi==111","perrsame");
     c->cd(2);
     gst_0->Draw("pyi","pdgi==111","");
     if(gst_1) gst_1->Draw("pyi","pdgi==111","perrsame");
     c->cd(3);
     gst_0->Draw("pzi","pdgi==111","");
     if(gst_1) gst_1->Draw("pzi","pdgi==111","perrsame");
     c->cd(4);
     gst_0->Draw("Ei","pdgi==111","");
     if(gst_1) gst_1->Draw("Ei","pdgi==111","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Primary Hadronic System: pi0's 4-momentum");
     ls->Draw();
     c->Update();

     //------ momentum of prim pi+
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxi","pdgi==211","");
     if(gst_1) gst_1->Draw("pxi","pdgi==211","perrsame");
     c->cd(2);
     gst_0->Draw("pyi","pdgi==211","");
     if(gst_1) gst_1->Draw("pyi","pdgi==211","perrsame");
     c->cd(3);
     gst_0->Draw("pzi","pdgi==211","");
     if(gst_1) gst_1->Draw("pzi","pdgi==211","perrsame");
     c->cd(4);
     gst_0->Draw("Ei","pdgi==211","");
     if(gst_1) gst_1->Draw("Ei","pdgi==211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Primary Hadronic System:pi+'s 4-momentum");
     ls->Draw();
     c->Update();
     
     //------ momentum of prim. pi+
     ps->NewPage();
     c->Clear();
     c->Divide(2,2);
     c->cd(1);
     gst_0->Draw("pxi","pdgi==-211","");
     if(gst_1) gst_1->Draw("pxi","pdgi==-211","perrsame");
     c->cd(2);
     gst_0->Draw("pyi","pdgi==-211","");
     if(gst_1) gst_1->Draw("pyi","pdgi==-211","perrsame");
     c->cd(3);
     gst_0->Draw("pzi","pdgi==-211","");
     if(gst_1) gst_1->Draw("pzi","pdgi==-211","perrsame");
     c->cd(4);
     gst_0->Draw("Ei","pdgi==-211","");
     if(gst_1) gst_1->Draw("Ei","pdgi==-211","perrsame");
     c->cd();
     ls->Clear();
     ls->SetHeader("Primary Hadronic System:  pi-'s 4-momentum");
     ls->Draw();
     c->Update();
  }//show?
       
  ps->Close();

  fin_0->Close();
  delete fin_0;
  if(fin_1) {
     fin_1->Close();
     delete fin_1;
  }
}
//_________________________________________________________________________________
string OutputFileName(string inpname, int mod)
{
// Builds the output filename based on the name of the input filename
// Perfors the following conversion: name.root -> name.gmctest.root

  assert(mod==0 || mod==1);
  unsigned int L = inpname.length();

  // if the last 4 characters are "root" (ROOT file extension) then
  // remove them
  if(inpname.substr(L-4, L).find("root") != string::npos) {
    inpname.erase(L-4, L);
  }

  // remove ghep.
  size_t pos = inpname.find("ghep.");
  if(pos != string::npos) {
    inpname.erase(pos, pos+4);
  }

  ostringstream name;

  if(mod==0)
  	name << inpname << "gst.root";
  else if(mod==1)
  	name << inpname << "gmctest.ps";

  return gSystem->BaseName(name.str().c_str());
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gmctest", pNOTICE) << "*** Parsing commad line arguments";

  //mode:
  try {
    LOG("gmctest", pINFO) << "Reading gmctest mode";
    gOptMode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'m');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctest", pFATAL) << "Unspecified gmctest mode";
      PrintSyntax();
      exit(1);
    }
  }

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

  //get GENIE event sample
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

  //get another (reference) GENIE event sample 
  try {
    LOG("gmctest", pINFO) << "Reading filename of event sample template";
    gOptInpFileRef = utils::clap::CmdLineArgAsString(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmctst", pNOTICE) << "Unspecified 'reference' event sample";
    }
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmctest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gmctest -m mode -f sample.root [-n nev] [-r reference_sample.root]\n";
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

