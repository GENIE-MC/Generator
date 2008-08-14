//____________________________________________________________________________
/*!

\program gntpc

\brief   Converts a GENIE neutrino events (from GHEP GENIE ROOT Trees) to a 
         variety of textual formats (typically used in legacy systems), in XML
         format or in summary ROOT ntuples.

         Syntax:
           gntpc -i input_filename [-o output_filename] -f format [-n nev]

         Options :
           [] denotes an optional argument
           -n number of events to convert
           -f specifies the output file format. 
	      ** 'Summary' ntuple formats **
                    0 : The 'definite' GENIE summary tree format (gst). 
                        The program will analyze the input event tree and save
                        usefull/detailed summary information on a flat ntuple 
                        that can be trivially used in bare-ROOT sessions.
	      ** T2K/GENIE formats **
   		    1 : NUANCE-style tracker text-based format 
   		   11 : Same as format 1 but with extra tweaks required by
                        skdetsim, the SuperK detector MC:
                        - K0, \bar{K0} -> KO_{long}, K0_{short}
                        - emulating 'NEUT' reaction codes
		    2 : A slightly tweaked NUANCE-style tracker text-based 
                        format - A fast & dirty way for getting GENIE event 
                        samples into the T2K (nd280/2km/SuperK) Monte Carlo.
   	            3 : A standardized bare-ROOT event-tree for getting GENIE
                        neutrino  & pass-through JPARC flux info into the
                        T2K (nd280/2km/SuperK) Monte Carlo
	      ** Generic GENIE XML / tabular or bare-ROOT formats **
                  100 : GENIE XML format 
	      ** GENIE test / cross-generator comparisons **
	          901 : NEUGEN-style text-based format for hadronization studies
	          902 : INTRANUKE summary ntuple for intranuclear-rescattering studies
           -o specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
                0 -> *.gst.root
                1 -> *.gtrac0.dat
               11 -> *.gtrac0.dat
                2 -> *.gtrac.dat
                3 -> *.gtrac.root
              100 -> *.gxml 
              901 -> *.ghad.dat
              902 -> *.ginuke.root
		
\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created September 23, 2005

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TBits.h>
#include <TObjString.h>
#include <TMath.h>

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#include "FluxDrivers/GJPARCNuFlux.h"
#endif

//define __GHAD_NTP__

using std::string;
using std::ostringstream;
using std::ofstream;
using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;
using std::vector;

using namespace genie;
using namespace genie::constants;

//func prototypes
void   ConvertToGST            (void);
void   ConvertToGT2KTracker    (void);
void   ConvertToGT2KRooTracker (void);
void   ConvertToGXML           (void);
void   ConvertToGHad           (void);
void   ConvertToGINuke         (void);
void   GetCommandLineArgs      (int argc, char ** argv);
void   PrintSyntax             (void);
string DefaultOutputFile       (void);
bool   CheckRootFilename       (string filename);

//input options (from command line arguments):
string   gOptInpFileName;
string   gOptOutFileName;
int      gOptOutFileFormat;
Long64_t gOptN;

//consts
const int kNPmax = 100;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);

  //-- call the appropriate conversion function
  switch(gOptOutFileFormat) {
   case (0)  :
	ConvertToGST();        
	break;  
   case (1)  :  
   case (2)  :  
   case (11) :  
	ConvertToGT2KTracker();        
	break;
   case (3) :  
	ConvertToGT2KRooTracker(); 
	break;
   case (100) :  
	ConvertToGXML();         
	break;
   case (901) :  
	ConvertToGHad();         
	break;
   case (902) :  
	ConvertToGINuke();         
	break;
   default:
     LOG("gntpc", pFATAL)
          << "Invalid output format [" << gOptOutFileFormat << "]";
     PrintSyntax();
     exit(3);
  }
  return 0;
}
//___________________________________________________________________
//    **** GENIE GHEP EVENT TREE -> GENIE SUMMARY NTUPLE ****
//___________________________________________________________________
void ConvertToGST(void)
{
  //-- define branch variables
  //
  int    brIev         = 0;      // Event number 
  int    brNeutrino    = 0;      // Neutrino pdg code
  int    brTarget      = 0;      // Nuclear target pdg code (10LZZZAAAI)
  int    brTargetZ     = 0;      // Nuclear target Z (extracted from pdg code above)
  int    brTargetA     = 0;      // Nuclear target A (extracted from pdg code above)
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
  int    brNfP         = 0;      // Nu. of final state p's + \bar{p}'s (after intranuclear rescattering)
  int    brNfN         = 0;      // Nu. of final state n's + \bar{n}'s
  int    brNfPip       = 0;      // Nu. of final state pi+'s
  int    brNfPim       = 0;      // Nu. of final state pi-'s
  int    brNfPi0       = 0;      // Nu. of final state pi0's (
  int    brNfKp        = 0;      // Nu. of final state K+'s
  int    brNfKm        = 0;      // Nu. of final state K-'s
  int    brNfK0        = 0;      // Nu. of final state K0's + \bar{K0}'s
  int    brNfEM        = 0;      // Nu. of final state gammas and e-/e+ (excluding pi0 decay products)
  int    brNfOther     = 0;      // Nu. of heavier final state hadrons (D+/-,D0,Ds+/-,Lamda,Sigma,Lamda_c,Sigma_c,...)
  int    brNiP         = 0;      // Nu. of 'primary' (: before intranuclear rescattering) p's + \bar{p}'s  
  int    brNiN         = 0;      // Nu. of 'primary' n's + \bar{n}'s  
  int    brNiPip       = 0;      // Nu. of 'primary' pi+'s 
  int    brNiPim       = 0;      // Nu. of 'primary' pi-'s 
  int    brNiPi0       = 0;      // Nu. of 'primary' pi0's 
  int    brNiKp        = 0;      // Nu. of 'primary' K+'s  
  int    brNiKm        = 0;      // Nu. of 'primary' K-'s  
  int    brNiK0        = 0;      // Nu. of 'primary' K0's + \bar{K0}'s 
  int    brNiEM        = 0;      // Nu. of 'primary' gammas and e-/e+ (eg from resonance decays)
  int    brNiOther     = 0;      // Nu. of 'primary' hadron shower particles
  int    brNf          = 0;      // Nu. of final state particles in hadronic system
  int    brPdgf[kNPmax];         // Pdg code of i^th final state particle in hadronic system
  double brEf  [kNPmax];         // Energy   of i^th final state particle in hadronic system
  double brPxf [kNPmax];         // Px       of i^th final state particle in hadronic system
  double brPyf [kNPmax];         // Py       of i^th final state particle in hadronic system
  double brPzf [kNPmax];         // Pz       of i^th final state particle in hadronic system
  int    brNi          = 0;      // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
  int    brPdgi[kNPmax];         // Pdg code of i^th particle in 'primary' hadronic system 
  double brEi  [kNPmax];         // Energy   of i^th particle in 'primary' hadronic system 
  double brPxi [kNPmax];         // Px       of i^th particle in 'primary' hadronic system 
  double brPyi [kNPmax];         // Py       of i^th particle in 'primary' hadronic system 
  double brPzi [kNPmax];         // Pz       of i^th particle in 'primary' hadronic system 
  double brVtxX;                 // Vertex x in detector coord system (SI)
  double brVtxY;                 // Vertex y in detector coord system (SI)
  double brVtxZ;                 // Vertex z in detector coord system (SI)
  double brVtxT;                 // Vertex t in detector coord system (SI)

  //-- open output file & create output summary tree & create the tree branches
  //
  LOG("gntpc", pNOTICE) 
       << "*** Saving summary tree to: " << gOptOutFileName;
  TFile fout(gOptOutFileName.c_str(),"recreate");

  TTree * s_tree = new TTree("gst","GENIE Summary Event Tree");

  //-- create tree branches
  //
  s_tree->Branch("iev",       &brIev,       "iev/I"       );
  s_tree->Branch("neu",	      &brNeutrino,  "neu/I"	  );
  s_tree->Branch("tgt",       &brTarget,    "tgt/I"	  );
  s_tree->Branch("Z",         &brTargetZ,   "Z/I"	  );
  s_tree->Branch("A",         &brTargetA,   "A/I"	  );
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
  s_tree->Branch("vtxx",      &brVtxX,	    "vtxx/D"      );
  s_tree->Branch("vtxy",      &brVtxY,	    "vtxy/D"      );
  s_tree->Branch("vtxz",      &brVtxZ,	    "vtxz/D"      );
  s_tree->Branch("vtxt",      &brVtxT,	    "vtxt/D"      );

  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           er_tree = 0;
  NtpMCTreeHeader * thdr    = 0;
  er_tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr    = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );
  if (!er_tree) {
    LOG("gntpc", pERROR) << "Null input ER tree";
    return;
  }
  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- get the mc record
  NtpMCEventRecord * mcrec = 0;
  er_tree->SetBranchAddress("gmcrec", &mcrec);

  if (!mcrec) {
    LOG("gntpc", pERROR) << "Null MC record";
    return;
  }
  
  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       er_tree->GetEntries() : TMath::Min( er_tree->GetEntries(), gOptN );
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  TLorentzVector pdummy(0,0,0,0);

  for(Long64_t iev = 0; iev < nmax; iev++) {

    er_tree->GetEntry(iev);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    // go further only if the event is physical
    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) {
      LOG("gntpc", pINFO) << "Skipping unphysical event";
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

    int tgtZ = 0;
    int tgtA = 0;
    if(pdg::IsIon(target->Pdg())) {
       tgtZ = pdg::IonPdgCodeToZ(target->Pdg());
       tgtA = pdg::IonPdgCodeToA(target->Pdg());
    } 
    if(target->Pdg() == kPdgProton ) { tgtZ = 1; tgtA = 1; }    
    if(target->Pdg() == kPdgNeutron) { tgtZ = 0; tgtA = 1; }    
  
    //summary info
    const Interaction * interaction = event.Summary();
    const InitialState & init_state = interaction->InitState();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    const Kinematics &   kine       = interaction->Kine();
    const XclsTag &      xcls       = interaction->ExclTag();
    const Target &       tgt        = init_state.Tgt();

    //vertex in detector coord system
    TLorentzVector * vtx = event.Vertex();

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

    LOG("gntpc", pDEBUG) 
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

    LOG("gntpc", pDEBUG) 
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

    LOG("gntpc", pDEBUG) << "Extracting final state hadronic system";

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

    LOG("gntpc", pDEBUG) << "Extracting primary hadronic system";

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
    brIev        = (int) iev;      
    brNeutrino   = neutrino->Pdg();      
    brTarget     = target->Pdg(); 
    brTargetZ    = tgtZ;
    brTargetA    = tgtA;   
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

      LOG("gntpc", pINFO) 
        << "Counting in primary hadronic system: idx = " << prim_had_syst[j]
        << " -> " << p->Name();
    }

    LOG("gntpc", pINFO) 
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

      LOG("gntpc", pINFO) 
        << "Counting in f/s system from hadronic vtx: idx = " << final_had_syst[j]
        << " -> " << p->Name();
    }

    LOG("gntpc", pINFO) 
     << "N(p):"             << brNfP
     << ", N(n):"           << brNfN
     << ", N(pi+):"         << brNfPip
     << ", N(pi-):"         << brNfPim
     << ", N(pi0):"         << brNfPi0
     << ", N(K+,K-,K0):"    << brNfKp+brNfKm+brNfK0
     << ", N(gamma,e-,e+):" << brNfEM
     << ", N(etc):"         << brNfOther << "\n";

    brVtxX = vtx->X();   
    brVtxY = vtx->Y();   
    brVtxZ = vtx->Z();   
    brVtxT = vtx->T();

    s_tree->Fill();

    mcrec->Clear();

  } // event loop

  fin.Close();

  fout.Write();
  fout.Close();
}
//___________________________________________________________________
//    **** GENIE GHEP EVENT TREE -> GENIE XML EVENT FILE ****
//___________________________________________________________________
void ConvertToGXML(void)
{
  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- get mc record
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- open the output stream
  ofstream output(gOptOutFileName.c_str(), ios::out);

  //-- add required header
  output << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  output << endl << endl;
  output << "<!-- generated by GENIE gntpc utility -->";   
  output << endl << endl;
  output << "<genie_event_list version=\"1.00\">" << endl;

  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       tree->GetEntries() : TMath::Min(tree->GetEntries(), gOptN);
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  //-- event loop
  for(Long64_t iev = 0; iev < nmax; iev++) {
    tree->GetEntry(iev);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    //
    // convert the current event
    //

    output << endl << endl;
    output << "  <!-- GENIE GHEP event -->" << endl;
    output << "  <ghep np=\"" << event.GetEntries() 
           << "\" unphysical=\"" 
           << (event.IsUnphysical() ? "true" : "false") << "\">" << endl;
    output << setiosflags(ios::scientific);

    // write-out the event-wide properties
    output << "   ";
    output << "  <!-- event weight   -->";
    output << " <wgt> " << event.Weight()   << " </wgt>";
    output << endl;
    output << "   ";
    output << "  <!-- cross sections -->";
    output << " <xsec_evnt> " << event.XSec()     << " </xsec_evnt>";
    output << " <xsec_kine> " << event.DiffXSec() << " </xsec_kine>";
    output << endl;
    output << "   ";
    output << "  <!-- event vertex   -->";
    output << " <vx> " << event.Vertex()->X() << " </vx>";
    output << " <vy> " << event.Vertex()->Y() << " </vy>";
    output << " <vz> " << event.Vertex()->Z() << " </vz>";
    output << " <vt> " << event.Vertex()->T() << " </vt>";
    output << endl;

    //  write-out the generated particle list
    output << "     <!-- particle list  -->" << endl;
    unsigned int i=0;
    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
      string type = "U";
      if(p->IsFake()) type = "F";
      else {
        if(p->IsParticle()) type = "P";
        if(p->IsNucleus() ) type = "N";
      }
      output << "     <p idx=\"" << i << "\" type=\"" << type << "\">" << endl;
      output << "        ";
      output << " <pdg> " << p->Pdg()       << " </pdg>";
      output << " <ist> " << p->Status()    << " </ist>";
      output << endl;
      output << "        ";
      output << " <mother>   "  
             << " <fst> " << setfill(' ') << setw(3) << p->FirstMother() << " </fst> "
             << " <lst> " << setfill(' ') << setw(3) << p->LastMother()  << " </lst> "
             << " </mother>";
      output << endl;
      output << "        ";
      output << " <daughter> "  
             << " <fst> " << setfill(' ') << setw(3) << p->FirstDaughter() << " </fst> "
             << " <lst> " << setfill(' ') << setw(3) << p->LastDaughter()  << " </lst> "
             << " </daughter>";
      output << endl;
      output << "        ";
      output << " <px> " << setfill(' ') << setw(20) << p->Px() << " </px>";
      output << " <py> " << setfill(' ') << setw(20) << p->Py() << " </py>";
      output << " <pz> " << setfill(' ') << setw(20) << p->Pz() << " </pz>";
      output << " <E>  " << setfill(' ') << setw(20) << p->E()  << " </E> ";
      output << endl;
      output << "        ";
      output << " <x>  " << setfill(' ') << setw(20) << p->Vx() << " </x> ";
      output << " <y>  " << setfill(' ') << setw(20) << p->Vy() << " </y> ";
      output << " <z>  " << setfill(' ') << setw(20) << p->Vz() << " </z> ";
      output << " <t>  " << setfill(' ') << setw(20) << p->Vt() << " </t> ";
      output << endl;

      if(p->PolzIsSet()) {
        output << "        ";
        output << " <ppolar> " << p->PolzPolarAngle()   << " </ppolar>";
        output << " <pazmth> " << p->PolzAzimuthAngle() << " </pazmth>";
        output << endl;
      }

      output << "     </p>" << endl;
      i++;
    }
    output << "  </ghep>" << endl;

    mcrec->Clear();
  } // event loop

  //-- add required footer
  output << endl << endl;
  output << "<genie_event_list version=\"1.00\">";

  output.close();
  fin.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//___________________________________________________________________
// **** GENIE GHEP EVENT TREE -> NUANCE-STYLE TRACKER TEXT FILE ****
//___________________________________________________________________
void ConvertToGT2KTracker(void)
{
  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- get mc record
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- open the output stream
  ofstream output(gOptOutFileName.c_str(), ios::out);

  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       tree->GetEntries() : TMath::Min(tree->GetEntries(), gOptN);
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  //-- event loop
  for(Long64_t iev = 0; iev < nmax; iev++) {
    tree->GetEntry(iev);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);
    Interaction * interaction = event.Summary();

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    GHepParticle * p = 0;
    TIter event_iter(&event);

    //
    // convert the current event
    //

    // add tracker begin tag
    output << "$ begin" << endl;

    // add 'NUANCE'-like event type
    //
    if(gOptOutFileFormat==1) {
    	int evtype = 0;
        /*
        conversion not tested yet
        const ProcessInfo &  proc = interaction->ProcInfo();
        const InitialState & init = interaction->InitState();
        if      (proc.IsQuasiElastic()   && proc.IsWeakCC()) evtype =  1;
        else if (proc.IsQuasiElastic()   && proc.IsWeakNC()) evtype =  2;
        else if (proc.IsDeepInelastic()  && proc.IsWeakCC()) evtype = 91;        
        else if (proc.IsDeepInelastic()  && proc.IsWeakNC()) evtype = 92;        
        else if (proc.IsCoherentPiProd() && proc.IsWeakNC()) evtype = 96;        
        else if (proc.IsCoherentPiProd() && proc.IsWeakCC()) evtype = 97;        
        else if (proc.IsNuElectronElastic())                 evtype = 98;
        else if (proc.IsInverseMuDecay())                    evtype = 99;
        else if (proc.IsResonant()) {
           int nn=0, np=0, npi0=0, npip=0, npim=0;
           bool nuclear_target = init.Tgt().IsNucleus();
           GHepStatus_t matched_ist = (nuclear_target) ? 
                     kIStHadronInTheNucleus : kIStStableFinalState;
           while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) 
           {
               GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();  
               if(ghep_ist != matched_ist) continue;
               int ghep_pdgc = p->Pdg();
               if(ghep_pdgc == kPdgProton ) np++;
               if(ghep_pdgc == kPdgNeutron) nn++;
               if(ghep_pdgc == kPdgPi0)     npi0++;
               if(ghep_pdgc == kPdgPiP)     npip++;
               if(ghep_pdgc == kPdgPiM)     npim++;
           }
           if(proc.IsWeakCC() && init.IsNuP()) {
             // v p -> l- p pi+ 
             if(np==1 && nn==0 && npip==1 && npi0==0 && npim==0) evtype = 3;  
           }
           if(proc.IsWeakCC() && init.IsNuN()) {
             // v n -> l- p pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 4;  
             // v n -> l- n pi+ 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 5;  
           }
           if(proc.IsWeakNC() && init.IsNuP()) {
             // v p -> v p pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 6;  
             // v p -> v n pi+ 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 7;  
           }
           if(proc.IsWeakNC() && init.IsNuN()) {
             // v n -> v n pi0 
             if(np==0 && nn==1 && npip==0 && npi0==1 && npim==0) evtype = 8;  
             // v n -> v p pi- 
             if(np==1 && nn==0 && npip==0 && npi0==0 && npim==1) evtype = 9;  
           }
           if(proc.IsWeakCC() && init.IsNuBarN()) {
             // \bar{v} n -> l+ n pi- 
             if(np==1 && nn==0 && npip==1 && npi0==0 && npim==0) evtype = 10; 
           }
           if(proc.IsWeakCC() && init.IsNuBarP()) {
             // \bar{v} p -> l+ n pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 11; 
             // \bar{v} p -> l+ p pi- 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 12; 
           }
           if(proc.IsWeakNC() && init.IsNuBarP()) {
             // \bar{v} p -> \bar{v} p pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 13;
             // \bar{v} p -> \bar{v} n pi+ 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 14; 
           }
           if(proc.IsWeakNC() && init.IsNuBarN()) {
             // \bar{v} n -> \bar{v} n pi0 
             if(np==0 && nn==1 && npip==0 && npi0==1 && npim==0) evtype = 15;
             // \bar{v} n -> \bar{v} p pi-  
             if(np==1 && nn==0 && npip==0 && npi0==0 && npim==1) evtype = 16;
           }
        }
        */
    	output << "$ nuance " << evtype << endl;
    } // nuance code

    // add 'NEUT'-like event type
    // a list of NEUT event types can be seen here: http://t2k.phy.duke.edu/bin/view/Main/NeutModes
    // 
    else if(gOptOutFileFormat==11) {
    	int evtype = 0;

        const ProcessInfo &  proc = interaction->ProcInfo();
        const InitialState & init = interaction->InitState();
        const XclsTag &      xcls = interaction->ExclTag();
        const Kinematics &   kine = interaction->Kine();
        const Target &       tgt  = init.Tgt();

        bool is_cc    = proc.IsWeakCC();
        bool is_nc    = proc.IsWeakNC();
        bool is_charm = xcls.IsCharmEvent();
        bool is_qel   = proc.IsQuasiElastic();
        bool is_dis   = proc.IsDeepInelastic();
        bool is_res   = proc.IsResonant();
        bool is_cohpi = proc.IsCoherentPiProd();
      //bool is_ve    = proc.IsNuElectronElastic();
      //bool is_imd   = proc.IsInverseMuDecay();
        bool is_p     = tgt.HitNucIsSet() ? tgt.HitNucPdg()==kPdgProton  : false;
        bool is_n     = tgt.HitNucIsSet() ? tgt.HitNucPdg()==kPdgNeutron : false;
        bool is_nu    = pdg::IsNeutrino    (init.ProbePdg());
        bool is_nubar = pdg::IsAntiNeutrino(init.ProbePdg());
        bool W_gt_2   = kine.KVSet(kKVW) ?  (kine.W() > 2.0) : false;

        // (quasi-)elastic, nc+cc, nu+nubar
        //
        if      (is_qel && !is_charm && is_cc && is_nu           ) evtype =   1;
        else if (is_qel && !is_charm && is_nc && is_nu && is_p   ) evtype =  51;
        else if (is_qel && !is_charm && is_nc && is_nu && is_n   ) evtype =  52;
        else if (is_qel && !is_charm && is_cc && is_nubar        ) evtype =  -1;
        else if (is_qel && !is_charm && is_nc && is_nubar && is_p) evtype = -51;
        else if (is_qel && !is_charm && is_nc && is_nubar && is_n) evtype = -52;

        // coherent pi, nc+cc, nu+nubar
        //
        else if (is_cohpi && is_cc && is_nu   ) evtype =  16;        
        else if (is_cohpi && is_cc && is_nubar) evtype = -16;        
        else if (is_cohpi && is_nc && is_nu   ) evtype =  36;
        else if (is_cohpi && is_nc && is_nubar) evtype = -36;

        // dis, W>2, nc+cc, nu+nubar 
        // (charm DIS not simulated by NEUT, will bundle GENIE charm DIS into this category)
        //
        else if (is_dis && W_gt_2 && is_cc && is_nu   ) evtype =  26;        
        else if (is_dis && W_gt_2 && is_nc && is_nu   ) evtype =  46;        
        else if (is_dis && W_gt_2 && is_cc && is_nubar) evtype = -26;        
        else if (is_dis && W_gt_2 && is_nc && is_nubar) evtype = -46;        

        // resonance or dis with W < 2 GeV
        //
        else if ( is_res || (is_dis && !W_gt_2) ) {

            LOG("gntpc", pNOTICE) << "Current event is RES or DIS with W<2";

           // check the number of pions and nucleons in the primary hadronic system
           // (_before_ intranuclear rescattering)
           //
           int nn=0, np=0, npi0=0, npip=0, npim=0, nKp=0, nKm=0, nK0=0, neta=0, nlambda=0, ngamma=0;
           bool nuclear_target = init.Tgt().IsNucleus();
           while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) 
           {
               GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();  
               int ghep_pdgc = p->Pdg();

               // For nuclear targets use hadrons marked as 'hadron in the nucleus'
               // which are the ones passed in the intranuclear rescattering
               // For free nucleon targets use particles marked as 'final state'
               // but make an exception for pi0 which must have been decayed by then

               bool count_it = 
                     ( nuclear_target && ghep_ist==kIStHadronInTheNucleus) ||
                     (!nuclear_target && ghep_ist==kIStStableFinalState  ) ||
                     (!nuclear_target && ghep_ist==kIStDecayedState && ghep_pdgc==kPdgPi0);
               if(!count_it) continue;

               if(ghep_pdgc == kPdgProton )    np++;            // p
               if(ghep_pdgc == kPdgNeutron)    nn++;            // n
               if(ghep_pdgc == kPdgPiP)        npip++;          // pi+
               if(ghep_pdgc == kPdgPiM)        npim++;          // pi-
               if(ghep_pdgc == kPdgPi0)        npi0++;          // pi0
               if(ghep_pdgc == kPdgEta)        neta++;          // eta0
               if(ghep_pdgc == kPdgKP)         nKp++;           // K+
               if(ghep_pdgc == kPdgKM)         nKm++;           // K-
               if(ghep_pdgc == kPdgK0)         nK0++;           // K0
               if(ghep_pdgc == kPdgAntiK0)     nK0++;           // K0
               if(ghep_pdgc == kPdgLambda)     nlambda++;       // Lamda
               if(ghep_pdgc == kPdgAntiLambda) nlambda++;       // Lamda
               if(ghep_pdgc == kPdgGamma)      ngamma++;        // photon
           }
           LOG("gntpc", pNOTICE) 
              << "Num of primary hadrons: p = " << np << ", n = " << nn
              << ", pi+ = " << npip << ", pi- = " << npim << ", pi0 = " << npi0;

           int nnuc = np + nn;
           int npi  = npi0 + npip + npim;
           int nK   = nK0 + nKp + nKm;
           int neKL = neta + nK + nlambda;
 
           bool is_single_pi_dis = (npi==1) && is_dis;
           bool is_radiative_dec = (nnuc==1) && (npi==0) && (ngamma==1);

           // res + non-res bkg (single pi dis, W < 2 GeV)
           //
           if(is_res || is_single_pi_dis) {

              //
              // single gamma from resonances
              //

              if      (is_res && is_nu    && is_cc && is_n && is_radiative_dec) evtype =  17;
              else if (is_res && is_nu    && is_nc && is_n && is_radiative_dec) evtype =  38;
              else if (is_res && is_nu    && is_nc && is_p && is_radiative_dec) evtype =  39;

              else if (is_res && is_nubar && is_cc && is_p && is_radiative_dec) evtype = -17;
              else if (is_res && is_nubar && is_nc && is_n && is_radiative_dec) evtype = -38;
              else if (is_res && is_nubar && is_nc && is_p && is_radiative_dec) evtype = -39;

              //
              // single pi (res + non-res bkg)
              //

              // nu CC
              else if (is_nu    && is_cc && is_p && np==1 && nn==0 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  11;
              else if (is_nu    && is_cc && is_n && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  12;
              else if (is_nu    && is_cc && is_n && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  13;

              // nubar CC
              else if (is_nu    && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  31;
              else if (is_nu    && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  32;
              else if (is_nu    && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype =  33;
              else if (is_nu    && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  34;

              //nubar CC
              else if (is_nubar && is_cc && is_n && np==0 && nn==1 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -11;
              else if (is_nubar && is_cc && is_p && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -12;
              else if (is_nubar && is_cc && is_p && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -13;

              //nubar NC
              else if (is_nubar && is_nc && is_n && np==0 && nn==1 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -31; 
              else if (is_nubar && is_nc && is_p && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype = -32;
              else if (is_nubar && is_nc && is_n && np==1 && nn==0 && npip==0 && npim==1 && npi0==0 && neKL==0) evtype = -33;
              else if (is_nubar && is_nc && is_p && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype = -34;

              //
              // single eta from res
              //

              else if (is_res &&  is_nu    && is_cc && is_n && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  22;
              else if (is_res &&  is_nu    && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  42;
              else if (is_res &&  is_nu    && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype =  43;

              else if (is_res &&  is_nubar && is_cc && is_p && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -22;
              else if (is_res &&  is_nubar && is_nc && is_n && np==0 && nn==1 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -42;
              else if (is_res &&  is_nubar && is_nc && is_p && np==1 && nn==0 && npi==0 && nK==0 && nlambda==0 && neta==1) evtype = -43;

              //
              // single K from res
              //

              else if (is_res &&  is_nu    && is_cc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  23;
              else if (is_res &&  is_nu    && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  44;
              else if (is_res &&  is_nu    && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype =  45;

              else if (is_res &&  is_nubar && is_cc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -23;
              else if (is_res &&  is_nubar && is_nc && is_n && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -44;
              else if (is_res &&  is_nubar && is_nc && is_p && nnuc==0 && npi==0 && nK==1 && nlambda==1 && neta==0) evtype = -45;
           }

           // multi-pi (1.3 GeV < W < 2.0 GeV)
           //
           else {
              if      (is_nu    && is_cc) evtype =  21;
              else if (is_nu    && is_nc) evtype =  41;
              else if (is_nubar && is_cc) evtype = -21;
              else if (is_nubar && is_nc) evtype = -41;
           }
        }
        LOG("gntpc", pNOTICE) << "NEUT-like event type = " << evtype;
    	output << "$ nuance " << evtype << endl;

    } //neut code

    // add genie "event type"
    //
    else {
    	output << "$ genie " << interaction->AsString() << endl;
    }

    // add tracker vertex info
    double vtxx = 0, vtxy = 0, vtxz = 0, vtxt = 0;
    output << "$ vertex " << vtxx << " " << vtxy
           << " " << vtxz << " " << vtxt << " " << endl;

    // add 'tracks' (GENIE's equivalent of GHepParticles)
    bool info_added  = false;
    event_iter.Reset();
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) 
    {
       // Neglect all GENIE pseudo-particles
       //
       if(p->IsFake()) continue;

       // Convert GENIE's GHEP pdgc & status to NUANCE's equivalent
       //
       GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();
       int ist;
       switch (ghep_ist) {
         case kIStInitialState:             ist = -1;   break;
         case kIStStableFinalState:         ist =  0;   break;
         case kIStIntermediateState:        ist = -2;   break;
         case kIStDecayedState:             ist = -2;   break;
         case kIStNucleonTarget:            ist = -1;   break;
         case kIStDISPreFragmHadronicState: ist = -2;   break;
         case kIStPreDecayResonantState:    ist = -2;   break;
         case kIStHadronInTheNucleus:       ist = -2;   break;
         case kIStUndefined:                ist = -999; break;
         default:                           ist = -999; break;
       }
       // Convert GENIE pdg code -> nuance PDG code
       // For most particles both generators use the standard PDG codes.
       // For nuclei GENIE follows the PDG-convention: 10LZZZAAAI
       // NUANCE is using: ZZZAAA
       int ghep_pdgc = p->Pdg();
       int pdgc = ghep_pdgc;
       if ( p->IsNucleus() ) {
         int Z = pdg::IonPdgCodeToZ(ghep_pdgc);
         int A = pdg::IonPdgCodeToA(ghep_pdgc);
         pdgc = 1000*Z + A;
       }
       //
       //
       if(gOptOutFileFormat==11) {
         if(pdgc==kPdgK0 || pdgc==kPdgAntiK0) {
            RandomGen * rnd = RandomGen::Instance();
            double R =  rnd->RndGen().Rndm();
            if(R>0.5) pdgc = kPdgK0L;
            else      pdgc = kPdgK0S;
         }
       }
       // Get particle's energy & momentum
       TLorentzVector * p4 = p->P4();
       double E  = p4->Energy() / units::MeV;
       double Px = p4->Px()     / units::MeV;
       double Py = p4->Py()     / units::MeV;
       double Pz = p4->Pz()     / units::MeV;
       double P  = p4->P()      / units::MeV;
       // Compute direction cosines
       double dcosx = (P>0) ? Px/P : -999;
       double dcosy = (P>0) ? Py/P : -999;
       double dcosz = (P>0) ? Pz/P : -999;

       GHepStatus_t gist = (GHepStatus_t) p->Status();
       bool is_init =
             (gist == kIStInitialState || gist == kIStNucleonTarget);

       if(!is_init && !info_added) {
         // Add nuance obsolete and flux info (not filled in by
         // GENIE here). Add it once after the initial state particles
         output << "$ info 2 949000 0.0000E+00" << endl;
         info_added = true;
       }
       output << "$ track " << pdgc << " " << E << " "
              << dcosx << " " << dcosy << " " << dcosz << " "
              << ist << endl;
    }
    //add  tracker end tag
    output << "$ end" << endl;

    mcrec->Clear();
  } // event loop

  // add tracker end-of-file tag
  output << "$ stop" << endl;

  output.close();
  fin.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//___________________________________________________________________
//  *** GENIE GHEP EVENT TREE -> A T2K ROOT TRACKER TYPE FORMAT ***
//___________________________________________________________________
void ConvertToGT2KRooTracker(void)
{
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__

  //-- get pdglib
  PDGLibrary * pdglib = PDGLibrary::Instance();

  //-- define the output tree branches
  TBits*      brEvtFlags = 0;             // generator-specific event flags
  TObjString* brEvtCode = 0;              // generator-specific string with 'event code'
  int         brEvtNum;                   // event num.
  double      brEvtXSec;                  // cross section for selected event (1E-38 cm2)
  double      brEvtDXSec;                 // cross section for selected event kinematics (1E-38 cm2 /{K^n})
  double      brEvtWght;                  // weight for that event
  double      brEvtProb;                  // probability for that event (given cross section, path lengths, etc)
  double      brEvtVtx[4];                // event vertex position in detector coord syst (SI)
  int         brStdHepN;                  // number of particles in particle array 
  // > stdhep-like particle array:
  int         brStdHepPdg   [kNPmax];     // pdg codes (& generator specific codes for pseudoparticles)
  int         brStdHepStatus[kNPmax];     // generator-specific status code
  double      brStdHepX4    [kNPmax][4];  // 4-x (x, y, z, t) of particle in hit nucleus frame (fm)
  double      brStdHepP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  double      brStdHepPolz  [kNPmax][3];  // polarization vector
  int         brStdHepFd    [kNPmax];     // first daughter
  int         brStdHepLd    [kNPmax];     // last  daughter 
  int         brStdHepFm    [kNPmax];     // first mother
  int         brStdHepLm    [kNPmax];     // last  mother
  // > neutrino parent info (passed-through from the beam-line MC / quantities in 'jnubeam' units)
  int         brNuParentPdg;              // parent hadron pdg code
  int         brNuParentDecMode;          // parent hadron decay mode
  double      brNuParentDecP4 [4];        // parent hadron 4-momentum at decay 
  double      brNuParentDecX4 [4];        // parent hadron 4-position at decay
  double      brNuParentProP4 [4];        // parent hadron 4-momentum at production
  double      brNuParentProX4 [4];        // parent hadron 4-position at production
  int         brNuParentProNVtx;          // parent hadron vtx id

  //-- open the output ROOT file
  TFile fout(gOptOutFileName.c_str(), "RECREATE");

  //-- create the output ROOT tree
  TTree * rootracker_tree = new TTree("gRooTracker","GENIE event tree for T2K in rootracker format");

  //-- create the output ROOT tree branches
  rootracker_tree->Branch("EvtFlags", "TBits",      &brEvtFlags, 32000, 1);           
  rootracker_tree->Branch("EvtCode",  "TObjString", &brEvtCode,  32000, 1);            
  rootracker_tree->Branch("EvtNum",          &brEvtNum,          "EvtNum/I");             
  rootracker_tree->Branch("EvtXSec",         &brEvtXSec,         "EvtXSec/D");            
  rootracker_tree->Branch("EvtDXSec",        &brEvtDXSec,        "EvtDXSec/D");           
  rootracker_tree->Branch("EvtWght",         &brEvtWght,         "EvtWght/D");            
  rootracker_tree->Branch("EvtProb",         &brEvtProb,         "EvtProb/D");            
  rootracker_tree->Branch("EvtVtx",           brEvtVtx,          "EvtVtx[4]/D");             
  rootracker_tree->Branch("StdHepN",         &brStdHepN,         "StdHepN/I");              
  rootracker_tree->Branch("StdHepPdg",        brStdHepPdg,       "StdHepPdg[StdHepN]/I");  
  rootracker_tree->Branch("StdHepStatus",     brStdHepStatus,    "StdHepStatus[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepX4",         brStdHepX4,        "StdHepX4[StdHepN][4]/D"); 
  rootracker_tree->Branch("StdHepP4",         brStdHepP4,        "StdHepP4[StdHepN][4]/D"); 
  rootracker_tree->Branch("StdHepPolz",       brStdHepPolz,      "StdHepPolz[StdHepN][3]/D"); 
  rootracker_tree->Branch("StdHepFd",         brStdHepFd,        "StdHepFd[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepLd",         brStdHepLd,        "StdHepLd[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepFm",         brStdHepFm,        "StdHepFm[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepLm",         brStdHepLm,        "StdHepLm[StdHepN]/I"); 
  rootracker_tree->Branch("NuParentPdg",     &brNuParentPdg,     "NuParentPdg/I");       
  rootracker_tree->Branch("NuParentDecMode", &brNuParentDecMode, "NuParentDecMode/I");   
  rootracker_tree->Branch("NuParentDecP4",    brNuParentDecP4,   "NuParentDecP4[4]/D");     
  rootracker_tree->Branch("NuParentDecX4",    brNuParentDecX4,   "NuParentDecX4[4]/D");     
  rootracker_tree->Branch("NuParentProP4",    brNuParentProP4,   "NuParentProP4[4]/D");     
  rootracker_tree->Branch("NuParentProX4",    brNuParentProX4,   "NuParentProX4[4]/D");     
  rootracker_tree->Branch("NuParentProNVtx", &brNuParentProNVtx, "NuParentProNVtx/I");   

  //-- open the input GENIE ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           gtree = 0;
  NtpMCTreeHeader * thdr  = 0;
  gtree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr  = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- get mc record
  NtpMCEventRecord * mcrec = 0;
  gtree->SetBranchAddress("gmcrec", &mcrec);

  flux::GJPARCNuFluxPassThroughInfo * flux_info = 0;
  gtree->SetBranchAddress("flux", &flux_info);

  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
      gtree->GetEntries() : TMath::Min(gtree->GetEntries(), gOptN);
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  //-- event loop
  for(Long64_t iev = 0; iev < nmax; iev++) {
    gtree->GetEntry(iev);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);
    Interaction * interaction = event.Summary();

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;
    LOG("gntpc", pINFO) << *interaction;
    if(flux_info) {
       LOG("gntpc", pINFO) << *flux_info;
    } else {
       LOG("gntpc", pINFO) << "No flux info associated with that event";
    }

    //
    // clear output tree branches
    //

    brEvtFlags  = 0;
    brEvtCode   = 0;
    brEvtNum    = 0;    
    brEvtXSec   = 0;
    brEvtDXSec  = 0;
    brEvtWght   = 0;
    brEvtProb   = 0;
    for(int k=0; k<4; k++) { 
      brEvtVtx[k] = 0;
    }
    brStdHepN = event.GetEntries(); 
    for(int i=0; i<kNPmax; i++) {
       brStdHepPdg   [i] = 0;  
       brStdHepStatus[i] = 0;  
       for(int k=0; k<4; k++) {
         brStdHepX4 [i][k] = 0;  
         brStdHepP4 [i][k] = 0;  
       }
       for(int k=0; k<3; k++) {
         brStdHepPolz [i][k] = 0;  
       }
       brStdHepFd    [i] = 0;  
       brStdHepLd    [i] = 0;  
       brStdHepFm    [i] = 0;  
       brStdHepLm    [i] = 0;  
    }
    brNuParentPdg     = 0;           
    brNuParentDecMode = 0;       
    for(int k=0; k<4; k++) {  
      brNuParentDecP4 [k] = 0;     
      brNuParentDecX4 [k] = 0;     
      brNuParentProP4 [k] = 0;     
      brNuParentProX4 [k] = 0;     
    }
    brNuParentProNVtx = 0;     

    //
    // copy current event info to output tree
    //

    brEvtFlags  = new TBits(*event.EventFlags());   
    brEvtCode   = new TObjString(event.Summary()->AsString().c_str());   
    brEvtNum    = (int) iev;    
    brEvtXSec   = (1E+38/units::cm2) * event.XSec();    
    brEvtDXSec  = (1E+38/units::cm2) * event.DiffXSec();    
    brEvtWght   = event.Weight();    
    brEvtProb   = event.Probability();    
    brEvtVtx[0] = event.Vertex()->X();    
    brEvtVtx[1] = event.Vertex()->Y();    
    brEvtVtx[2] = event.Vertex()->Z();    
    brEvtVtx[3] = event.Vertex()->T();    

    brStdHepN = event.GetEntries(); 

    int iparticle=0;
    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
        assert(p);
        brStdHepPdg   [iparticle] = p->Pdg(); 
        brStdHepStatus[iparticle] = (int) p->Status(); 
        brStdHepX4    [iparticle][0] = p->X4()->X(); 
        brStdHepX4    [iparticle][1] = p->X4()->Y(); 
        brStdHepX4    [iparticle][2] = p->X4()->Z(); 
        brStdHepX4    [iparticle][3] = p->X4()->T(); 
        brStdHepP4    [iparticle][0] = p->P4()->Px(); 
        brStdHepP4    [iparticle][1] = p->P4()->Py(); 
        brStdHepP4    [iparticle][2] = p->P4()->Pz(); 
        brStdHepP4    [iparticle][3] = p->P4()->E(); 
        if(p->PolzIsSet()) {
          brStdHepPolz  [iparticle][0] = TMath::Sin(p->PolzPolarAngle()) * TMath::Cos(p->PolzAzimuthAngle());
          brStdHepPolz  [iparticle][1] = TMath::Sin(p->PolzPolarAngle()) * TMath::Sin(p->PolzAzimuthAngle());
          brStdHepPolz  [iparticle][2] = TMath::Cos(p->PolzPolarAngle());
        }
        brStdHepFd    [iparticle] = p->FirstDaughter(); 
        brStdHepLd    [iparticle] = p->LastDaughter(); 
        brStdHepFm    [iparticle] = p->FirstMother(); 
        brStdHepLm    [iparticle] = p->LastMother(); 
        iparticle++;
    }
    // Copy flux info - that may not be available, eg if events were generated using 
    // plain flux histograms - not the beam simulation's output flux ntuples
    if(flux_info) {
       brNuParentPdg       = flux_info->pdg;        
       brNuParentDecMode   = flux_info->decayMode;        

       brNuParentDecP4 [0] = flux_info->decayP * flux_info->decayDirX; // px
       brNuParentDecP4 [1] = flux_info->decayP * flux_info->decayDirY; // py
       brNuParentDecP4 [2] = flux_info->decayP * flux_info->decayDirZ; // px
       brNuParentDecP4 [3] = TMath::Sqrt(
                                 TMath::Power(pdglib->Find(flux_info->pdg)->Mass(), 2.)
                               + TMath::Power(flux_info->decayP, 2.)
                              ); // E
       brNuParentDecX4 [0] = flux_info->decayX; // x
       brNuParentDecX4 [1] = flux_info->decayY; // y       
       brNuParentDecX4 [2] = flux_info->decayZ; // x   
       brNuParentDecX4 [3] = 0;                 // t

       brNuParentProP4 [0] = flux_info->prodP * flux_info->prodDirX; // px
       brNuParentProP4 [1] = flux_info->prodP * flux_info->prodDirY; // py
       brNuParentProP4 [2] = flux_info->prodP * flux_info->prodDirZ; // px
       brNuParentProP4 [3] = TMath::Sqrt(
                                TMath::Power(pdglib->Find(flux_info->pdg)->Mass(), 2.)
                              + TMath::Power(flux_info->prodP, 2.)
                              ); // E
       brNuParentProX4 [0] = flux_info->prodX; // x
       brNuParentProX4 [1] = flux_info->prodY; // y       
       brNuParentProX4 [2] = flux_info->prodZ; // x   
       brNuParentProX4 [3] = 0;                // t

       brNuParentProNVtx   = flux_info->prodNVtx;
    }
    rootracker_tree->Fill();
    mcrec->Clear();

  } // event loop

  // Copy POT normalization for the generated sample
  double pot = gtree->GetWeight();
  rootracker_tree->SetWeight(pot);

  fin.Close();

  fout.Write();
  fout.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";

#else

  LOG("gntpc", pWARN) 
    << "\n You should enable --with-flux-drivers during the GENIE"
    << " installation step so as to access the JPARC neutrino flux"
    << " pass-through info for the simulated neutrino interactions.";

#endif
}
//___________________________________________________________________
// * GENIE GHEP EVENT TREE -> NEUGEN-style format for AGKY studies *
//___________________________________________________________________
void ConvertToGHad(void)
{
// Neugen-style text format for the AGKY hadronization model studies
// Format:
// (blank line) 
// event number, neutrino particle code, CCNC, IM, A, Z
// int_type, x, y, w, ihadmod 
// neutrino particle code, 5 vec
// lepton particle code, 5-vec
// outgoing hadronic system, 5-vec
// number of stable daughters of hadronic system
// ... then for each stable daughter
// particle id, 5 vec 

  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- get mc record
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- open the output stream
  ofstream output(gOptOutFileName.c_str(), ios::out);

  //-- open output root file and create ntuple -- if required
#ifdef __GHAD_NTP__
  TFile fout("ghad.root","recreate");  
  TTree * ghad = new TTree("ghad","");   
  ghad->Branch("i",       &brIev,          "i/I " );
  ghad->Branch("W",       &brW,            "W/D " );
  ghad->Branch("n",       &brN,            "n/I " );
  ghad->Branch("pdg",      brPdg,          "pdg[n]/I " );
  ghad->Branch("E",        brE,            "E[n]/D"    );
  ghad->Branch("px",       brPx,           "px[n]/D"   );
  ghad->Branch("py",       brPy,           "py[n]/D"   );
  ghad->Branch("pz",       brPz,           "pz[n]/D"   );
#endif

  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
      tree->GetEntries() : TMath::Min(tree->GetEntries(), gOptN);
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  //-- event loop
  for(Long64_t iev = 0; iev < nmax; iev++) {
    tree->GetEntry(iev);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

#ifdef __GHAD_NTP__
    brN = 0;  
    for(int k=0; k<kNPmax; k++) {
      brPdg[k]=0;       
      brE  [k]=0;  
      brPx [k]=0; 
      brPy [k]=0; 
      brPz [k]=0;  
    }
#endif

    //
    // convert the current event
    //
    const Interaction * interaction = event.Summary();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    const InitialState & init_state = interaction->InitState();

    bool is_dis = proc_info.IsDeepInelastic();
    bool is_res = proc_info.IsResonant();
    bool is_cc  = proc_info.IsWeakCC();

    bool pass   = is_cc && (is_dis || is_res);
    if(!pass) {
      mcrec->Clear();
      continue;
    }

    int ccnc   = is_cc ? 1 : 0;
    int inttyp = 3; 

    int im     = -1;
    if      (init_state.IsNuP    ()) im = 1; 
    else if (init_state.IsNuN    ()) im = 2; 
    else if (init_state.IsNuBarP ()) im = 3; 
    else if (init_state.IsNuBarN ()) im = 4; 
    else return;

    GHepParticle * neutrino = event.Probe();
    assert(neutrino);
    GHepParticle * target = event.Particle(1);
    assert(target);
    GHepParticle * fsl = event.FinalStatePrimaryLepton();
    assert(fsl);
    GHepParticle * hitnucl = event.HitNucleon();
    assert(hitnucl);

    int nupdg  = neutrino->Pdg();
    int fslpdg = fsl->Pdg();
    int A      = target->A();
    int Z      = target->Z();

    const TLorentzVector & k1 = *(neutrino->P4());  // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());       // l 4-p (k2)
//  const TLorentzVector & p1 = *(hitnucl->P4());   // N 4-p (p1)      
//  const TLorentzVector & ph = *(hadsyst->P4());   // had-syst 4-p 
  
    TLorentzVector ph;
    if(is_dis) {
      GHepParticle * hadsyst = event.FinalStateHadronicSystem();
      assert(hadsyst);
      ph = *(hadsyst->P4());
    }   
    if(is_res) {
      GHepParticle * hadres = event.Particle(hitnucl->FirstDaughter());
      ph = *(hadres->P4());
    }
   
    const Kinematics & kine = interaction->Kine();
    bool get_selected = true;
    double x  = kine.x (get_selected);
    double y  = kine.y (get_selected);
    double W  = kine.W (get_selected);

    int hadmod  = -1;
    int ihadmom = -1;
    TIter event_iter(&event);
    GHepParticle * p = 0;
    int i=-1;
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
      i++;
      int pdg = p->Pdg();
      if (pdg == kPdgHadronicSyst )  { hadmod= 2; ihadmom=i; }
      if (pdg == kPdgString       )  { hadmod=11; ihadmom=i; }
      if (pdg == kPdgCluster      )  { hadmod=12; ihadmom=i; }
      if (pdg == kPdgIndep        )  { hadmod=13; ihadmom=i; }
    }

    output << endl;
    output << iev    << "\t"  
           << nupdg  << "\t"  << ccnc << "\t"  << im << "\t"  
           << A      << "\t"  << Z << endl;
    output << inttyp << "\t" << x << "\t" << y << "\t" << W << "\t" 
           << hadmod << endl;
    output << nupdg       << "\t"
           << k1.Px()     << "\t" << k1.Py() << "\t" << k1.Pz() << "\t"
           << k1.Energy() << "\t" << k1.M()  << endl;
    output << fslpdg      << "\t"
           << k2.Px()     << "\t" << k2.Py() << "\t" << k2.Pz() << "\t"
           << k2.Energy() << "\t" << k2.M()  << endl;
    output << 111111 << "\t"
           << ph.Px()     << "\t" << ph.Py() << "\t" << ph.Pz() << "\t"
           << ph.Energy() << "\t" << ph.M()  << endl;

    vector<int> hadv;

    event_iter.Reset();
    i=-1;
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
      i++;
      if(i<ihadmom) continue;

      GHepStatus_t ist = p->Status();
      int pdg = p->Pdg();

      if(ist == kIStDISPreFragmHadronicState) continue;

      if(ist == kIStStableFinalState) {
        GHepParticle * mom = event.Particle(p->FirstMother());
        GHepStatus_t mom_ist = mom->Status();
        int mom_pdg = mom->Pdg();
        bool skip = (mom_pdg == kPdgPi0 && mom_ist== kIStDecayedState);
        if(!skip) { hadv.push_back(i); }
      }

      if(pdg==kPdgPi0 && ist==kIStDecayedState) { hadv.push_back(i); }
    }

    output << hadv.size() << endl;

#ifdef __GHAD_NTP__
    brIev = (int) iev;   
    brW   = W;  
    brN   = hadv.size();
    int k=0;
#endif

    vector<int>::const_iterator hiter = hadv.begin();
    for( ; hiter != hadv.end(); ++hiter) {
      int id = *hiter;
      GHepParticle * p = event.Particle(id);
      int pdg = p->Pdg();
      double px = p->P4()->Px();
      double py = p->P4()->Py();
      double pz = p->P4()->Pz();
      double E  = p->P4()->Energy();
      double m  = p->P4()->M();
      output << pdg << "\t" 
             << px  << "\t" << py << "\t" << pz << "\t"
             << E   << "\t" << m  << endl;

#ifdef __GHAD_NTP__
      brPx[k]  = px;
      brPy[k]  = py;
      brPz[k]  = pz;
      brE[k]   = E;
      brPdg[k] = pdg;
      k++;
#endif
    }

#ifdef __GHAD_NTP__
    ghad->Fill();
#endif

    mcrec->Clear();

  } // event loop

  output.close();
  fin.Close();

#ifdef __GHAD_NTP__
  ghad->Write("ghad");
  fout.Write();
  fout.Close();
#endif

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//___________________________________________________________________
// * GENIE GHEP EVENT TREE -> Summary tree for INTRANUKE studies *
//___________________________________________________________________
void ConvertToGINuke(void)
{
  //-- open output file & create output summary tree & create the tree branches
  //
  LOG("gntpc", pNOTICE)
       << "*** Saving summary tree to: " << gOptOutFileName;
  TFile fout(gOptOutFileName.c_str(),"recreate");
  
  TTree * s_tree = new TTree("ginuke","GENIE INuke Summary Tree");
  assert(s_tree);

  //-- define summary ntuple
  //

  //
  // ... 
  // ... add code here
  // ... 
  //

  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           er_tree = 0;
  NtpMCTreeHeader * thdr    = 0;
  er_tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr    = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );
  if (!er_tree) {
    LOG("gntpc", pERROR) << "Null input tree";
    return;
  }
  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- get the mc record
  NtpMCEventRecord * mcrec = 0;
  er_tree->SetBranchAddress("gmcrec", &mcrec);
  if (!mcrec) {
    LOG("gntpc", pERROR) << "Null MC record";
    return;
  }
  
  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       er_tree->GetEntries() : TMath::Min( er_tree->GetEntries(), gOptN );
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  for(Long64_t iev = 0; iev < nmax; iev++) {

    er_tree->GetEntry(iev);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;


    // analyze current event and fill the summary ntuple

    //
    // ...
    // ... add code here
    //
    // ...


    mcrec->Clear();

  } // event loop

  fin.Close();
    
  fout.Write();
  fout.Close();
       
  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//___________________________________________________________________
//            FUNCTIONS FOR PARSING CMD-LINE ARGUMENTS 
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  //get input ROOT file (containing a GENIE ER ntuple)
  try {
    LOG("gntpc", pINFO) << "Reading input filename";
    gOptInpFileName = utils::clap::CmdLineArgAsString(argc,argv,'i');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gntpc", pFATAL)
               << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //check input GENIE ROOT file
  bool inpok = !(gSystem->AccessPathName(gOptInpFileName.c_str()));
  if (!inpok) {
    LOG("gntpc", pFATAL)
           << "The input ROOT file ["
                       << gOptInpFileName << "] is not accessible";
    exit(2);
  }

  //get output file format
  try {
    LOG("gntpc", pINFO) << "Reading output file format";
    gOptOutFileFormat = utils::clap::CmdLineArgAsInt(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gntpc", pFATAL) << "Unspecified output file format";
      exit(3);
    }
  }

  //get output file name 
  try {
    LOG("gntpc", pINFO) << "Reading output filename";
    gOptOutFileName = utils::clap::CmdLineArgAsString(argc,argv,'o');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gntpc", pINFO)
                << "Unspecified output filename - Using default";
      gOptOutFileName = DefaultOutputFile();
    }
  }

  //get number of events to convert
  try {
    LOG("gntpc", pINFO) << "Reading number of events to analyze";
    gOptN = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gntpc", pINFO)
          << "Unspecified number of events to analyze - Use all";
      gOptN = -1;
    }
  }
}
//___________________________________________________________________
string DefaultOutputFile(void)
{
  // filename extension - depending on file format
  string ext="";
  if      (gOptOutFileFormat ==   0)  ext = "gst.root";
  else if (gOptOutFileFormat ==   1)  ext = "gtrac0.dat";
  else if (gOptOutFileFormat ==  11)  ext = "gtrac0.dat";
  else if (gOptOutFileFormat ==   2)  ext = "gtrac.dat";
  else if (gOptOutFileFormat ==   3)  ext = "gtrac.root";
  else if (gOptOutFileFormat == 100)  ext = "gxml";
  else if (gOptOutFileFormat == 901)  ext = "ghad.dat";
  else if (gOptOutFileFormat == 902)  ext = "ginuke.root";

  string inpname = gOptInpFileName;
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
  name << inpname << ext;

  return gSystem->BaseName(name.str().c_str());
}
//___________________________________________________________________
void PrintSyntax(void)
{
  string basedir  = string( gSystem->Getenv("GENIE") );
  string thisfile = basedir + string("/src/stdapp/gNtpConv.cxx");
  string cmd      = "less " + thisfile;

  gSystem->Exec(cmd.c_str());
}
//___________________________________________________________________
