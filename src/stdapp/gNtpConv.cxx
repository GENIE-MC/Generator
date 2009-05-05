//____________________________________________________________________________
/*!

\program gntpc

\brief   Converts a native GENIE (GHEP/ROOT) event tree file to a host of
         plain text, XML or bare-ROOT formats.

         Syntax:
           gntpc -i input_filename [-o output_filename] -f format [-n nev]

         Options :
           [] denotes an optional argument
           -n number of events to convert
           -f is a string that specifies the output file format. 

	      [T2K/GENIE formats]
   	       * `t2k_rootracker':
                     The standardized bare-ROOT GENIE event tree used by the 
                     nd280, INGRID and 2km MC.
                     Includes full information about the generated neutrino
                     event and pass-through JPARC flux info.
   	       * `t2k_tracker': 
                     A tracker-type format with tweaks required by the SuperK
                     detector MC (SKDETSIM):
                        - Converting K0, \bar{K0} -> KO_{long}, K0_{short}
                        - Emulating 'NEUT' reaction codes
                        - Appropriate $track ordering for SKDETSIM
                        - Passing detailed GENIE MC truth and JPARC flux info
                          using the tracker $info lines. This information, 
                          propaged by SKDETSIM to the DSTs, is identical with the 
                          one used at the near detectors and can be used for 
                          global systematic studies.

	      [Generic formats]
               * `gst': 
                    The 'definite' GENIE summary tree format (gst).
   	       * `gxml': 
                     GENIE XML event format 
   	       * `rootracker': 
                     Similar to `t2k_rootracker' but without the T2K-specific
                     pass-through JPARC flux-info

	      [GENIE test / cross-generator comparison formats]
   	       * `ghad': 
	             NEUGEN-style text-based format for hadronization studies
   	       * `ginuke': 
	             INTRANUKE summary ntuple for intranuclear-rescattering studies

	      [Other (depreciated) formats]
   	       * `nuance_tracker': 
   		     NUANCE-style tracker text-based format 

           -o specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
               `t2k_tracker'    -> *.gtrac.dat
               `t2k_rootracker' -> *.gtrac.root
               `gst'            -> *.gst.root
               `rootracker'     -> *.gtrac.root
               `gxml'           -> *.gxml 
               `ghad'           -> *.ghad.dat
               `ginuke'         -> *.ginuke.root
               `nuance_tracker' -> *.gtrac_legacy.dat
		
         Examples:
           (1)  shell% gntpc -i myfile.ghep.root -f t2k_rootracker

                Converts all events in the GHEP file myfile.ghep.root into the
                t2k_rootracker format. 
                The output file is named myfile.gtrac.root

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created September 23, 2005

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
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
#include "GHEP/GHepUtils.h"
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
void   ConvertToGTracker       (void);
void   ConvertToGRooTracker    (void);
void   ConvertToGXML           (void);
void   ConvertToGHad           (void);
void   ConvertToGINuke         (void);
void   GetCommandLineArgs      (int argc, char ** argv);
void   PrintSyntax             (void);
string DefaultOutputFile       (void);
bool   CheckRootFilename       (string filename);

//format enum
typedef enum EGNtpcFmt {
  kConvFmt_undef = 0,
  kConvFmt_t2k_rootracker,
  kConvFmt_t2k_tracker,
  kConvFmt_gst,
  kConvFmt_rootracker,
  kConvFmt_gxml,
  kConvFmt_ghad,
  kConvFmt_ginuke,
  kConvFmt_nuance_tracker
} GNtpcFmt_t;

//input options (from command line arguments):
string     gOptInpFileName;
string     gOptOutFileName;
GNtpcFmt_t gOptOutFileFormat;
Long64_t   gOptN; 

//consts
const int kNPmax = 100;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);

  //-- call the appropriate conversion function
  switch(gOptOutFileFormat) {
   case (kConvFmt_t2k_rootracker) :  
	ConvertToGRooTracker(); 
	break;
   case (kConvFmt_t2k_tracker)  :  
	ConvertToGTracker();        
	break;
   case (kConvFmt_gst)  :
	ConvertToGST();        
	break;  
   case (kConvFmt_rootracker) :  
	ConvertToGRooTracker(); 
	break;
   case (kConvFmt_gxml) :  
	ConvertToGXML();         
	break;
   case (kConvFmt_ghad) :  
	ConvertToGHad();         
	break;
   case (kConvFmt_ginuke) :  
	ConvertToGINuke();         
	break;
   case (kConvFmt_nuance_tracker)  :  
	ConvertToGTracker();        
	break;
   default:
     LOG("gntpc", pFATAL)
          << "Invalid output format [" << gOptOutFileFormat << "]";
     PrintSyntax();
     gAbortingInErr = true;
     exit(3);
  }
  return 0;
}
//___________________________________________________________________
// ***** GENIE GHEP EVENT TREE FORMAT -> GENIE SUMMARY NTUPLE *****
//___________________________________________________________________
void ConvertToGST(void)
{
  //-- some constants
  const double e_h = 1.3; // typical e/h ratio used for computing mean `calorimetric response'

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
  bool   brIsCoh       = false;  // Is COH?
  bool   brIsImd       = false;  // Is IMD?
  bool   brIsNuEL      = false;  // Is ve elastic?
  bool   brIsCC        = false;  // Is CC?
  bool   brIsNC        = false;  // Is NC?
  bool   brIsCharmPro  = false;  // Produces charm?
  int    brCodeNeut    = 0;      // The equivalent NEUT reaction code (if any)
  int    brCodeNuance  = 0;      // The equivalent NUANCE reaction code (if any)
  double brWeight      = 0;      // Event weight
  double brKineXs      = 0;      // Bjorken x as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineYs      = 0;      // Inelasticity y as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineTs      = 0;      // Energy transfer to nucleus at COH events as was generated during kinematical selection
  double brKineQ2s     = 0;      // Momentum transfer Q^2 as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineWs      = 0;      // Hadronic invariant mass W as was generated during kinematical selection; takes fermi momentum / off-shellness into account
  double brKineX       = 0;      // Experimental-like Bjorken x; neglects fermi momentum / off-shellness 
  double brKineY       = 0;      // Experimental-like inelasticity y; neglects fermi momentum / off-shellness 
  double brKineT       = 0;      // Experimental-like energy transfer to nucleus at COH events 
  double brKineQ2      = 0;      // Experimental-like momentum transfer Q^2; neglects fermi momentum / off-shellness
  double brKineW       = 0;      // Experimental-like hadronic invariant mass W; neglects fermi momentum / off-shellness 
  double brEv          = 0;      // Neutrino energy @ LAB
  double brPxv         = 0;      // Neutrino px @ LAB
  double brPyv         = 0;      // Neutrino py @ LAB
  double brPzv         = 0;      // Neutrino pz @ LAB
  double brEn          = 0;      // Initial state hit nucleon energy @ LAB
  double brPxn         = 0;      // Initial state hit nucleon px @ LAB
  double brPyn         = 0;      // Initial state hit nucleon py @ LAB
  double brPzn         = 0;      // Initial state hit nucleon pz @ LAB
  double brEl          = 0;      // Final state primary lepton energy @ LAB
  double brPxl         = 0;      // Final state primary lepton px @ LAB
  double brPyl         = 0;      // Final state primary lepton py @ LAB
  double brPzl         = 0;      // Final state primary lepton pz  @ LAB
  int    brNfP         = 0;      // Nu. of final state p's + \bar{p}'s (after intranuclear rescattering)
  int    brNfN         = 0;      // Nu. of final state n's + \bar{n}'s
  int    brNfPip       = 0;      // Nu. of final state pi+'s
  int    brNfPim       = 0;      // Nu. of final state pi-'s
  int    brNfPi0       = 0;      // Nu. of final state pi0's (
  int    brNfKp        = 0;      // Nu. of final state K+'s
  int    brNfKm        = 0;      // Nu. of final state K-'s
  int    brNfK0        = 0;      // Nu. of final state K0's + \bar{K0}'s
  int    brNfEM        = 0;      // Nu. of final state gammas and e-/e+ 
  int    brNfOther     = 0;      // Nu. of heavier final state hadrons (D+/-,D0,Ds+/-,Lamda,Sigma,Lamda_c,Sigma_c,...)
  int    brNiP         = 0;      // Nu. of `primary' (: before intranuclear rescattering) p's + \bar{p}'s  
  int    brNiN         = 0;      // Nu. of `primary' n's + \bar{n}'s  
  int    brNiPip       = 0;      // Nu. of `primary' pi+'s 
  int    brNiPim       = 0;      // Nu. of `primary' pi-'s 
  int    brNiPi0       = 0;      // Nu. of `primary' pi0's 
  int    brNiKp        = 0;      // Nu. of `primary' K+'s  
  int    brNiKm        = 0;      // Nu. of `primary' K-'s  
  int    brNiK0        = 0;      // Nu. of `primary' K0's + \bar{K0}'s 
  int    brNiEM        = 0;      // Nu. of `primary' gammas and e-/e+ 
  int    brNiOther     = 0;      // Nu. of other `primary' hadron shower particles
  int    brNf          = 0;      // Nu. of final state particles in hadronic system
  int    brPdgf[kNPmax];         // Pdg code of k^th final state particle in hadronic system
  double brEf  [kNPmax];         // Energy   of k^th final state particle in hadronic system @ LAB
  double brPxf [kNPmax];         // Px       of k^th final state particle in hadronic system @ LAB
  double brPyf [kNPmax];         // Py       of k^th final state particle in hadronic system @ LAB
  double brPzf [kNPmax];         // Pz       of k^th final state particle in hadronic system @ LAB
  int    brNi          = 0;      // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
  int    brPdgi[kNPmax];         // Pdg code of k^th particle in 'primary' hadronic system 
  double brEi  [kNPmax];         // Energy   of k^th particle in 'primary' hadronic system @ LAB
  double brPxi [kNPmax];         // Px       of k^th particle in 'primary' hadronic system @ LAB
  double brPyi [kNPmax];         // Py       of k^th particle in 'primary' hadronic system @ LAB
  double brPzi [kNPmax];         // Pz       of k^th particle in 'primary' hadronic system @ LAB
  double brVtxX;                 // Vertex x in detector coord system (SI)
  double brVtxY;                 // Vertex y in detector coord system (SI)
  double brVtxZ;                 // Vertex z in detector coord system (SI)
  double brVtxT;                 // Vertex t in detector coord system (SI)
  double brCalResp0;             // Approximate calorimetric response to the hadronic system computed as sum of
				 //  - (kinetic energy) for pi+, pi-, p, n 
                                 //  - (energy + mass)  for antiproton, antineutron
                                 //  - ((e/h) * energy)   for pi0, gamma, e-, e+, where e/h is set to 1.3
                                 //  - (kinetic energy) for other particles

  //-- open output file & create output summary tree & create the tree branches
  //
  LOG("gntpc", pNOTICE) 
       << "*** Saving summary tree to: " << gOptOutFileName;
  TFile fout(gOptOutFileName.c_str(),"recreate");

  TTree * s_tree = new TTree("gst","GENIE Summary Event Tree");

  //-- create tree branches
  //
  s_tree->Branch("iev",           &brIev,           "iev/I"         );
  s_tree->Branch("neu",	          &brNeutrino,      "neu/I"	    );
  s_tree->Branch("tgt",           &brTarget,        "tgt/I"	    );
  s_tree->Branch("Z",             &brTargetZ,       "Z/I"	    );
  s_tree->Branch("A",             &brTargetA,       "A/I"	    );
  s_tree->Branch("hitnuc",        &brHitNuc,        "hitnuc/I"      );
  s_tree->Branch("hitqrk",        &brHitQrk,        "hitqrk/I"      );
  s_tree->Branch("resid",         &brResId,	    "resid/I"	    );
  s_tree->Branch("sea",	          &brFromSea,       "sea/O"	    );
  s_tree->Branch("qel",	          &brIsQel,	    "qel/O"	    );
  s_tree->Branch("res",	          &brIsRes,	    "res/O"	    );
  s_tree->Branch("dis",	          &brIsDis,	    "dis/O"	    );
  s_tree->Branch("coh",           &brIsCoh,         "coh/O"	    );
  s_tree->Branch("imd",	          &brIsImd,	    "imd/O"	    );
  s_tree->Branch("nuel",          &brIsNuEL,        "nuel/O"	    );
  s_tree->Branch("cc",	          &brIsCC,	    "cc/O"	    );
  s_tree->Branch("nc",	          &brIsNC,	    "nc/O"	    );
  s_tree->Branch("charm",         &brIsCharmPro,    "charm/O"	    );
  s_tree->Branch("neut_code",     &brCodeNeut,      "neut_code/I"   );
  s_tree->Branch("nuance_code",   &brCodeNuance,    "nuance_code/I" );
  s_tree->Branch("wght",          &brWeight,        "wght/D"	    );
  s_tree->Branch("xs",	          &brKineXs,        "xs/D"	    );
  s_tree->Branch("ys",	          &brKineYs,        "ys/D"	    );
  s_tree->Branch("ts",	          &brKineTs,        "ts/D"	    );
  s_tree->Branch("Q2s",	          &brKineQ2s,       "Q2s/D"	    );
  s_tree->Branch("Ws",	          &brKineWs,        "Ws/D"	    );
  s_tree->Branch("x",	          &brKineX,	    "x/D"	    );
  s_tree->Branch("y",	          &brKineY,	    "y/D"	    );
  s_tree->Branch("t",	          &brKineT,	    "t/D"	    );
  s_tree->Branch("Q2",	          &brKineQ2,        "Q2/D"	    );
  s_tree->Branch("W",	          &brKineW,	    "W/D"	    );
  s_tree->Branch("Ev",	          &brEv,	    "Ev/D"	    );
  s_tree->Branch("pxv",	          &brPxv,	    "pxv/D"	    );
  s_tree->Branch("pyv",	          &brPyv,	    "pyv/D"	    );
  s_tree->Branch("pzv",	          &brPzv,	    "pzv/D"	    );
  s_tree->Branch("En",	          &brEn,	    "En/D"	    );
  s_tree->Branch("pxn",	          &brPxn,	    "pxn/D"	    );
  s_tree->Branch("pyn",	          &brPyn,	    "pyn/D"	    );
  s_tree->Branch("pzn",	          &brPzn,	    "pzn/D"	    );
  s_tree->Branch("El",	          &brEl,	    "El/D"	    );
  s_tree->Branch("pxl",	          &brPxl,	    "pxl/D"	    );
  s_tree->Branch("pyl",	          &brPyl,	    "pyl/D"	    );
  s_tree->Branch("pzl",	          &brPzl,	    "pzl/D"	    );
  s_tree->Branch("nfp",	          &brNfP,	    "nfp/I"	    );
  s_tree->Branch("nfn",	          &brNfN,	    "nfn/I"	    );
  s_tree->Branch("nfpip",         &brNfPip,	    "nfpip/I"	    );
  s_tree->Branch("nfpim",         &brNfPim,	    "nfpim/I"	    );
  s_tree->Branch("nfpi0",         &brNfPi0,	    "nfpi0/I"	    );
  s_tree->Branch("nfkp",          &brNfKp,	    "nfkp/I"	    );
  s_tree->Branch("nfkm",          &brNfKm,	    "nfkm/I"	    );
  s_tree->Branch("nfk0",          &brNfK0,	    "nfk0/I"	    );
  s_tree->Branch("nfem",          &brNfEM,	    "nfem/I"	    );
  s_tree->Branch("nfother",       &brNfOther,       "nfother/I"     );
  s_tree->Branch("nip",	          &brNiP,	    "np/I"	    );
  s_tree->Branch("nin",	          &brNiN,	    "nn/I"	    );
  s_tree->Branch("nipip",         &brNiPip,	    "npip/I"	    );
  s_tree->Branch("nipim",         &brNiPim,	    "npim/I"	    );
  s_tree->Branch("nipi0",         &brNiPi0,	    "npi0/I"	    );
  s_tree->Branch("nikp",          &brNiKp,	    "nkp/I"	    );
  s_tree->Branch("nikm",          &brNiKm,	    "nkm/I"	    );
  s_tree->Branch("nik0",          &brNiK0,	    "nk0/I"	    );
  s_tree->Branch("niem",          &brNiEM,	    "niem/I"	    );
  s_tree->Branch("niother",       &brNiOther,       "niother/I"     );
  s_tree->Branch("ni",	         &brNi,	            "ni/I"	    );
  s_tree->Branch("pdgi",          brPdgi,	    "pdgi[ni]/I "   );
  s_tree->Branch("Ei",	          brEi,	            "Ei[ni]/D"      );
  s_tree->Branch("pxi",	          brPxi,	    "pxi[ni]/D"     );
  s_tree->Branch("pyi",	          brPyi,	    "pyi[ni]/D"     );
  s_tree->Branch("pzi",	          brPzi,	    "pzi[ni]/D"     );
  s_tree->Branch("nf",	         &brNf,	            "nf/I"	    );
  s_tree->Branch("pdgf",          brPdgf,	    "pdgf[nf]/I "   );
  s_tree->Branch("Ef",	          brEf,	            "Ef[nf]/D"      );
  s_tree->Branch("pxf",	          brPxf,	    "pxf[nf]/D"     );
  s_tree->Branch("pyf",	          brPyf,	    "pyf[nf]/D"     );
  s_tree->Branch("pzf",	          brPzf,	    "pzf[nf]/D"     );
  s_tree->Branch("vtxx",         &brVtxX,	    "vtxx/D"        );
  s_tree->Branch("vtxy",         &brVtxY,	    "vtxy/D"        );
  s_tree->Branch("vtxz",         &brVtxZ,	    "vtxz/D"        );
  s_tree->Branch("vtxt",         &brVtxT,	    "vtxt/D"        );
  s_tree->Branch("calresp0",     &brCalResp0,	    "calresp0/D"    );

  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           er_tree = 0;
  NtpMCTreeHeader * thdr    = 0;
  er_tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr    = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );
  if (!er_tree) {
    LOG("gntpc", pERROR) << "Null input GHEP event tree";
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
    bool is_coh    = proc_info.IsCoherent();
    bool is_imd    = proc_info.IsInverseMuDecay();
    bool is_nuel   = proc_info.IsNuElectronElastic();
    bool is_weakcc = proc_info.IsWeakCC();
    bool is_weaknc = proc_info.IsWeakNC();

    if(!hitnucl) { assert(is_coh || is_imd || is_nuel); }
  
    // hit quark 
    // set only for DIS events
    int  qrk  = (is_dis) ? tgt.HitQrkPdg() : 0;     
    bool seaq = (is_dis) ? tgt.HitSeaQrk() : false; 

    // resonance id ($GENIE/src/BaryonResonance/BaryonResonance.h)
    // set only for resonance neutrinoproduction
    int resid = (is_res) ? xcls.Resonance() : -99;

    // (qel or dis) charm production?
    bool charm = xcls.IsCharmEvent();

    // get neut and nuance equivalent reaction codes (if any)
    brCodeNeut    = utils::ghep::NeutReactionCode(&event);
    brCodeNuance  = utils::ghep::NuanceReactionCode(&event);

    //weight
    double weight = event.Weight();

    // Access kinematical params _exactly_ as they were selected internally
    // (at the hit nucleon rest frame; 
    // for bound nucleons: taking into account fermi momentum and off-shell kinematics)
    //
    bool get_selected = true;
    double xs  = kine.x (get_selected);
    double ys  = kine.y (get_selected);
    double ts  = (is_coh) ? kine.t (get_selected) : -1;
    double Q2s = kine.Q2(get_selected);
    double Ws  = kine.W (get_selected);

    LOG("gntpc", pDEBUG) 
       << "[Select] Q2 = " << Q2s << ", W = " << Ws 
       << ", x = " << xs << ", y = " << ys << ", t = " << ts;

    // Calculate the same kinematical params but now as an experimentalist would 
    // measure them by neglecting the fermi momentum and off-shellness of bound nucleons
    //

    const TLorentzVector & k1 = *(neutrino->P4());                     // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());                          // l 4-p (k2)
    const TLorentzVector & p1 = (hitnucl) ? *(hitnucl->P4()) : pdummy; // N 4-p (p1)      

    double M  = kNucleonMass; 
    TLorentzVector q  = k1-k2;                     // q=k1-k2, 4-p transfer
    double Q2 = -1 * q.M2();                       // momemtum transfer
    double v  = (hitnucl) ? q.Energy()       : -1; // v (E transfer to the nucleus)
    double x  = (hitnucl) ? 0.5*Q2/(M*v)     : -1; // Bjorken x
    double y  = (hitnucl) ? v/k1.Energy()    : -1; // Inelasticity, y = q*P1/k1*P1
    double W2 = (hitnucl) ? M*M + 2*M*v - Q2 : -1; // Hadronic Invariant mass ^ 2
    double W  = (hitnucl) ? TMath::Sqrt(W2)  : -1; 
    double t  = 0;

    LOG("gntpc", pDEBUG) 
       << "[Calc] Q2 = " << Q2 << ", W = " << W 
       << ", x = " << x << ", y = " << y << ", t = " << t;

    // Extract more info on the hadronic system
    // Only for QEL/RES/DIS/COH events
    //
    bool study_hadsyst = (is_qel || is_res || is_dis || is_coh);
    
    //
    TObjArrayIter piter(&event);
    GHepParticle * p = 0;
    int ip=-1;

    //
    // Extract the final state system originating from the hadronic vertex 
    // (after the intranuclear rescattering step)
    //

    LOG("gntpc", pDEBUG) << "Extracting final state hadronic system";

    vector<int> final_had_syst;
    while( (p = (GHepParticle *) piter.Next()) && study_hadsyst)
    {
      ip++;
      // don't count final state lepton as part hadronic system 
      //if(!is_coh && event.Particle(ip)->FirstMother()==0) continue;
      if(event.Particle(ip)->FirstMother()==0) continue;
      if(p->IsFake()) continue;
      int pdgc = p->Pdg();
      int ist  = p->Status();
      if(ist==kIStStableFinalState) {
         if (pdgc == kPdgGamma || pdgc == kPdgElectron || pdgc == kPdgPositron)  {
            int igmom = p->FirstMother();
            if(igmom!=-1) {
	      // only count e+'s e-'s or gammas not from decay of pi0
	      if(event.Particle(igmom)->Pdg() != kPdgPi0) { final_had_syst.push_back(ip); }
            }
         } else {
            final_had_syst.push_back(ip);
         }
      }
      // now add pi0's that were decayed as short lived particles
      else if(pdgc == kPdgPi0){
	int ifd = p->FirstDaughter();
	int fd_pdgc = event.Particle(ifd)->Pdg();
	// just require that first daughter is one of gamma, e+ or e-  
	if(fd_pdgc == kPdgGamma || fd_pdgc == kPdgElectron || fd_pdgc == kPdgPositron){
	  final_had_syst.push_back(ip);
	}
      }
    }//particle-loop

    if( count(final_had_syst.begin(), final_had_syst.end(), -1) > 0) {
        mcrec->Clear();
 	continue;
    }

    //
    // Extract info on the primary hadronic system (before any intranuclear rescattering)
    // looking for particles with status_code == kIStHadronInTheNucleus 
    // An exception is the coherent production and scattering off free nucleon targets 
    // (no intranuclear rescattering) in which case primary hadronic system is set to be 
    // 'identical' with the final  state hadronic system
    //

    LOG("gntpc", pDEBUG) << "Extracting primary hadronic system";
    
    ip = -1;
    TObjArrayIter piter_prim(&event);

    vector<int> prim_had_syst;
    if(study_hadsyst) {
      // if coherent or free nucleon target set primary states equal to final states
      if(!target->IsNucleus() || (is_coh)) {
         vector<int>::const_iterator hiter = final_had_syst.begin();
         for( ; hiter != final_had_syst.end(); ++hiter) {
           prim_had_syst.push_back(*hiter);
         }
      } 
      // otherwise loop over all particles and store indices of those which are hadrons
      // created within the nucleus
      else {
	while( (p = (GHepParticle *) piter_prim.Next()) ){
	  ip++;      
	  int ist_comp  = p->Status();
	  if(ist_comp==kIStHadronInTheNucleus) {
	    prim_had_syst.push_back(ip); 
	  }
	}//particle-loop   
	//
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
    brIsCoh      = is_coh;  
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
    brPxv        = k1.Px();  
    brPyv        = k1.Py();  
    brPzv        = k1.Pz();  
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

    brCalResp0   = 0;

    brNf = final_had_syst.size();
    for(int j=0; j<brNf; j++) {
      p = event.Particle(final_had_syst[j]);
      assert(p);

      int    hpdg = p->Pdg();     
      double hE   = p->Energy();     
      double hKE  = p->KinE();     
      double hpx  = p->Px();     
      double hpy  = p->Py();     
      double hpz  = p->Pz();     
      double hm   = p->Mass();     

      brPdgf[j] = hpdg;
      brEf  [j] = hE;
      brPxf [j] = hpx;
      brPyf [j] = hpy;
      brPzf [j] = hpz;

      if      ( hpdg == kPdgProton      )  { brNfP++;     brCalResp0 += hKE;        }
      else if ( hpdg == kPdgAntiProton  )  { brNfP++;     brCalResp0 += (hE + hm);  }
      else if ( hpdg == kPdgNeutron     )  { brNfN++;     brCalResp0 += hKE;        }
      else if ( hpdg == kPdgAntiNeutron )  { brNfN++;     brCalResp0 += (hE + hm);  }
      else if ( hpdg == kPdgPiP         )  { brNfPip++;   brCalResp0 += hKE;        }
      else if ( hpdg == kPdgPiM         )  { brNfPim++;   brCalResp0 += hKE;        }
      else if ( hpdg == kPdgPi0         )  { brNfPi0++;   brCalResp0 += (e_h * hE); }
      else if ( hpdg == kPdgKP          )  { brNfKp++;    brCalResp0 += hKE;        }
      else if ( hpdg == kPdgKM          )  { brNfKm++;    brCalResp0 += hKE;        }
      else if ( hpdg == kPdgK0          )  { brNfK0++;    brCalResp0 += hKE;        }
      else if ( hpdg == kPdgAntiK0      )  { brNfK0++;    brCalResp0 += hKE;        }
      else if ( hpdg == kPdgGamma       )  { brNfEM++;    brCalResp0 += (e_h * hE); }
      else if ( hpdg == kPdgElectron    )  { brNfEM++;    brCalResp0 += (e_h * hE); }
      else if ( hpdg == kPdgPositron    )  { brNfEM++;    brCalResp0 += (e_h * hE); }
      else                                 { brNfOther++; brCalResp0 += hKE;        }

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
//  ** GENIE GHEP EVENT TREE FORMAT -> GENIE XML EVENT FILE FORMAT **
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
// ******** GENIE GHEP EVENT TREE FORMAT -> TRACKER FORMATS ********
//___________________________________________________________________
void ConvertToGTracker(void)
{
  //-- get pdglib
  PDGLibrary * pdglib = PDGLibrary::Instance();

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

  flux::GJPARCNuFluxPassThroughInfo * flux_info = 0;
  tree->SetBranchAddress("flux", &flux_info);

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
    int iparticle = -1;

    // **** Convert the current event:

    //
    // -- Add tracker begin tag
    //
    output << "$ begin" << endl;

    //
    // -- Add the appropriate reaction code
    //

    // add 'NEUT'-like event type
    if(gOptOutFileFormat == kConvFmt_t2k_tracker) {
    	int evtype = utils::ghep::NeutReactionCode(&event);
        LOG("gntpc", pNOTICE) << "NEUT-like event type = " << evtype;
    	output << "$ genie " << evtype << endl;
    } //neut code

    // add 'NUANCE'-like event type
    else if(gOptOutFileFormat == kConvFmt_nuance_tracker) {
    	int evtype = utils::ghep::NuanceReactionCode(&event);
        LOG("gntpc", pNOTICE) << "NUANCE-like event type = " << evtype;
    	output << "$ nuance " << evtype << endl;
    } // nuance code

    else {
        gAbortingInErr = true;
        exit(1);
    }

    //
    // -- Add '$vertex' line
    //
    output << "$ vertex " 
           << event.Vertex()->X() << " "
           << event.Vertex()->Y() << " "
           << event.Vertex()->Z() << " "
           << event.Vertex()->T() << endl;

    //
    // -- Add '$track' lines
    //

    // Loop over the generated GHEP particles and decide which ones 
    // to write-out in $track lines
    vector<int> tracks;

    event_iter.Reset();
    iparticle = -1;
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) 
    {
       iparticle++;

       // Neglect all GENIE pseudo-particles
       if(p->IsFake()) continue;

       int          ghep_pdgc   = p->Pdg();
       GHepStatus_t ghep_ist    = (GHepStatus_t) p->Status();

       //  
       // Keep 'initial state', 'nucleon target', 'hadron in the nucleus' and 'final state' particles.
       // Neglect pi0 decays if they were performed within GENIE (write out the decayed pi0 and neglect 
       // the {gamma + gamma} or {gamma + e- + e+} final state
       //

       // is pi0 decay?
       bool is_pi0_dec = false;
       if(ghep_ist == kIStDecayedState && ghep_pdgc == kPdgPi0) {
         vector<int> pi0dv; // daughters vector
         int ghep_fd = p->FirstDaughter();
         int ghep_ld = p->LastDaughter();
         for(int jd = ghep_fd; jd <= ghep_ld; jd++) {
           if(jd!=-1) {
              pi0dv.push_back(event.Particle(jd)->Pdg());
           }
         }
         sort(pi0dv.begin(), pi0dv.end());
         is_pi0_dec = (pi0dv.size()==2 && pi0dv[0]==kPdgGamma && pi0dv[1]==kPdgGamma) ||
                      (pi0dv.size()==3 && pi0dv[0]==kPdgPositron && pi0dv[1]==kPdgElectron && pi0dv[2]==kPdgGamma);
       }

       // is pi0 decay product?
       int ghep_fm     = p->FirstMother();
       int ghep_fmpdgc = (ghep_fm==-1) ? 0 : event.Particle(ghep_fm)->Pdg();
       bool is_pi0_dpro = (ghep_pdgc == kPdgGamma    && ghep_fmpdgc == kPdgPi0) ||
                          (ghep_pdgc == kPdgElectron && ghep_fmpdgc == kPdgPi0) ||
                          (ghep_pdgc == kPdgPositron && ghep_fmpdgc == kPdgPi0);

       bool keep = (ghep_ist == kIStInitialState)       ||
                   (ghep_ist == kIStNucleonTarget)      ||
                   (ghep_ist == kIStHadronInTheNucleus) ||
                   (ghep_ist == kIStDecayedState     &&  is_pi0_dec ) ||
                   (ghep_ist == kIStStableFinalState && !is_pi0_dpro);
       if(!keep) continue;

       // Apparently SKDETSIM chokes with O16 - Neglect the nuclear target in this case
       //
       if (gOptOutFileFormat == kConvFmt_t2k_tracker && p->IsNucleus()) continue;

       tracks.push_back(iparticle);
    }

    //bool info_added  = false;

    // Looping twice to ensure that all final state particle are grouped together.
    // On the second loop add only f/s particles. On the first loop add all but f/s particles
    for(int iloop=0; iloop<=1; iloop++) 
    {
      for(vector<int>::const_iterator ip = tracks.begin(); ip != tracks.end(); ++ip) 
      {
         iparticle = *ip;
         p = event.Particle(iparticle);

         int ghep_pdgc = p->Pdg();
         GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();

         bool fs = (ghep_ist==kIStStableFinalState) || 
                   (ghep_ist==kIStDecayedState && ghep_pdgc==kPdgPi0);

         if(iloop==0 &&  fs) continue;
         if(iloop==1 && !fs) continue;

         // Convert GENIE's GHEP pdgc & status to NUANCE's equivalent
         //
         int ist;
         switch (ghep_ist) {
           case kIStInitialState:             ist = -1;                              break;
           case kIStStableFinalState:         ist =  0;                              break;
           case kIStIntermediateState:        ist = -2;                              break;
           case kIStDecayedState:             ist = (ghep_pdgc==kPdgPi0) ? 0 : -2;   break;
           case kIStNucleonTarget:            ist = -1;                              break;
           case kIStDISPreFragmHadronicState: ist = -999;                            break;
           case kIStPreDecayResonantState:    ist = -999;                            break;
           case kIStHadronInTheNucleus:       ist = -2;                              break;
           case kIStUndefined:                ist = -999;                            break;
           default:                           ist = -999;                            break;
         }
         // Convert GENIE pdg code -> nuance PDG code
         // For most particles both generators use the standard PDG codes.
         // For nuclei GENIE follows the PDG-convention: 10LZZZAAAI
         // NUANCE is using: ZZZAAA
         int pdgc = ghep_pdgc;
         if ( p->IsNucleus() ) {
           int Z = pdg::IonPdgCodeToZ(ghep_pdgc);
           int A = pdg::IonPdgCodeToA(ghep_pdgc);
           pdgc = 1000*Z + A;
         }

         // The SK detector MC expects K0_Long, K0_Short - not K0, \bar{K0}
         // Do the conversion here:
         if(gOptOutFileFormat == kConvFmt_t2k_tracker) {
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

// <obsolte/>
//         GHepStatus_t gist = (GHepStatus_t) p->Status();
//         bool is_init =
//               (gist == kIStInitialState || gist == kIStNucleonTarget);
//
//         if(!is_init && !info_added) {
//           // Add nuance obsolete and flux info (not filled in by
//           // GENIE here). Add it once after the initial state particles
//           output << "$ info 2 949000 0.0000E+00" << endl;
//           info_added = true;
//         }
// </obsolte>

         LOG("gntpc", pNOTICE) 
           << "Adding $track corrsponding to GHEP particle at position: " << iparticle
           << " (tracker status code: " << ist << ")";

         output << "$ track " << pdgc << " " << E << " "
                << dcosx << " " << dcosy << " " << dcosz << " "
                << ist << endl;

      }//tracks
    }//iloop

    //
    // -- Add $info lines as necessary
    //

    if(gOptOutFileFormat == kConvFmt_t2k_tracker) {
      //
      // Writing $info lines with information identical to the one saved at the rootracker-format 
      // files for the nd280MC. SKDETSIM can propagate all that complete MC truth information into 
      // friend event trees that can be 'linked' with the SK DSTs.
      // Having identical generator info for both SK and nd280 will enable global studies
      //
      // The $info lines are formatted as follows:
      //
      // $ info event_num err_flag string_event_code
      // $ info xsec_event diff_xsec_kinematics weight prob
      // $ info vtxx vtxy vtxz vtxt
      // $ info nparticles
      // $ info 0 pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz
      // $ info 1 pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz
      // ... ... ...
      // $ info k pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz
      // ... ... ...
      // $ info jnubeam_parent_pdg_code jnubeam_parent_decay_mode
      // $ info jnubeam_parent_decay_px jnubeam_parent_decay_py jnubeam_parent_decay_pz jnubeam_parent_decay_E
      // $ info jnubeam_parent_decay_x jnubeam_parent_decay_y jnubeam_parent_decay_z jnubeam_parent_decay_t
      // $ info jnubeam_parent_prod_px jnubeam_parent_prod_py jnubeam_parent_prod_pz jnubeam_parent_prod_E
      // $ info jnubeam_parent_prod_x jnubeam_parent_prod_y jnubeam_parent_prod_z jnubeam_parent_prod_t
      // $ info jnubeam_parent_nvtx
      //
      // Comments:
      // - The jnubeam lines may not always be available (eg if event generation used histogram-based flux descriptions)
      // - The jnubeam variables are in whatever units are used by jnubeam.
      // - The err_flag is a bit field (16 bits)
      // - The string_event_code is a rather long string which encapsulates lot of summary info on the event
      //   (neutrino/nuclear target/hit nucleon/hit quark(if any)/process type/...).
      //   Information on how to parse that string code is available at the T2K event reweighting package.
      // - event_xsec is the event cross section in 1E-38cm^2
      // - diff_event_xsec is the cross section for the selected in 1E-38cm^2/{K^n}
      // - weight is the event weight (1 for unweighted MC)
      // - prob is the event probability (given cross sectios and density-weighted path-length)
      // - vtxx,y,z,t is the vertex position/time in SI units 
      // - nparticles is the number of particles in the GHEP record (number of $info lines to follow before the start of the JNUBEAM block)
      // - first_/last_daughter first_/last_mother indicate the particle
      // - px,py,pz,E is the particle 4-momentum at the LAB frame (in GeV)
      // - x,y,z,t is the particle 4-position at the hit nucleus coordinate system (in fm, t is not set)
      // - polx,y,z is the particle polarization vector
      // See also ConvertToGRooTracker() for further descriptions of the variables stored at
      // the rootracker files.
      //
      // event info
      //
      output << "$ info " << (int) iev << " " << *(event.EventFlags()) << " " << interaction->AsString() << endl;
      output << "$ info " << (1E+38/units::cm2) * event.XSec() << " "
                          << (1E+38/units::cm2) * event.DiffXSec() << " "  
                          << event.Weight() << " "
                          << event.Probability()
                          << endl;
      output << "$ info " << event.Vertex()->X() << " "
                          << event.Vertex()->Y() << " "
                          << event.Vertex()->Z() << " "
                          << event.Vertex()->T() 
                          << endl;
      //
      // copy stdhep-like particle list
      //
      iparticle = 0;
      event_iter.Reset();
      output << "$ info " << event.GetEntries() << endl;
      while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) 
      {
        assert(p);
        output << "$ info " 
               << iparticle << " " 
               << p->Pdg() << " " << (int) p->Status() << " "
               << p->FirstDaughter() << " " << p->LastDaughter() << " " 
               << p->FirstMother() << " " << p->LastMother() << " "
               << p->X4()->X()  << " " << p->X4()->Y()  << " " << p->X4()->Z()  << " " << p->X4()->T() << " "
               << p->P4()->Px() << " " << p->P4()->Py() << " " << p->P4()->Pz() << " " << p->P4()->E() << " ";
        if(p->PolzIsSet()) {
            output << TMath::Sin(p->PolzPolarAngle()) * TMath::Cos(p->PolzAzimuthAngle()) << " "
                   << TMath::Sin(p->PolzPolarAngle()) * TMath::Sin(p->PolzAzimuthAngle()) << " "
                   << TMath::Cos(p->PolzPolarAngle());
        } else {
            output << "0. 0. 0.";
        }
        output << endl;
        iparticle++;
      }
      //
      // JNUBEAM flux info - that may not be available, eg if events were generated using 
      // plain flux histograms - not the beam simulation's output flux ntuples
      //
      if(flux_info) {
         // parent hadron pdg code and decay mode
         output << "$ info " << flux_info->pdg << " " << flux_info->decayMode << endl;
         // parent hadron px,py,pz,E at decay
         output << "$ info " << flux_info->decayP * flux_info->decayDirX << " " 
                             << flux_info->decayP * flux_info->decayDirY << " " 
                             << flux_info->decayP * flux_info->decayDirZ << " " 
                             << TMath::Sqrt(
                                   TMath::Power(pdglib->Find(flux_info->pdg)->Mass(), 2.)
                                 + TMath::Power(flux_info->decayP, 2.)
                                )  << endl;
         // parent hadron x,y,z,t at decay
         output << "$ info " << flux_info->decayX << " "
                             << flux_info->decayY << " "
                             << flux_info->decayZ << " "
                             << "0." 
                             << endl;
         // parent hadron px,py,pz,E at production
         output << "$ info " << flux_info->prodP * flux_info->prodDirX << " "
                             << flux_info->prodP * flux_info->prodDirY << " "
                             << flux_info->prodP * flux_info->prodDirZ << " "
                             << TMath::Sqrt(
                                   TMath::Power(pdglib->Find(flux_info->pdg)->Mass(), 2.)
                                 + TMath::Power(flux_info->prodP, 2.)
                                ) << endl;
         // parent hadron x,y,z,t at production
         output << "$ info " << flux_info->prodX << " "
                             << flux_info->prodY << " "
                             << flux_info->prodZ << " "
                             << "0." 
                             << endl;
         // nvtx
         output << "$ info " << output << "$info " << endl;
     }
    }//fmt==kConvFmt_t2k_tracker

    //
    // -- Add  tracker end tag
    //
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
// **** GENIE GHEP EVENT TREE FORMAT -> A T2K ROOTRACKER FORMAT ****
//___________________________________________________________________
void ConvertToGRooTracker(void)
{
  //-- get pdglib
  PDGLibrary * pdglib = PDGLibrary::Instance();

  //-- define the output rootracker tree branches
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
  // > etc 
  int         brNeutCode;                 // > NEUT-like reaction code for the GENIE event

  //-- open the output ROOT file
  TFile fout(gOptOutFileName.c_str(), "RECREATE");

  //-- create the output ROOT tree
  TTree * rootracker_tree = new TTree("gRooTracker","GENIE event tree rootracker format");

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
  if(gOptOutFileFormat == kConvFmt_t2k_rootracker) {
    // JNUBEAM pass-through info -- available ony on the t2k version of the rootracker format
    rootracker_tree->Branch("NuParentPdg",     &brNuParentPdg,     "NuParentPdg/I");       
    rootracker_tree->Branch("NuParentDecMode", &brNuParentDecMode, "NuParentDecMode/I");   
    rootracker_tree->Branch("NuParentDecP4",    brNuParentDecP4,   "NuParentDecP4[4]/D");     
    rootracker_tree->Branch("NuParentDecX4",    brNuParentDecX4,   "NuParentDecX4[4]/D");     
    rootracker_tree->Branch("NuParentProP4",    brNuParentProP4,   "NuParentProP4[4]/D");     
    rootracker_tree->Branch("NuParentProX4",    brNuParentProX4,   "NuParentProX4[4]/D");     
    rootracker_tree->Branch("NuParentProNVtx", &brNuParentProNVtx, "NuParentProNVtx/I");   
    // NEUT-like reaction code -- available only on the t2k version of the rootracker format
    rootracker_tree->Branch("G2NeutEvtCode",   &brNeutCode,        "G2NeutEvtCode/I");   
  }

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

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
  flux::GJPARCNuFluxPassThroughInfo * flux_info = 0;
  if(gOptOutFileFormat == kConvFmt_t2k_rootracker) {
     gtree->SetBranchAddress("flux", &flux_info);
  }
#else
  LOG("gntpc", pWARN) 
    << "\n Flux drivers are not enabled." 
    << "\n No flux pass-through information will be written-out in the rootracker file"
    << "\n If this isn't what you are supposed to be doing then build GENIE by adding "
    << "--with-flux-drivers in the configuration step.";
#endif

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
    if(gOptOutFileFormat == kConvFmt_t2k_rootracker) {
       if(flux_info) {
          LOG("gntpc", pINFO) << *flux_info;
       } else {
          LOG("gntpc", pINFO) << "No JNUBEAM flux info associated with this event";
       }
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
    brNeutCode = 0;     

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

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
    // Copy flux info if this is the t2k rootracker variance.
    // The flux may not be available, eg if events were generated using plain flux 
    // histograms and not the JNUBEAM simulation's output flux ntuples.
    if(gOptOutFileFormat == kConvFmt_t2k_rootracker) {
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
    }
#endif

    // map GENIE event to NEUT reaction codes
    brNeutCode = utils::ghep::NeutReactionCode(&event);

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
  //get input ROOT file (containing a GENIE GHEP event tree)
  try {
    LOG("gntpc", pINFO) << "Reading input filename";
    gOptInpFileName = utils::clap::CmdLineArgAsString(argc,argv,'i');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gntpc", pFATAL)
               << "Unspecified input filename - Exiting";
      PrintSyntax();
      gAbortingInErr = true;
      exit(1);
    }
  }

  //check input GENIE ROOT file
  bool inpok = !(gSystem->AccessPathName(gOptInpFileName.c_str()));
  if (!inpok) {
    LOG("gntpc", pFATAL)
           << "The input ROOT file ["
                       << gOptInpFileName << "] is not accessible";
    gAbortingInErr = true;
    exit(2);
  }

  //get output file format
  try {
    LOG("gntpc", pINFO) << "Reading output file format";
    string fmt = utils::clap::CmdLineArgAsString(argc,argv,'f');

         if (fmt == "t2k_rootracker") { gOptOutFileFormat = kConvFmt_t2k_rootracker; }
    else if (fmt == "t2k_tracker")    { gOptOutFileFormat = kConvFmt_t2k_tracker;    }
    else if (fmt == "gst")            { gOptOutFileFormat = kConvFmt_gst;            }
    else if (fmt == "rootracker")     { gOptOutFileFormat = kConvFmt_rootracker;     }
    else if (fmt == "gxml")           { gOptOutFileFormat = kConvFmt_gxml;           }
    else if (fmt == "ghad")           { gOptOutFileFormat = kConvFmt_ghad;           }
    else if (fmt == "ginuke")         { gOptOutFileFormat = kConvFmt_ginuke;         }
    else if (fmt == "nuance_tracker") { gOptOutFileFormat = kConvFmt_nuance_tracker; }
    else                              { gOptOutFileFormat = kConvFmt_undef;          }

    if(gOptOutFileFormat == kConvFmt_undef) {
      LOG("gntpc", pFATAL) << "Unknown output file format (" << fmt << ")";
      gAbortingInErr = true;
      exit(3);
    }

  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gntpc", pFATAL) << "Unspecified output file format";
      gAbortingInErr = true;
      exit(4);
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
  if      (gOptOutFileFormat == kConvFmt_t2k_rootracker) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_t2k_tracker   ) { ext = "gtrac.dat";        }
  else if (gOptOutFileFormat == kConvFmt_gst           ) { ext = "gst.root";         }
  else if (gOptOutFileFormat == kConvFmt_rootracker    ) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_gxml          ) { ext = "gxml";             }
  else if (gOptOutFileFormat == kConvFmt_ghad          ) { ext = "ghad.dat";         }
  else if (gOptOutFileFormat == kConvFmt_ginuke        ) { ext = "ginuke.root";      }
  else if (gOptOutFileFormat == kConvFmt_nuance_tracker) { ext = "gtrac_legacy.dat"; }

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
