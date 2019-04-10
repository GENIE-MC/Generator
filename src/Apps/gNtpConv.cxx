//_____________________________________________________________________________________________
/*!

\program gntpc

\brief   Converts a native GENIE (GHEP/ROOT) event tree file to a host of
         plain text, XML or bare-ROOT formats.

         Syntax:
           gntpc -i input_file [-o output_file] -f format [-n nev] [-v vrs] [-c] 
                 [--seed random_number_seed]
                 [--message-thresholds xml_file]
                 [--event-record-print-level level]


         Options :

           [] denotes an optional argument

           -n 
              Number of events to convert 
              (optional, default: convert all events)
           -v 
              Output format version, if multiple versions are supported
              (optional, default: use latest version of each format)
           -c 
              Copy MC job metadata (gconfig and genv TFolders) from the input GHEP file.
           -f 
              A string that specifies the output file format. 
              >>
	      >> Generic formats:
              >>
               * `gst': 
                    The 'definite' GENIE summary tree format (gst).
   	       * `gxml': 
                     GENIE XML event format 
   	       * `ghep_mock_data': 
                     Output file has the same format as the input file (GHEP) but
                     all information other than final state particles is hidden
   	       * `rootracker': 
                     A bare-ROOT STDHEP-like GENIE event tree.
   	       * `rootracker_mock_data': 
                     As the `rootracker' format but hiddes all information
                     except the final state particles.
              >>
	      >> Experiment-specific formats:
              >>
   	       * `t2k_rootracker':
                     A variance of the `rootracker' format used by the nd280, INGRID and 2km. 
                     Includes, in addition, JPARC flux pass-through info.
   	       * `numi_rootracker':
                     A variance of the `rootracker' format for the NuMI expts.
                     Includes, in addition, NuMI flux pass-through info.
   	       * `t2k_tracker': 
                     A tracker-type format with tweaks required by the SuperK MC (SKDETSIM):
                        - Converting K0, \bar{K0} -> KO_{long}, K0_{short}
                        - Emulating 'NEUT' reaction codes
                        - Appropriate $track ordering for SKDETSIM
                        - Passing detailed GENIE MC truth and JPARC flux info
                          using the tracker $info lines. This information, 
                          propaged by SKDETSIM to the DSTs, is identical with the 
                          one used at the near detectors and can be used for 
                          global systematic studies.
              >>
	      >> GENIE test / cross-generator comparison formats:
              >>
   	       * `ghad': 
	             NEUGEN-style text-based format for hadronization studies
   	       * `ginuke': 
	             A summary ntuple for intranuclear-rescattering studies using simulated
                     hadron-nucleus samples
              >>
	      >> Other (depreciated) formats:
              >>
   	       * `nuance_tracker': 
   		     NUANCE-style tracker text-based format 
           -o  
              Specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
               `gst'                  -> *.gst.root
               `gxml'                 -> *.gxml 
               `ghep_mock_data'       -> *.mockd.ghep.root
               `rootracker'           -> *.gtrac.root
               `rootracker_mock_data' -> *.mockd.gtrac.root
               `t2k_rootracker'       -> *.gtrac.root
               `numi_rootracker'      -> *.gtrac.root
               `t2k_tracker'          -> *.gtrac.dat
               `nuance_tracker'       -> *.gtrac_legacy.dat
               `ghad'                 -> *.ghad.dat
               `ginuke'               -> *.ginuke.root
           --seed
              Random number seed.
         --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
           --event-record-print-level
              Allows users to set the level of information shown when the event
              record is printed in the screen. See GHepRecord::Print().
		
         Examples:
           (1)  shell% gntpc -i myfile.ghep.root -f t2k_rootracker

                Converts all events in the GHEP file myfile.ghep.root into the
                t2k_rootracker format. 
                The output file is named myfile.gtrac.root

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created September 23, 2005

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_____________________________________________________________________________________________

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
#include <TFolder.h>
#include <TBits.h>
#include <TObjString.h>
#include <TMath.h>
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/T2KEvGenMetaData.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#include "Tools/Flux/GJPARCNuFlux.h"
#include "Tools/Flux/GNuMIFlux.h"
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
using std::setiosflags;
using std::vector;

using namespace genie;
using namespace genie::constants;

//func prototypes
void   ConvertToGST              (void);
void   ConvertToGXML             (void);
void   ConvertToGHepMock         (void);
void   ConvertToGTracker         (void);
void   ConvertToGRooTracker      (void);
void   ConvertToGHad             (void);
void   ConvertToGINuke           (void);
void   GetCommandLineArgs        (int argc, char ** argv);
void   PrintSyntax               (void);
string DefaultOutputFile         (void);
int    LatestFormatVersionNumber (void);
bool   CheckRootFilename         (string filename);
int    HAProbeFSI                (int, int, int, double [], int [], int, int, int); //Test code
//format enum
typedef enum EGNtpcFmt {
  kConvFmt_undef = 0,
  kConvFmt_gst,
  kConvFmt_gxml,
  kConvFmt_ghep_mock_data,
  kConvFmt_rootracker,
  kConvFmt_rootracker_mock_data,
  kConvFmt_t2k_rootracker,
  kConvFmt_numi_rootracker,
  kConvFmt_t2k_tracker,
  kConvFmt_nuance_tracker,
  kConvFmt_ghad,
  kConvFmt_ginuke
} GNtpcFmt_t;

//input options (from command line arguments):
string     gOptInpFileName;         ///< input file name
string     gOptOutFileName;         ///< output file name
GNtpcFmt_t gOptOutFileFormat;       ///< output file format id
int        gOptVersion;             ///< output file format version
Long64_t   gOptN;                   ///< number of events to process
bool       gOptCopyJobMeta = false; ///< copy MC job metadata (gconfig, genv TFolders)
long int   gOptRanSeed;             ///< random number seed

//genie version used to generate the input event file 
int gFileMajorVrs = -1;
int gFileMinorVrs = -1;
int gFileRevisVrs = -1;

//consts
const int kNPmax = 250;
//____________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);

  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);

  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  PDGLibrary::Instance()->AddDarkMatter( 1.0, 0.5 ) ;
  
  // Call the appropriate conversion function
  switch(gOptOutFileFormat) {

   case (kConvFmt_gst)  :

	ConvertToGST();        
	break;  

   case (kConvFmt_gxml) :  

	ConvertToGXML();         
	break;

   case (kConvFmt_ghep_mock_data) :  

	ConvertToGHepMock();         
	break;

   case (kConvFmt_rootracker          ) :  
   case (kConvFmt_rootracker_mock_data) :  
   case (kConvFmt_t2k_rootracker      ) :  
   case (kConvFmt_numi_rootracker     ) :  

	ConvertToGRooTracker(); 
	break;

   case (kConvFmt_t2k_tracker   )  :  
   case (kConvFmt_nuance_tracker)  :  

	ConvertToGTracker();        
	break;

   case (kConvFmt_ghad) :  

	ConvertToGHad();         
	break;

   case (kConvFmt_ginuke) :  

	ConvertToGINuke();         
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
//____________________________________________________________________________________
// GENIE GHEP EVENT TREE FORMAT -> GENIE SUMMARY NTUPLE 
//____________________________________________________________________________________
void ConvertToGST(void)
{
  // Some constants
  const double e_h = 1.3; // typical e/h ratio used for computing mean `calorimetric response'

  // Define branch variables
  //
  int    brIev         = 0;      // Event number 
  int    brNeutrino    = 0;      // Neutrino pdg code
  int    brFSPrimLept  = 0;      // Final state primary lepton pdg code
  int    brTarget      = 0;      // Nuclear target pdg code (10LZZZAAAI)
  int    brTargetZ     = 0;      // Nuclear target Z (extracted from pdg code above)
  int    brTargetA     = 0;      // Nuclear target A (extracted from pdg code above)
  int    brHitNuc      = 0;      // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
  int    brHitQrk      = 0;      // Hit quark pdg code        (set for DIS events only)
  bool   brFromSea     = false;  // Hit quark is from sea     (set for DIS events only)
  int    brResId       = 0;      // Produced baryon resonance (set for resonance events only)
  bool   brIsQel       = false;  // Is QEL?
  bool   brIsRes       = false;  // Is RES?
  bool   brIsDis       = false;  // Is DIS?
  bool   brIsCoh       = false;  // Is Coherent?
  bool   brIsMec       = false;  // Is MEC?
  bool   brIsDfr       = false;  // Is Diffractive?
  bool   brIsImd       = false;  // Is IMD?
  bool   brIsSingleK   = false;  // Is single kaon?  
  bool   brIsImdAnh    = false;  // Is IMD annihilation?
  bool   brIsNuEL      = false;  // Is ve elastic?
  bool   brIsEM        = false;  // Is EM process?
  bool   brIsCC        = false;  // Is Weak CC process?
  bool   brIsNC        = false;  // Is Weak NC process?
  bool   brIsCharmPro  = false;  // Produces charm?
  bool   brIsAMNuGamma = false;  // is anomaly mediated nu gamma
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
  double brEvRF        = 0;      // Neutrino energy @ the rest-frame of the hit-object (eg nucleon for CCQE, e- for ve- elastic,...)
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
  double brPzl         = 0;      // Final state primary lepton pz @ LAB
  double brPl          = 0;      // Final state primary lepton p  @ LAB
  double brCosthl      = 0;      // Final state primary lepton cos(theta) wrt to neutrino direction
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
  int    brPdgf  [kNPmax];       // Pdg code of k^th final state particle in hadronic system
  double brEf    [kNPmax];       // Energy     of k^th final state particle in hadronic system @ LAB
  double brPxf   [kNPmax];       // Px         of k^th final state particle in hadronic system @ LAB
  double brPyf   [kNPmax];       // Py         of k^th final state particle in hadronic system @ LAB
  double brPzf   [kNPmax];       // Pz         of k^th final state particle in hadronic system @ LAB
  double brPf    [kNPmax];       // P          of k^th final state particle in hadronic system @ LAB
  double brCosthf[kNPmax];       // cos(theta) of k^th final state particle in hadronic system @ LAB wrt to neutrino direction
  int    brNi          = 0;      // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
  int    brPdgi[kNPmax];         // Pdg code of k^th particle in 'primary' hadronic system 
  int    brResc[kNPmax];         // FSI code of k^th particle in 'primary' hadronic system 
  double brEi  [kNPmax];         // Energy   of k^th particle in 'primary' hadronic system @ LAB
  double brPxi [kNPmax];         // Px       of k^th particle in 'primary' hadronic system @ LAB
  double brPyi [kNPmax];         // Py       of k^th particle in 'primary' hadronic system @ LAB
  double brPzi [kNPmax];         // Pz       of k^th particle in 'primary' hadronic system @ LAB
  double brVtxX;                 // Vertex x in detector coord system (SI)
  double brVtxY;                 // Vertex y in detector coord system (SI)
  double brVtxZ;                 // Vertex z in detector coord system (SI)
  double brVtxT;                 // Vertex t in detector coord system (SI)
  double brSumKEf;               // Sum of kinetic energies of all final state particles
  double brCalResp0;             // Approximate calorimetric response to the hadronic system computed as sum of
				 //  - (kinetic energy) for pi+, pi-, p, n 
                                 //  - (energy + 2*mass) for antiproton, antineutron
                                 //  - ((e/h) * energy)   for pi0, gamma, e-, e+, where e/h is set to 1.3
                                 //  - (kinetic energy) for other particles

  // Open output file & create output summary tree & create the tree branches
  //
  LOG("gntpc", pNOTICE) 
       << "*** Saving summary tree to: " << gOptOutFileName;
  TFile fout(gOptOutFileName.c_str(),"recreate");

  TTree * s_tree = new TTree("gst","GENIE Summary Event Tree");

  // Create tree branches
  //
  s_tree->Branch("iev",           &brIev,           "iev/I"         );
  s_tree->Branch("neu",	          &brNeutrino,      "neu/I"	    );
  s_tree->Branch("fspl",	      &brFSPrimLept,    "fspl/I"	    );
  s_tree->Branch("tgt",           &brTarget,        "tgt/I"	    );
  s_tree->Branch("Z",             &brTargetZ,       "Z/I"	    );
  s_tree->Branch("A",             &brTargetA,       "A/I"	    );
  s_tree->Branch("hitnuc",        &brHitNuc,        "hitnuc/I"      );
  s_tree->Branch("hitqrk",        &brHitQrk,        "hitqrk/I"      );
  s_tree->Branch("resid",         &brResId,	    "resid/I"	    );
  s_tree->Branch("sea",	          &brFromSea,       "sea/O"	    );
  s_tree->Branch("qel",	          &brIsQel,	    "qel/O"	    );
  s_tree->Branch("mec",	          &brIsMec,	    "mec/O"	    );
  s_tree->Branch("res",	          &brIsRes,	    "res/O"	    );
  s_tree->Branch("dis",	          &brIsDis,	    "dis/O"	    );
  s_tree->Branch("coh",           &brIsCoh,         "coh/O"	    );
  s_tree->Branch("dfr",           &brIsDfr,         "dfr/O"	    );
  s_tree->Branch("imd",	          &brIsImd,	    "imd/O"	    );
  s_tree->Branch("imdanh",        &brIsImdAnh,	    "imdanh/O"	    );
  s_tree->Branch("singlek",       &brIsSingleK,     "singlek/O"     );  
  s_tree->Branch("nuel",          &brIsNuEL,        "nuel/O"	    );
  s_tree->Branch("em",	          &brIsEM,	    "em/O"	    );
  s_tree->Branch("cc",	          &brIsCC,	    "cc/O"	    );
  s_tree->Branch("nc",	          &brIsNC,	    "nc/O"	    );
  s_tree->Branch("charm",         &brIsCharmPro,    "charm/O"	    );
  s_tree->Branch("amnugamma",     &brIsAMNuGamma,    "amnugamma/O"	    );
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
  s_tree->Branch("EvRF",	      &brEvRF,	    "EvRF/D"	    );
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
  s_tree->Branch("pl",            &brPl,            "pl/D"          );
  s_tree->Branch("cthl",          &brCosthl,        "cthl/D"        );
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
  s_tree->Branch("nip",	          &brNiP,	    "nip/I"	    );
  s_tree->Branch("nin",	          &brNiN,	    "nin/I"	    );
  s_tree->Branch("nipip",         &brNiPip,	    "nipip/I"	    );
  s_tree->Branch("nipim",         &brNiPim,	    "nipim/I"	    );
  s_tree->Branch("nipi0",         &brNiPi0,	    "nipi0/I"	    );
  s_tree->Branch("nikp",          &brNiKp,	    "nikp/I"	    );
  s_tree->Branch("nikm",          &brNiKm,	    "nikm/I"	    );
  s_tree->Branch("nik0",          &brNiK0,	    "nik0/I"	    );
  s_tree->Branch("niem",          &brNiEM,	    "niem/I"	    );
  s_tree->Branch("niother",       &brNiOther,       "niother/I"     );
  s_tree->Branch("ni",	         &brNi,	            "ni/I"	    );
  s_tree->Branch("pdgi",          brPdgi,	    "pdgi[ni]/I "   );
  s_tree->Branch("resc",          brResc,	    "resc[ni]/I "   );
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
  s_tree->Branch("pf",            brPf,             "pf[nf]/D"      );
  s_tree->Branch("cthf",          brCosthf,         "cthf[nf]/D"    );
  s_tree->Branch("vtxx",         &brVtxX,	    "vtxx/D"        );
  s_tree->Branch("vtxy",         &brVtxY,	    "vtxy/D"        );
  s_tree->Branch("vtxz",         &brVtxZ,	    "vtxz/D"        );
  s_tree->Branch("vtxt",         &brVtxT,	    "vtxt/D"        );
  s_tree->Branch("sumKEf",       &brSumKEf,	    "sumKEf/D"      );
  s_tree->Branch("calresp0",     &brCalResp0,	    "calresp0/D"    );

  // Open the ROOT file and get the TTree & its header
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

  // Get the mc record
  NtpMCEventRecord * mcrec = 0;
  er_tree->SetBranchAddress("gmcrec", &mcrec);
  if (!mcrec) {
    LOG("gntpc", pERROR) << "Null MC record";
    return;
  }
  
  // Figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ? 
       er_tree->GetEntries() : TMath::Min( er_tree->GetEntries(), gOptN );
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }

  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  TLorentzVector pdummy(0,0,0,0);

  // Event loop
  for(Long64_t iev = 0; iev < nmax; iev++) {
    er_tree->GetEntry(iev);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    // Go further only if the event is physical
    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) {
      LOG("gntpc", pINFO) << "Skipping unphysical event";
      mcrec->Clear();
      continue;
    }

    // Clean-up arrays
    //
    for(int j=0; j<kNPmax; j++) {
       brPdgi   [j] =  0;     
       brResc   [j] = -1;     
       brEi     [j] =  0;     
       brPxi    [j] =  0;     
       brPyi    [j] =  0;     
       brPzi    [j] =  0;     
       brPdgf   [j] =  0;     
       brEf     [j] =  0;     
       brPxf    [j] =  0;     
       brPyf    [j] =  0;     
       brPzf    [j] =  0;     
       brPf     [j] =  0;     
       brCosthf [j] =  0;     
    }

    // Computing event characteristics
    //

    //input particles
    GHepParticle * neutrino = event.Probe();
    GHepParticle * target = event.Particle(1);
    assert(target);
    GHepParticle * fsl = event.FinalStatePrimaryLepton();
    GHepParticle * hitnucl = event.HitNucleon();

    int tgtZ = 0;
    int tgtA = 0;
    if(pdg::IsIon(target->Pdg())) {
       tgtZ = pdg::IonPdgCodeToZ(target->Pdg());
       tgtA = pdg::IonPdgCodeToA(target->Pdg());
    } 
    if(target->Pdg() == kPdgProton   ) { tgtZ = 1; tgtA = 1; }    
    if(target->Pdg() == kPdgNeutron  ) { tgtZ = 0; tgtA = 1; }    
  
    // Summary info
    const Interaction * interaction = event.Summary();
    const InitialState & init_state = interaction->InitState();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    const Kinematics &   kine       = interaction->Kine();
    const XclsTag &      xcls       = interaction->ExclTag();
    const Target &       tgt        = init_state.Tgt();

    // Vertex in detector coord system
    TLorentzVector * vtx = event.Vertex();

    // Process id
    bool is_qel    = proc_info.IsQuasiElastic();
    bool is_res    = proc_info.IsResonant();
    bool is_dis    = proc_info.IsDeepInelastic();
    bool is_coh    = proc_info.IsCoherent();
    bool is_dfr    = proc_info.IsDiffractive();
    bool is_imd    = proc_info.IsInverseMuDecay();
    bool is_imdanh = proc_info.IsIMDAnnihilation();
    bool is_singlek = proc_info.IsSingleKaon();    
    bool is_nuel      = proc_info.IsNuElectronElastic();
    bool is_em        = proc_info.IsEM();
    bool is_weakcc    = proc_info.IsWeakCC();
    bool is_weaknc    = proc_info.IsWeakNC();
    bool is_mec       = proc_info.IsMEC();
    bool is_amnugamma = proc_info.IsAMNuGamma();

    if (!hitnucl && neutrino) {
        assert(is_coh || is_imd || is_imdanh || is_nuel | is_amnugamma);
    }
  
    // Hit quark - set only for DIS events
    int  qrk  = (is_dis) ? tgt.HitQrkPdg() : 0;     
    bool seaq = (is_dis) ? tgt.HitSeaQrk() : false; 

    // Resonance id ($GENIE/src/BaryonResonance/BaryonResonance.h) -
    // set only for resonance neutrinoproduction
    int resid = (is_res) ? EResonance(xcls.Resonance()) : -99;

    // (qel or dis) charm production?
    bool charm = xcls.IsCharmEvent();

    // Get NEUT and NUANCE equivalent reaction codes (if any)
    brCodeNeut    = utils::ghep::NeutReactionCode(&event);
    brCodeNuance  = utils::ghep::NuanceReactionCode(&event);

    // Get event weight
    double weight = event.Weight();

    // Access kinematical params _exactly_ as they were selected internally
    // (at the hit nucleon rest frame; 
    // for bound nucleons: taking into account fermi momentum and off-shell kinematics)
    //
    bool get_selected = true;
    double xs  = kine.x (get_selected);
    double ys  = kine.y (get_selected);
    double ts  = (is_coh || is_dfr) ? kine.t (get_selected) : -1;
    double Q2s = kine.Q2(get_selected);
    double Ws  = kine.W (get_selected);

    LOG("gntpc", pDEBUG) 
       << "[Select] Q2 = " << Q2s << ", W = " << Ws 
       << ", x = " << xs << ", y = " << ys << ", t = " << ts;

    // Calculate the same kinematical params but now as an experimentalist would 
    // measure them by neglecting the fermi momentum and off-shellness of bound nucleons
    //

    const TLorentzVector & k1 = (neutrino) ? *(neutrino->P4()) : pdummy;  // v 4-p (k1)
    const TLorentzVector & k2 = (fsl)      ? *(fsl->P4())      : pdummy;  // l 4-p (k2)
    const TLorentzVector & p1 = (hitnucl)  ? *(hitnucl->P4())  : pdummy;  // N 4-p (p1)      

    double M  = kNucleonMass; 
    TLorentzVector q  = k1-k2;                     // q=k1-k2, 4-p transfer
    double Q2 = -1 * q.M2();                       // momemtum transfer
    
    double v  = (hitnucl) ? q.Energy()       : -1; // v (E transfer to the nucleus)
    double x, y, W2, W;
    if(!is_coh){ 
    
       x  = (hitnucl) ? 0.5*Q2/(M*v)     : -1; // Bjorken x
       y  = (hitnucl) ? v/k1.Energy()    : -1; // Inelasticity, y = q*P1/k1*P1

       W2 = (hitnucl) ? M*M + 2*M*v - Q2 : -1; // Hadronic Invariant mass ^ 2
       W  = (hitnucl) ? TMath::Sqrt(W2)  : -1; 
    } else{

       v = q.Energy();
       x  =  0.5*Q2/(M*v);      // Bjorken x
       y  = v/k1.Energy();    // Inelasticity, y = q*P1/k1*P1

       W2 = M*M + 2*M*v - Q2;  // Hadronic Invariant mass ^ 2
       W  = TMath::Sqrt(W2); 

    }
  
    double t  = (is_coh || is_dfr) ? kine.t (get_selected) : -1;

    // Get v 4-p at hit nucleon rest-frame
    TLorentzVector k1_rf = k1;         
    if(hitnucl) {
       k1_rf.Boost(-1.*p1.BoostVector());
    }

//    if(is_mec){
//      v = q.Energy();
//      x = 0.5*Q2/(M*v);
//      y = v/k1.Energy();
//      W2 = M*M + 2*M*v - Q2;
//      W = TMath::Sqrt(W2);
//    }

    LOG("gntpc", pDEBUG) 
       << "[Calc] Q2 = " << Q2 << ", W = " << W 
       << ", x = " << x << ", y = " << y << ", t = " << t;

    // Extract more info on the hadronic system
    // Only for QEL/RES/DIS/COH/MEC events
    //
    bool study_hadsyst = (is_qel || is_res || is_dis || is_coh || is_dfr || is_mec || is_singlek);
    
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
      if(pdg::IsPseudoParticle(p->Pdg())) continue;
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
      if(!pdg::IsIon(target->Pdg()) || (is_coh)) {
         vector<int>::const_iterator hiter = final_had_syst.begin();
         for( ; hiter != final_had_syst.end(); ++hiter) {
           prim_had_syst.push_back(*hiter);
         }
      } 
      //to find the true particles emitted from the principal vertex,
      // looping over all Ist=14 particles ok for hA, but doesn't
      // work for hN.  We must now look specifically for these particles.
      int ist_store = -10;
      if(is_res){
	while( (p = (GHepParticle *) piter_prim.Next()) ){
	  ip++;      
	  int ist_comp  = p->Status();
	  if(ist_comp==kIStDecayedState) {
	    ist_store = ip;    //store this mother
	    continue;
	  }
	  //	  LOG("gntpc",pNOTICE) << p->FirstMother()<< "  "<<ist_store;
	  if(p->FirstMother()==ist_store) {
	      prim_had_syst.push_back(ip);
	    }
	}
      }
      if(is_dis){
	while( (p = (GHepParticle *) piter_prim.Next()) ){
	  ip++;      
	  int ist_comp  = p->Status();
	  if(ist_comp==kIStDISPreFragmHadronicState) {
	    ist_store = ip;    //store this mother
	    continue;
	  }
	  if(p->FirstMother()==ist_store) {
	      prim_had_syst.push_back(ip);
	    }
	}
      }
      if(is_qel){
	while( (p = (GHepParticle *) piter_prim.Next()) ){
	  ip++;      
	  int ist_comp  = p->Status();
	  if(ist_comp==kIStNucleonTarget) {
	    ist_store = ip;    //store this mother
	    continue;
	  }
	  //	  LOG("gntpc",pNOTICE) << p->FirstMother()<< "  "<<ist_store;
	  if(p->FirstMother()==ist_store) {
	      prim_had_syst.push_back(ip);
	    }
	}
      }      
      if(is_mec){
	while( (p = (GHepParticle *) piter_prim.Next()) ){
	  ip++;      
	  int ist_comp  = p->Status();
	  if(ist_comp==kIStDecayedState) {
	    ist_store = ip;    //store this mother
	    continue;
	  }
	  //	  LOG("gntpc",pNOTICE) << "MEC: " << p->FirstMother()<< "  "<<ist_store;
	  if(p->FirstMother()==ist_store) {
	      prim_had_syst.push_back(ip);
	    }
	}
      }
      // otherwise loop over all particles and store indices of those which are hadrons
      // created within the nucleus
      /*      else {
	while( (p = (GHepParticle *) piter_prim.Next()) ){
	  ip++;      
	  int ist_comp  = p->Status();
	  if(ist_comp==kIStHadronInTheNucleus) {
	    prim_had_syst.push_back(ip); 
	  }
	  }//particle-loop   */
	//
	// also include gammas from nuclear de-excitations (appearing in the daughter list of the 
	// hit nucleus, earlier than the primary hadronic system extracted above)
	for(int i = target->FirstDaughter(); i <= target->LastDaughter(); i++) {
	  if(i<0) continue;
	  if(event.Particle(i)->Status()==kIStStableFinalState) { prim_had_syst.push_back(i); }
	}      
	//      }//freenuc?
    }//study_hadsystem?
    
    if( count(prim_had_syst.begin(), prim_had_syst.end(), -1) > 0) {
        mcrec->Clear();
 	continue;
    }

    //
    // Al information has been assembled -- Start filling up the tree branches
    //
    brIev        = (int) iev;      
    brNeutrino   = (neutrino) ? neutrino->Pdg() : 0;      
    brFSPrimLept = (fsl) ? fsl->Pdg() : 0;
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
    brIsDfr      = is_dfr;  
    brIsImd      = is_imd;
    brIsSingleK  = is_singlek;    
    brIsNuEL     = is_nuel;  
    brIsEM       = is_em;  
    brIsMec      = is_mec;
    brIsCC       = is_weakcc;  
    brIsNC       = is_weaknc;  
    brIsCharmPro = charm;
    brIsAMNuGamma= is_amnugamma;
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
    brEvRF       = k1_rf.Energy();      
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
    brPl         = k2.P();
    brCosthl     = TMath::Cos( k2.Vect().Angle(k1.Vect()) );

    // Primary hadronic system (from primary neutrino interaction, before FSI)
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
      brResc[j] = p->RescatterCode();     
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

    // Final state (visible) hadronic system
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

    brSumKEf     = (fsl) ? fsl->KinE() : 0;
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
      double hp   = TMath::Sqrt(hpx*hpx + hpy*hpy + hpz*hpz);
      double hm   = p->Mass();     
      double hcth = TMath::Cos( p->P4()->Vect().Angle(k1.Vect()) );

      brPdgf  [j] = hpdg;
      brEf    [j] = hE;
      brPxf   [j] = hpx;
      brPyf   [j] = hpy;
      brPzf   [j] = hpz;
      brPf    [j] = hp;
      brCosthf[j] = hcth;

      brSumKEf += hKE;

      if      ( hpdg == kPdgProton      )  { brNfP++;     brCalResp0 += hKE;        }
      else if ( hpdg == kPdgAntiProton  )  { brNfP++;     brCalResp0 += (hE + 2*hm);}
      else if ( hpdg == kPdgNeutron     )  { brNfN++;     brCalResp0 += hKE;        }
      else if ( hpdg == kPdgAntiNeutron )  { brNfN++;     brCalResp0 += (hE + 2*hm);}
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


  // Copy MC job metadata (gconfig and genv TFolders)
  if(gOptCopyJobMeta) {
    TFolder * genv    = (TFolder*) fin.Get("genv");
    TFolder * gconfig = (TFolder*) fin.Get("gconfig");
    fout.cd();       
    genv    -> Write("genv");
    gconfig -> Write("gconfig");
  }

  fin.Close();

  fout.Write();
  fout.Close();
}
//____________________________________________________________________________________
// GENIE GHEP EVENT TREE FORMAT -> GENIE XML EVENT FILE FORMAT 
//____________________________________________________________________________________
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
      if      (pdg::IsPseudoParticle(p->Pdg())) type = "F";
      else if (pdg::IsParticle      (p->Pdg())) type = "P";
      else if (pdg::IsIon           (p->Pdg())) type = "N";
      
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

      if(p->RescatterCode() != -1) {
        output << "        ";
        output << " <rescatter> " << p->RescatterCode()   << " </rescatter>";
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
//____________________________________________________________________________________
// GENIE GHEP FORMAT -> GHEP MOCK DATA FORMAT
//____________________________________________________________________________________
void ConvertToGHepMock(void)
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
        
  //-- figure out how many events to analyze
  Long64_t nmax = (gOptN<0) ?
       tree->GetEntries() : TMath::Min(tree->GetEntries(), gOptN);
  if (nmax<0) {
    LOG("gntpc", pERROR) << "Number of events = 0";
    return;
  }
  LOG("gntpc", pNOTICE) << "*** Analyzing: " << nmax << " events";

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(kNFGHEP, thdr->runnu);
  ntpw.CustomizeFilename(gOptOutFileName);
  ntpw.Initialize();

  //-- event loop
  for(Long64_t iev = 0; iev < nmax; iev++) {
    tree->GetEntry(iev);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    EventRecord * stripped_event = new EventRecord;
    Interaction * nullint = new Interaction;

    stripped_event -> AttachSummary (nullint);
    stripped_event -> SetWeight     (event.Weight());
    stripped_event -> SetVertex     (*event.Vertex());

    GHepParticle * p = 0;
    TIter iter(&event);
    while( (p = (GHepParticle *)iter.Next()) ) {
       if(!p) continue;
       GHepStatus_t ist = p->Status();
       if(ist!=kIStStableFinalState) continue;
       stripped_event->AddParticle(
          p->Pdg(), ist, -1,-1,-1,-1, *p->P4(), *p->X4());
    }//p

    ntpw.AddEventRecord(iev,stripped_event);

    mcrec->Clear();
  } // event loop

  //-- save the generated MC events
  ntpw.Save();

  //-- rename the output file
      
  fin.Close();
      
  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//____________________________________________________________________________________
// GENIE GHEP EVENT TREE FORMAT -> TRACKER FORMATS
//____________________________________________________________________________________
void ConvertToGTracker(void)
{
  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  gFileMajorVrs = utils::system::GenieMajorVrsNum(thdr->cvstag.GetString().Data());
  gFileMinorVrs = utils::system::GenieMinorVrsNum(thdr->cvstag.GetString().Data());
  gFileRevisVrs = utils::system::GenieRevisVrsNum(thdr->cvstag.GetString().Data());

  //-- get mc record
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
  flux::GJPARCNuFluxPassThroughInfo * flux_info = 0;
  tree->SetBranchAddress("flux", &flux_info);
#else
  LOG("gntpc", pWARN) 
    << "\n Flux drivers are not enabled." 
    << "\n No flux pass-through information will be written-out in the rootracker file"
    << "\n If this isn't what you are supposed to be doing then build GENIE by adding "
    << "--with-flux-drivers in the configuration step.";
#endif

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

       int          ghep_pdgc   = p->Pdg();
       GHepStatus_t ghep_ist    = (GHepStatus_t) p->Status();

       // Neglect all GENIE pseudo-particles
       if(pdg::IsPseudoParticle(ghep_pdgc)) continue;

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
       if (gOptOutFileFormat == kConvFmt_t2k_tracker && pdg::IsIon(p->Pdg())) continue;

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
         if ( pdg::IsIon(p->Pdg()) ) {
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
         const TLorentzVector * p4 = p->P4();
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
      // version 1: 
      //
      // $ info event_num err_flag string_event_code
      // $ info xsec_event diff_xsec_kinematics weight prob
      // $ info vtxx vtxy vtxz vtxt
      // $ info nparticles
      // $ info 0 pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz 
      // $ info 1 pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz 
      // ... ... ...
      // $ info n pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz 
      //
      // version 2:
      //
      // $ info event_num err_flag string_event_code
      // $ info xsec_event diff_xsec_kinematics weight prob
      // $ info vtxx vtxy vtxz vtxt
      // $ info etc
      // $ info nparticles
      // $ info 0 pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz rescatter_code
      // $ info 1 pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz rescatter_code
      // ... ... ...
      // $ info n pdg_code status_code first_daughter last_daughter first_mother last_mother px py pz E x y z t polx poly polz rescatter_code
      //
      // Comments:
      // - The err_flag is a bit field (16 bits)
      // - The string_event_code is a rather long string which encapsulates lot of summary info on the event
      //   (neutrino/nuclear target/hit nucleon/hit quark(if any)/process type/...).
      //   Information on how to parse that string code is available at the T2K event reweighting package.
      // - event_xsec is the event cross section in 1E-38cm^2
      // - diff_event_xsec is the cross section for the selected in 1E-38cm^2/{K^n}
      // - weight is the event weight (1 for unweighted MC)
      // - prob is the event probability (given cross sectios and density-weighted path-length)
      // - vtxx,y,z,t is the vertex position/time in SI units 
      // - etc (added in format vrs >= 2) is used to pass any additional information with event-scope. 
      //   For the time being it is being used to pass the hit quark id (for DIS events) that was lost before 
      //   as SKDETSIM doesn't read the string_event_code where this info is nominally contained.
      //   The quark id is set as (quark_pdg_code) x 10 + i, where i=0 for valence and i=1 for sea quarks. Set to -1 for non-DIS events.
      // - nparticles is the number of particles in the GHEP record (number of $info lines to follow before the start of the JNUBEAM block)
      // - first_/last_daughter first_/last_mother indicate the particle
      // - px,py,pz,E is the particle 4-momentum at the LAB frame (in GeV)
      // - x,y,z,t is the particle 4-position at the hit nucleus coordinate system (in fm, t is not set)
      // - polx,y,z is the particle polarization vector
      // - rescatter_code (added in format vrs >= 2) is a model-dependent intranuclear rescattering code
      //   added to simplify the event analysis (although, in principle, it is recoverable from the particle record).
      //   See $GENIE/src/HadronTransport/INukeHadroFates.h for the meaning of various codes when INTRANUKE is in use.
      //   The rescattering code is stored at the GHEP event record for files generated with GENIE vrs >= 2.5.1.
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

      // insert etc info line for format versions >= 2
      if(gOptVersion >= 2) {
         int quark_id = -1;
         if( interaction->ProcInfo().IsDeepInelastic() && interaction->InitState().Tgt().HitQrkIsSet() ) {
            int quark_pdg = interaction->InitState().Tgt().HitQrkPdg();
            int sorv      = ( interaction->InitState().Tgt().HitSeaQrk() ) ? 1 : 0; // sea q: 1, valence q: 0
            quark_id = 10 * quark_pdg + sorv;
         }
         output << "$ info " << quark_id << endl;
      }

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

        // append rescattering code for format versions >= 2 
        if(gOptVersion >= 2) {
           int rescat_code = -1;
           bool have_rescat_code = false;
           if(gFileMajorVrs >= 2) {
             if(gFileMinorVrs >= 5) {
                if(gFileRevisVrs >= 1) {
                    have_rescat_code = true;
                }
             }
           }
           if(have_rescat_code) {
             rescat_code = p->RescatterCode();
           }
           output << " ";
           output << rescat_code;
        }

        output << endl;
        iparticle++;
      }
      //
      // JNUBEAM flux info - this info will only be available if events were generated 
      // by gT2Kevgen using JNUBEAM flux ntuples as inputs
      //
/*
The T2K/SK collaboration produces MC based on JNUBEAM flux histograms, not flux ntuples.
Therefore JNUBEAM flux pass-through info is never available for generated events.
Commented-out the following info so as not to maintain/support unused code.
If this section is ever re-instated the JNUBEAM passed-through info needs to be matched 
to the latest version of JNUBEAM and an appropriate updated t2k_tracker format needs to 
be agreed with the SKDETSIM maintainers.

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
      PDGLibrary * pdglib = PDGLibrary::Instance();
      if(flux_info) {
         // parent hadron pdg code and decay mode
         output << "$ info " << pdg::GeantToPdg(flux_info->ppid) << " " << flux_info->mode << endl;
         // parent hadron px,py,pz,E at decay
         output << "$ info " << flux_info->ppi * flux_info->npi[0] << " " 
                             << flux_info->ppi * flux_info->npi[1] << " " 
                             << flux_info->ppi * flux_info->npi[2] << " " 
                             << TMath::Sqrt(
                                   TMath::Power(pdglib->Find(pdg::GeantToPdg(flux_info->ppid))->Mass(), 2.)
                                 + TMath::Power(flux_info->ppi, 2.)
                                )  << endl;
         // parent hadron x,y,z,t at decay
         output << "$ info " << flux_info->xpi[0] << " "
                             << flux_info->xpi[1] << " "
                             << flux_info->xpi[2] << " "
                             << "0." 
                             << endl;
         // parent hadron px,py,pz,E at production
         output << "$ info " << flux_info->ppi0 * flux_info->npi0[0] << " "
                             << flux_info->ppi0 * flux_info->npi0[1] << " "
                             << flux_info->ppi0 * flux_info->npi0[2] << " "
                             << TMath::Sqrt(
                                   TMath::Power(pdglib->Find(pdg::GeantToPdg(flux_info->ppid))->Mass(), 2.)
                                 + TMath::Power(flux_info->ppi0, 2.)
                                ) << endl;
         // parent hadron x,y,z,t at production
         output << "$ info " << flux_info->xpi0[0] << " "
                             << flux_info->xpi0[1] << " "
                             << flux_info->xpi0[2] << " "
                             << "0." 
                             << endl;
         // nvtx
         output << "$ info " << output << "$info " << endl;
     }
#endif
*/
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
//____________________________________________________________________________________
// GENIE GHEP EVENT TREE FORMAT -> ROOTRACKER FORMATS 
//____________________________________________________________________________________
void ConvertToGRooTracker(void)
{
  //-- define the output rootracker tree branches

  // event info

  TBits*      brEvtFlags = 0;             // Generator-specific event flags
  TObjString* brEvtCode = 0;              // Generator-specific string with 'event code'
  int         brEvtNum;                   // Event num.
  double      brEvtXSec;                  // Cross section for selected event (1E-38 cm2)
  double      brEvtDXSec;                 // Cross section for selected event kinematics (1E-38 cm2 /{K^n})
  double      brEvtWght;                  // Weight for that event
  double      brEvtProb;                  // Probability for that event (given cross section, path lengths, etc)
  double      brEvtVtx[4];                // Event vertex position in detector coord syst (SI)
  int         brStdHepN;                  // Number of particles in particle array 
  // stdhep-like particle array:
  int         brStdHepPdg   [kNPmax];     // Pdg codes (& generator specific codes for pseudoparticles)
  int         brStdHepStatus[kNPmax];     // Generator-specific status code
  int         brStdHepRescat[kNPmax];     // Hadron transport model - specific rescattering code
  double      brStdHepX4    [kNPmax][4];  // 4-x (x, y, z, t) of particle in hit nucleus frame (fm)
  double      brStdHepP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  double      brStdHepPolz  [kNPmax][3];  // Polarization vector
  int         brStdHepFd    [kNPmax];     // First daughter
  int         brStdHepLd    [kNPmax];     // Last  daughter 
  int         brStdHepFm    [kNPmax];     // First mother
  int         brStdHepLm    [kNPmax];     // Last  mother

  //
  // >> info available at the t2k rootracker variance only
  //
  TObjString* brNuFileName = 0;           // flux file name
  long        brNuFluxEntry;              // entry number from flux file

  // neutrino parent info (passed-through from the beam-line MC / quantities in 'jnubeam' units)
  int         brNuParentPdg;              // parent hadron pdg code
  int         brNuParentDecMode;          // parent hadron decay mode
  double      brNuParentDecP4 [4];        // parent hadron 4-momentum at decay 
  double      brNuParentDecX4 [4];        // parent hadron 4-position at decay
  double      brNuParentProP4 [4];        // parent hadron 4-momentum at production
  double      brNuParentProX4 [4];        // parent hadron 4-position at production
  int         brNuParentProNVtx;          // parent hadron vtx id
  // variables added since 10a flux compatibility changes
  int         brNuIdfd;                   // detector location id
  float       brNuCospibm;                // cosine of the angle between the parent particle direction and the beam direction
  float       brNuCospi0bm;               // same as above except at the production of the parent particle 
  int         brNuGipart;                 // primary particle ID
  float       brNuGpos0[3];               // primary particle starting point
  float       brNuGvec0[3];               // primary particle direction at the starting point
  float       brNuGamom0;                 // momentum of the primary particle at the starting point
  // variables added since 10d and 11a flux compatibility changes
  float       brNuRnu;                    // neutrino r position at ND5/6 plane
  float       brNuXnu[2];                 // neutrino (x,y) position at ND5/6 plane
  // interation history information
  int         brNuNg;                     // number of parents (number of generations) 
  int         brNuGpid[flux::fNgmax];     // particle ID of each ancestor particles
  int         brNuGmec[flux::fNgmax];     // particle production mechanism of each ancestor particle
  float       brNuGcosbm[flux::fNgmax];   // ancestor particle cos(theta) relative to beam
  float       brNuGv[flux::fNgmax][3];    // X,Y and Z vertex position of each ancestor particle
  float       brNuGp[flux::fNgmax][3];    // Px,Px and Pz directional momentum of each ancestor particle
  // out-of-target secondary interactions
  int         brNuGmat[flux::fNgmax];     // material in which the particle originates 
  float       brNuGdistc[flux::fNgmax];   // distance traveled through carbon
  float       brNuGdistal[flux::fNgmax];  // distance traveled through aluminum
  float       brNuGdistti[flux::fNgmax];  // distance traveled through titanium
  float       brNuGdistfe[flux::fNgmax];  // distance traveled through iron
  
  float       brNuNorm;                   // normalisation weight (makes no sense to apply this when generating unweighted events) 
  float       brNuEnusk;                  // "Enu" for SK
  float       brNuNormsk;                 // "norm" for SK
  float       brNuAnorm;                  // Norm component from ND acceptance calculation
  float       brNuVersion;                // Jnubeam version
  int         brNuNtrig;                  // Number of Triggers in simulation
  int         brNuTuneid;                 // Parameter set identifier
  int         brNuPint;                   // Interaction model ID
  float       brNuBpos[2];                // Beam center position
  float       brNuBtilt[2];               // Beam Direction
  float       brNuBrms[2];                // Beam RMS Width
  float       brNuEmit[2];                // Beam Emittance
  float       brNuAlpha[2];               // Beam alpha parameter
  float       brNuHcur[3];                // Horns 1, 2 and 3 Currents 
  int         brNuRand;                   // Random seed
  // codes for T2K cross-generator comparisons 
  int         brNeutCode;                 // NEUT-like reaction code for the GENIE event

  //
  // >> info available at the numi rootracker variance only
  //

  // neutrino parent info (GNuMI passed-through info)
  // see http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/[/v19/output_gnumi.html]
  int        brNumiFluxRun;               // Run number 
  int        brNumiFluxEvtno;             // Event number (proton on target)
  double     brNumiFluxNdxdz;             // Neutrino direction slope (dx/dz) for a random decay
  double     brNumiFluxNdydz;             // Neutrino direction slope (dy/dz) for a random decay
  double     brNumiFluxNpz;               // Neutrino momentum (GeV/c) along z direction (beam axis)
  double     brNumiFluxNenergy;           // Neutrino energy (GeV/c) for a random decay
  double     brNumiFluxNdxdznea;          // Neutrino direction slope (dx/dz) for a decay forced at center of near detector 
  double     brNumiFluxNdydznea;          // Neutrino direction slope (dy/dz) for a decay forced at center of near detector
  double     brNumiFluxNenergyn;          // Neutrino energy for a decay forced at center of near detector 
  double     brNumiFluxNwtnear;           // Neutrino weight for a decay forced at center of near detector 
  double     brNumiFluxNdxdzfar;          // Neutrino direction slope (dx/dz) for a decay forced at center of far detector
  double     brNumiFluxNdydzfar;          // Neutrino direction slope (dy/dz) for a decay forced at center of far detector
  double     brNumiFluxNenergyf;          // Neutrino energy for a decay forced at center of far detector
  double     brNumiFluxNwtfar;            // Neutrino weight for a decay forced at center of far detector
  int        brNumiFluxNorig;             // Obsolete
  int        brNumiFluxNdecay;            // Decay mode that produced neutrino:
                                          // -  1  K0L -> nue pi- e+
                                          // -  2  K0L -> nuebar pi+ e-
                                          // -  3  K0L -> numu pi- mu+
                                          // -  4  K0L -> numubar pi+ mu-
                                          // -  5  K+  -> numu mu+
                                          // -  6  K+  -> nue pi0 e+
                                          // -  7  K+  -> numu pi0 mu+
                                          // -  8  K-  -> numubar mu-
                                          // -  9  K-  -> nuebar pi0 e-
                                          // - 10  K-  -> numubar pi0 mu-
                                          // - 11  mu+ -> numubar nue e+
                                          // - 12  mu- -> numu nuebar e-
                                          // - 13  pi+ -> numu mu+
                                          // - 14  pi- -> numubar mu-
  int        brNumiFluxNtype;             // Neutrino flavor
  double     brNumiFluxVx;                // Position of hadron/muon decay, X coordinate
  double     brNumiFluxVy;                // Position of hadron/muon decay, Y coordinate
  double     brNumiFluxVz;                // Position of hadron/muon decay, Z coordinate
  double     brNumiFluxPdpx;              // Parent momentum at decay point, X - component
  double     brNumiFluxPdpy;              // Parent momentum at decay point, Y - component 
  double     brNumiFluxPdpz;              // Parent momentum at decay point, Z - component 
  double     brNumiFluxPpdxdz;            // Parent dx/dz direction at production 
  double     brNumiFluxPpdydz;            // Parent dy/dz direction at production 
  double     brNumiFluxPppz;              // Parent Z momentum at production 
  double     brNumiFluxPpenergy;          // Parent energy at production 
  int        brNumiFluxPpmedium;          // Tracking medium number where parent was produced 
  int        brNumiFluxPtype;             // Parent particle ID (PDG)
  double     brNumiFluxPpvx;              // Parent production vertex, X coordinate (cm)
  double     brNumiFluxPpvy;              // Parent production vertex, Y coordinate (cm)
  double     brNumiFluxPpvz;              // Parent production vertex, Z coordinate (cm)
  double     brNumiFluxMuparpx;           // Repeat of information above, but for muon neutrino parents 
  double     brNumiFluxMuparpy;           // ...
  double     brNumiFluxMuparpz;           // ...
  double     brNumiFluxMupare;            // ...
  double     brNumiFluxNecm;              // Neutrino energy in COM frame 
  double     brNumiFluxNimpwt;            // Weight of neutrino parent 
  double     brNumiFluxXpoint;            // Unused
  double     brNumiFluxYpoint;            // Unused
  double     brNumiFluxZpoint;            // Unused
  double     brNumiFluxTvx;               // Exit point of parent particle at the target, X coordinate 
  double     brNumiFluxTvy;               // Exit point of parent particle at the target, Y coordinate
  double     brNumiFluxTvz;               // Exit point of parent particle at the target, Z coordinate
  double     brNumiFluxTpx;               // Parent momentum exiting the target, X - component
  double     brNumiFluxTpy;               // Parent momentum exiting the target, Y - component
  double     brNumiFluxTpz;               // Parent momentum exiting the target, Z - component
  double     brNumiFluxTptype;            // Parent particle ID exiting the target
  double     brNumiFluxTgen;              // Parent generation in cascade
                                          // -  1  primary proton 
                                          // -  2  particles produced by proton interaction
                                          // -  3  particles produced by interactions of the 2's, ... 
  double     brNumiFluxTgptype;           // Type of particle that created a particle flying of the target 
  double     brNumiFluxTgppx;             // Momentum of a particle, that created a particle that flies off 
                                          //  the target (at the interaction point), X - component
  double     brNumiFluxTgppy;             // Momentum of a particle, that created a particle that flies off 
                                          //  the target (at the interaction point), Y - component
  double     brNumiFluxTgppz;             // Momentum of a particle, that created a particle that flies off 
                                          //  the target (at the interaction point), Z - component
  double     brNumiFluxTprivx;            // Primary particle interaction vertex, X coordinate
  double     brNumiFluxTprivy;            // Primary particle interaction vertex, Y coordinate
  double     brNumiFluxTprivz;            // Primary particle interaction vertex, Z coordinate
  double     brNumiFluxBeamx;             // Primary proton origin, X coordinate
  double     brNumiFluxBeamy;             // Primary proton origin, Y coordinate
  double     brNumiFluxBeamz;             // Primary proton origin, Z coordinate
  double     brNumiFluxBeampx;            // Primary proton momentum, X - component
  double     brNumiFluxBeampy;            // Primary proton momentum, Y - component
  double     brNumiFluxBeampz;            // Primary proton momentum, Z - component

  //-- open the output ROOT file
  TFile fout(gOptOutFileName.c_str(), "RECREATE");

  //-- create the output ROOT tree
  TTree * rootracker_tree = new TTree("gRooTracker","GENIE event tree rootracker format");

  //-- is it a `mock data' variance?
  bool hide_truth = (gOptOutFileFormat == kConvFmt_rootracker_mock_data);

  //-- create the output ROOT tree branches

  // branches common to all rootracker(_mock_data) formats
  if(!hide_truth) {
    // full version
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
    rootracker_tree->Branch("StdHepRescat",     brStdHepRescat,    "StdHepRescat[StdHepN]/I"); 
    rootracker_tree->Branch("StdHepX4",         brStdHepX4,        "StdHepX4[StdHepN][4]/D"); 
    rootracker_tree->Branch("StdHepP4",         brStdHepP4,        "StdHepP4[StdHepN][4]/D"); 
    rootracker_tree->Branch("StdHepPolz",       brStdHepPolz,      "StdHepPolz[StdHepN][3]/D"); 
    rootracker_tree->Branch("StdHepFd",         brStdHepFd,        "StdHepFd[StdHepN]/I"); 
    rootracker_tree->Branch("StdHepLd",         brStdHepLd,        "StdHepLd[StdHepN]/I"); 
    rootracker_tree->Branch("StdHepFm",         brStdHepFm,        "StdHepFm[StdHepN]/I"); 
    rootracker_tree->Branch("StdHepLm",         brStdHepLm,        "StdHepLm[StdHepN]/I"); 
  } else {
    // for mock_data variances
    rootracker_tree->Branch("EvtNum",          &brEvtNum,          "EvtNum/I");             
    rootracker_tree->Branch("EvtWght",         &brEvtWght,         "EvtWght/D");            
    rootracker_tree->Branch("EvtVtx",           brEvtVtx,          "EvtVtx[4]/D");             
    rootracker_tree->Branch("StdHepN",         &brStdHepN,         "StdHepN/I");              
    rootracker_tree->Branch("StdHepPdg",        brStdHepPdg,       "StdHepPdg[StdHepN]/I");  
    rootracker_tree->Branch("StdHepX4",         brStdHepX4,        "StdHepX4[StdHepN][4]/D"); 
    rootracker_tree->Branch("StdHepP4",         brStdHepP4,        "StdHepP4[StdHepN][4]/D"); 
  }

  // extra branches of the t2k rootracker variance
  if(gOptOutFileFormat == kConvFmt_t2k_rootracker) 
  {
    // NEUT-like reaction code
    rootracker_tree->Branch("G2NeutEvtCode",   &brNeutCode,        "G2NeutEvtCode/I");   
    // JNUBEAM pass-through info
    rootracker_tree->Branch("NuFileName", "TObjString", &brNuFileName, 32000, 1); 
    rootracker_tree->Branch("NuParentPdg",     &brNuParentPdg,     "NuParentPdg/I");       
    rootracker_tree->Branch("NuParentDecMode", &brNuParentDecMode, "NuParentDecMode/I");   
    rootracker_tree->Branch("NuParentDecP4",    brNuParentDecP4,   "NuParentDecP4[4]/D");     
    rootracker_tree->Branch("NuParentDecX4",    brNuParentDecX4,   "NuParentDecX4[4]/D");     
    rootracker_tree->Branch("NuParentProP4",    brNuParentProP4,   "NuParentProP4[4]/D");     
    rootracker_tree->Branch("NuParentProX4",    brNuParentProX4,   "NuParentProX4[4]/D");     
    rootracker_tree->Branch("NuParentProNVtx", &brNuParentProNVtx, "NuParentProNVtx/I");   
    // Branches added since JNUBEAM '10a' compatibility changes
    rootracker_tree->Branch("NuFluxEntry",     &brNuFluxEntry,     "NuFluxEntry/L");
    rootracker_tree->Branch("NuIdfd",          &brNuIdfd,          "NuIdfd/I");
    rootracker_tree->Branch("NuCospibm",       &brNuCospibm,       "NuCospibm/F");
    rootracker_tree->Branch("NuCospi0bm",      &brNuCospi0bm,      "NuCospi0bm/F");
    rootracker_tree->Branch("NuGipart",        &brNuGipart,        "NuGipart/I");
    rootracker_tree->Branch("NuGpos0",          brNuGpos0,         "NuGpos0[3]/F");
    rootracker_tree->Branch("NuGvec0",          brNuGvec0,         "NuGvec0[3]/F");
    rootracker_tree->Branch("NuGamom0",        &brNuGamom0,        "NuGamom0/F");
    // Branches added since JNUBEAM '10d' compatibility changes
    rootracker_tree->Branch("NuXnu",        brNuXnu, "NuXnu[2]/F");
    rootracker_tree->Branch("NuRnu",       &brNuRnu,      "NuRnu/F");
    rootracker_tree->Branch("NuNg",        &brNuNg,       "NuNg/I");
    rootracker_tree->Branch("NuGpid",       brNuGpid,     "NuGpid[NuNg]/I");
    rootracker_tree->Branch("NuGmec",       brNuGmec,     "NuGmec[NuNg]/I");
    rootracker_tree->Branch("NuGv",         brNuGv,       "NuGv[NuNg][3]/F");
    rootracker_tree->Branch("NuGp",         brNuGp,       "NuGp[NuNg][3]/F");
    rootracker_tree->Branch("NuGcosbm",     brNuGcosbm,   "NuGcosbm[NuNg]/F");
    rootracker_tree->Branch("NuGmat",       brNuGmat,     "NuGmat[NuNg]/I");
    rootracker_tree->Branch("NuGdistc",     brNuGdistc,   "NuGdistc[NuNg]/F");
    rootracker_tree->Branch("NuGdistal",    brNuGdistal,  "NuGdistal[NuNg]/F");
    rootracker_tree->Branch("NuGdistti",    brNuGdistti,  "NuGdistti[NuNg]/F");
    rootracker_tree->Branch("NuGdistfe",    brNuGdistfe,  "NuGdistfe[NuNg]/F");
    rootracker_tree->Branch("NuNorm",      &brNuNorm,     "NuNorm/F");
    rootracker_tree->Branch("NuEnusk",     &brNuEnusk,    "NuEnusk/F");
    rootracker_tree->Branch("NuNormsk",    &brNuNormsk,   "NuNormsk/F");
    rootracker_tree->Branch("NuAnorm",     &brNuAnorm,    "NuAnorm/F");
    rootracker_tree->Branch("NuVersion",   &brNuVersion,  "NuVersion/F");
    rootracker_tree->Branch("NuNtrig",     &brNuNtrig,    "NuNtrig/I");
    rootracker_tree->Branch("NuTuneid",    &brNuTuneid,   "NuTuneid/I");
    rootracker_tree->Branch("NuPint",      &brNuPint,     "NuPint/I");
    rootracker_tree->Branch("NuBpos",       brNuBpos,     "NuBpos[2]/F");
    rootracker_tree->Branch("NuBtilt",      brNuBtilt,    "NuBtilt[2]/F");
    rootracker_tree->Branch("NuBrms",       brNuBrms,     "NuBrms[2]/F");
    rootracker_tree->Branch("NuEmit",       brNuEmit,     "NuEmit[2]/F");
    rootracker_tree->Branch("NuAlpha",      brNuAlpha,    "NuAlpha[2]/F");
    rootracker_tree->Branch("NuHcur",       brNuHcur,     "NuHcur[3]/F");
    rootracker_tree->Branch("NuRand",      &brNuRand,     "NuRand/I");

  }

  // extra branches of the numi rootracker variance
  if(gOptOutFileFormat == kConvFmt_numi_rootracker) 
  {
   // GNuMI pass-through info
   rootracker_tree->Branch("NumiFluxRun",      &brNumiFluxRun,       "NumiFluxRun/I");
   rootracker_tree->Branch("NumiFluxEvtno",    &brNumiFluxEvtno,     "NumiFluxEvtno/I");
   rootracker_tree->Branch("NumiFluxNdxdz",    &brNumiFluxNdxdz,     "NumiFluxNdxdz/D");
   rootracker_tree->Branch("NumiFluxNdydz",    &brNumiFluxNdydz,     "NumiFluxNdydz/D");
   rootracker_tree->Branch("NumiFluxNpz",      &brNumiFluxNpz,       "NumiFluxNpz/D");
   rootracker_tree->Branch("NumiFluxNenergy",  &brNumiFluxNenergy,   "NumiFluxNenergy/D");
   rootracker_tree->Branch("NumiFluxNdxdznea", &brNumiFluxNdxdznea,  "NumiFluxNdxdznea/D");
   rootracker_tree->Branch("NumiFluxNdydznea", &brNumiFluxNdydznea,  "NumiFluxNdydznea/D");
   rootracker_tree->Branch("NumiFluxNenergyn", &brNumiFluxNenergyn,  "NumiFluxNenergyn/D");
   rootracker_tree->Branch("NumiFluxNwtnear",  &brNumiFluxNwtnear,   "NumiFluxNwtnear/D");
   rootracker_tree->Branch("NumiFluxNdxdzfar", &brNumiFluxNdxdzfar,  "NumiFluxNdxdzfar/D");
   rootracker_tree->Branch("NumiFluxNdydzfar", &brNumiFluxNdydzfar,  "NumiFluxNdydzfar/D");
   rootracker_tree->Branch("NumiFluxNenergyf", &brNumiFluxNenergyf,  "NumiFluxNenergyf/D");
   rootracker_tree->Branch("NumiFluxNwtfar",   &brNumiFluxNwtfar,    "NumiFluxNwtfar/D");
   rootracker_tree->Branch("NumiFluxNorig",    &brNumiFluxNorig,     "NumiFluxNorig/I");
   rootracker_tree->Branch("NumiFluxNdecay",   &brNumiFluxNdecay,    "NumiFluxNdecay/I");
   rootracker_tree->Branch("NumiFluxNtype",    &brNumiFluxNtype,     "NumiFluxNtype/I");
   rootracker_tree->Branch("NumiFluxVx",       &brNumiFluxVx,        "NumiFluxVx/D");
   rootracker_tree->Branch("NumiFluxVy",       &brNumiFluxVy,        "NumiFluxVy/D");
   rootracker_tree->Branch("NumiFluxVz",       &brNumiFluxVz,        "NumiFluxVz/D");
   rootracker_tree->Branch("NumiFluxPdpx",     &brNumiFluxPdpx,      "NumiFluxPdpx/D");
   rootracker_tree->Branch("NumiFluxPdpy",     &brNumiFluxPdpy,      "NumiFluxPdpy/D");
   rootracker_tree->Branch("NumiFluxPdpz",     &brNumiFluxPdpz,      "NumiFluxPdpz/D");
   rootracker_tree->Branch("NumiFluxPpdxdz",   &brNumiFluxPpdxdz,    "NumiFluxPpdxdz/D");
   rootracker_tree->Branch("NumiFluxPpdydz",   &brNumiFluxPpdydz,    "NumiFluxPpdydz/D");
   rootracker_tree->Branch("NumiFluxPppz",     &brNumiFluxPppz,      "NumiFluxPppz/D");
   rootracker_tree->Branch("NumiFluxPpenergy", &brNumiFluxPpenergy,  "NumiFluxPpenergy/D");
   rootracker_tree->Branch("NumiFluxPpmedium", &brNumiFluxPpmedium,  "NumiFluxPpmedium/I");
   rootracker_tree->Branch("NumiFluxPtype",    &brNumiFluxPtype,     "NumiFluxPtype/I");
   rootracker_tree->Branch("NumiFluxPpvx",     &brNumiFluxPpvx,      "NumiFluxPpvx/D");
   rootracker_tree->Branch("NumiFluxPpvy",     &brNumiFluxPpvy,      "NumiFluxPpvy/D");
   rootracker_tree->Branch("NumiFluxPpvz",     &brNumiFluxPpvz,      "NumiFluxPpvz/D");
   rootracker_tree->Branch("NumiFluxMuparpx",  &brNumiFluxMuparpx,   "NumiFluxMuparpx/D");
   rootracker_tree->Branch("NumiFluxMuparpy",  &brNumiFluxMuparpy,   "NumiFluxMuparpy/D");
   rootracker_tree->Branch("NumiFluxMuparpz",  &brNumiFluxMuparpz,   "NumiFluxMuparpz/D");
   rootracker_tree->Branch("NumiFluxMupare",   &brNumiFluxMupare,    "NumiFluxMupare/D");
   rootracker_tree->Branch("NumiFluxNecm",     &brNumiFluxNecm,      "NumiFluxNecm/D");
   rootracker_tree->Branch("NumiFluxNimpwt",   &brNumiFluxNimpwt,    "NumiFluxNimpwt/D");
   rootracker_tree->Branch("NumiFluxXpoint",   &brNumiFluxXpoint,    "NumiFluxXpoint/D");
   rootracker_tree->Branch("NumiFluxYpoint",   &brNumiFluxYpoint,    "NumiFluxYpoint/D");
   rootracker_tree->Branch("NumiFluxZpoint",   &brNumiFluxZpoint,    "NumiFluxZpoint/D");
   rootracker_tree->Branch("NumiFluxTvx",      &brNumiFluxTvx,       "NumiFluxTvx/D");
   rootracker_tree->Branch("NumiFluxTvy",      &brNumiFluxTvy,       "NumiFluxTvy/D");
   rootracker_tree->Branch("NumiFluxTvz",      &brNumiFluxTvz,       "NumiFluxTvz/D");
   rootracker_tree->Branch("NumiFluxTpx",      &brNumiFluxTpx,       "NumiFluxTpx/D");
   rootracker_tree->Branch("NumiFluxTpy",      &brNumiFluxTpy,       "NumiFluxTpy/D");
   rootracker_tree->Branch("NumiFluxTpz",      &brNumiFluxTpz,       "NumiFluxTpz/D");
   rootracker_tree->Branch("NumiFluxTptype",   &brNumiFluxTptype,    "NumiFluxTptype/I");
   rootracker_tree->Branch("NumiFluxTgen",     &brNumiFluxTgen,      "NumiFluxTgen/I");
   rootracker_tree->Branch("NumiFluxTgptype",  &brNumiFluxTgptype,   "NumiFluxTgptype/I");
   rootracker_tree->Branch("NumiFluxTgppx",    &brNumiFluxTgppx,     "NumiFluxTgppx/D");
   rootracker_tree->Branch("NumiFluxTgppy",    &brNumiFluxTgppy,     "NumiFluxTgppy/D");
   rootracker_tree->Branch("NumiFluxTgppz",    &brNumiFluxTgppz,     "NumiFluxTgppz/D");
   rootracker_tree->Branch("NumiFluxTprivx",   &brNumiFluxTprivx,    "NumiFluxTprivx/D");
   rootracker_tree->Branch("NumiFluxTprivy",   &brNumiFluxTprivy,    "NumiFluxTprivy/D");
   rootracker_tree->Branch("NumiFluxTprivz",   &brNumiFluxTprivz,    "NumiFluxTprivz/D");
   rootracker_tree->Branch("NumiFluxBeamx",    &brNumiFluxBeamx,     "NumiFluxBeamx/D");
   rootracker_tree->Branch("NumiFluxBeamy",    &brNumiFluxBeamy,     "NumiFluxBeamy/D");
   rootracker_tree->Branch("NumiFluxBeamz",    &brNumiFluxBeamz,     "NumiFluxBeamz/D");
   rootracker_tree->Branch("NumiFluxBeampx",   &brNumiFluxBeampx,    "NumiFluxBeampx/D");
   rootracker_tree->Branch("NumiFluxBeampy",   &brNumiFluxBeampy,    "NumiFluxBeampy/D");
   rootracker_tree->Branch("NumiFluxBeampz",   &brNumiFluxBeampz,    "NumiFluxBeampz/D");
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

  //-- print-out metadata associated with the input event file in case the
  //   event file was generated using the gT2Kevgen driver
  //   (assuming this is the case if the requested output format is the t2k_rootracker format)
  if(gOptOutFileFormat == kConvFmt_t2k_rootracker) 
  {
    // Check can find the MetaData
    genie::utils::T2KEvGenMetaData * metadata = NULL;
    metadata = (genie::utils::T2KEvGenMetaData *) gtree->GetUserInfo()->At(0);
    if(metadata){
      LOG("gntpc", pINFO) << "Found T2KMetaData!";
      LOG("gntpc", pINFO) << *metadata;
    }
    else { 
      LOG("gntpc", pWARN) 
        << "Could not find T2KMetaData attached to the event tree!";
    }
  }

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
  flux::GJPARCNuFluxPassThroughInfo * jnubeam_flux_info = 0;
  if(gOptOutFileFormat == kConvFmt_t2k_rootracker) {
     gtree->SetBranchAddress("flux", &jnubeam_flux_info);
  }
  flux::GNuMIFluxPassThroughInfo * gnumi_flux_info = 0;
  if(gOptOutFileFormat == kConvFmt_numi_rootracker) {
     gtree->SetBranchAddress("flux", &gnumi_flux_info);
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
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
    if(gOptOutFileFormat == kConvFmt_t2k_rootracker) {
       if(jnubeam_flux_info) {
          LOG("gntpc", pINFO) << *jnubeam_flux_info;
       } else {
          LOG("gntpc", pINFO) << "No JNUBEAM flux info associated with this event";
       }
    }
#endif

    //
    // clear output tree branches
    //
    if(brEvtFlags) delete brEvtFlags;
    brEvtFlags  = 0;
    if(brEvtCode) delete brEvtCode;
    brEvtCode   = 0;
    brEvtNum    = 0;    
    brEvtXSec   = 0;
    brEvtDXSec  = 0;
    brEvtWght   = 0;
    brEvtProb   = 0;
    for(int k=0; k<4; k++) { 
      brEvtVtx[k] = 0;
    }
    brStdHepN = 0; 
    for(int i=0; i<kNPmax; i++) {
       brStdHepPdg   [i] =  0;  
       brStdHepStatus[i] = -1;  
       brStdHepRescat[i] = -1;  
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
    brNuFluxEntry = -1;
    brNuIdfd = -999999;
    brNuCospibm = -999999.;
    brNuCospi0bm = -999999.;
    brNuGipart = -1;
    brNuGamom0 = -999999.;   
    for(int k=0; k< 3; k++){
      brNuGvec0[k] = -999999.;
      brNuGpos0[k] = -999999.;
    }    
    // variables added since 10d flux compatibility changes
    for(int k=0; k<2; k++) {
      brNuXnu[k] = brNuBpos[k] = brNuBtilt[k] = brNuBrms[k] = brNuEmit[k] = brNuAlpha[k] = -999999.; 
    }
    for(int k=0; k<3; k++) brNuHcur[k] = -999999.; 
    for(int np = 0; np < flux::fNgmax; np++){
        for(int  k=0; k<3; k++){
          brNuGv[np][k] = -999999.;
          brNuGp[np][k] = -999999.;
        }
      brNuGpid[np] = -999999;
      brNuGmec[np] = -999999;
      brNuGmat[np] = -999999;
      brNuGcosbm[np]  = -999999.;
      brNuGdistc[np]  = -999999.;
      brNuGdistal[np] = -999999.;
      brNuGdistti[np] = -999999.;
      brNuGdistfe[np] = -999999.;
    }  
    brNuNg     = -999999;
    brNuRnu    = -999999.;
    brNuNorm   = -999999.;
    brNuEnusk  = -999999.;
    brNuNormsk = -999999.;
    brNuAnorm  = -999999.;
    brNuVersion= -999999.;
    brNuNtrig  = -999999;
    brNuTuneid = -999999;
    brNuPint   = -999999;
    brNuRand = -999999;
    if(brNuFileName) delete brNuFileName;
    brNuFileName = 0;

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

    int iparticle=0;
    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
        assert(p);

        // for mock_data variances write out only stable final state particles
        if(hide_truth && p->Status() != kIStStableFinalState) continue;

        brStdHepPdg   [iparticle] = p->Pdg(); 
        brStdHepStatus[iparticle] = (int) p->Status(); 
        brStdHepRescat[iparticle] = p->RescatterCode(); 
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
    brStdHepN = iparticle; 

    //
    // fill in additional info for the t2k_rootracker format
    //
    if(gOptOutFileFormat == kConvFmt_t2k_rootracker) {

      // map GENIE event to NEUT reaction codes
      brNeutCode = utils::ghep::NeutReactionCode(&event);

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
      // Copy flux info if this is the t2k rootracker variance.
      // The flux may not be available, eg if events were generated using plain flux 
      // histograms and not the JNUBEAM simulation's output flux ntuples.
      PDGLibrary * pdglib = PDGLibrary::Instance();
      if(jnubeam_flux_info) {
        brNuParentPdg       = pdg::GeantToPdg(jnubeam_flux_info->ppid);
        brNuParentDecMode   = jnubeam_flux_info->mode;

        brNuParentDecP4 [0] = jnubeam_flux_info->ppi * jnubeam_flux_info->npi[0]; // px
        brNuParentDecP4 [1] = jnubeam_flux_info->ppi * jnubeam_flux_info->npi[1]; // py
        brNuParentDecP4 [2] = jnubeam_flux_info->ppi * jnubeam_flux_info->npi[2]; // px
        brNuParentDecP4 [3] = TMath::Sqrt(
                                 TMath::Power(pdglib->Find(brNuParentPdg)->Mass(), 2.)
                               + TMath::Power(jnubeam_flux_info->ppi, 2.)
                              ); // E
        brNuParentDecX4 [0] = jnubeam_flux_info->xpi[0]; // x
        brNuParentDecX4 [1] = jnubeam_flux_info->xpi[1]; // y       
        brNuParentDecX4 [2] = jnubeam_flux_info->xpi[2]; // x   
        brNuParentDecX4 [3] = 0;                 // t

        brNuParentProP4 [0] = jnubeam_flux_info->ppi0 * jnubeam_flux_info->npi0[0]; // px
        brNuParentProP4 [1] = jnubeam_flux_info->ppi0 * jnubeam_flux_info->npi0[1]; // py
        brNuParentProP4 [2] = jnubeam_flux_info->ppi0 * jnubeam_flux_info->npi0[2]; // px
        brNuParentProP4 [3] = TMath::Sqrt(
                                TMath::Power(pdglib->Find(brNuParentPdg)->Mass(), 2.)
                              + TMath::Power(jnubeam_flux_info->ppi0, 2.)
                              ); // E
        brNuParentProX4 [0] = jnubeam_flux_info->xpi0[0]; // x
        brNuParentProX4 [1] = jnubeam_flux_info->xpi0[1]; // y       
        brNuParentProX4 [2] = jnubeam_flux_info->xpi0[2]; // x   
        brNuParentProX4 [3] = 0;                // t

        brNuParentProNVtx   = jnubeam_flux_info->nvtx0;

        // Copy info added post JNUBEAM '10a' compatibility changes 
        brNuFluxEntry = jnubeam_flux_info->fluxentry;
        brNuIdfd = jnubeam_flux_info->idfd;
        brNuCospibm = jnubeam_flux_info->cospibm;
        brNuCospi0bm = jnubeam_flux_info->cospi0bm;
        brNuGipart = jnubeam_flux_info->gipart;
        brNuGamom0 = jnubeam_flux_info->gamom0;
        for(int k=0; k<3; k++){
            brNuGpos0[k] = (double) jnubeam_flux_info->gpos0[k];
            brNuGvec0[k] = (double) jnubeam_flux_info->gvec0[k];
        }
        // Copy info added post JNUBEAM '10d' compatibility changes 
        brNuXnu[0] = (double) jnubeam_flux_info->xnu;
        brNuXnu[1] = (double) jnubeam_flux_info->ynu;
        brNuRnu    = (double) jnubeam_flux_info->rnu; 
        for(int k=0; k<2; k++){
          brNuBpos[k] = (double) jnubeam_flux_info->bpos[k];
          brNuBtilt[k] = (double) jnubeam_flux_info->btilt[k];
          brNuBrms[k] = (double) jnubeam_flux_info->brms[k];
          brNuEmit[k] = (double) jnubeam_flux_info->emit[k];
          brNuAlpha[k] = (double) jnubeam_flux_info->alpha[k];
        } 
        for(int k=0; k<3; k++) brNuHcur[k] = jnubeam_flux_info->hcur[k]; 
        for(int np = 0; np < flux::fNgmax; np++){
          brNuGv[np][0] = jnubeam_flux_info->gvx[np];
          brNuGv[np][1] = jnubeam_flux_info->gvy[np];
          brNuGv[np][2] = jnubeam_flux_info->gvz[np];
          brNuGp[np][0] = jnubeam_flux_info->gpx[np];
          brNuGp[np][1] = jnubeam_flux_info->gpy[np];
          brNuGp[np][2] = jnubeam_flux_info->gpz[np];
          brNuGpid[np]  = jnubeam_flux_info->gpid[np];
          brNuGmec[np]  = jnubeam_flux_info->gmec[np];
          brNuGcosbm[np]  = jnubeam_flux_info->gcosbm[np];
          brNuGmat[np]    = jnubeam_flux_info->gmat[np];
          brNuGdistc[np]  = jnubeam_flux_info->gdistc[np];
          brNuGdistal[np] = jnubeam_flux_info->gdistal[np];
          brNuGdistti[np] = jnubeam_flux_info->gdistti[np];
          brNuGdistfe[np] = jnubeam_flux_info->gdistfe[np];
        }  
        brNuNg     = jnubeam_flux_info->ng;
        brNuNorm   = jnubeam_flux_info->norm;
        brNuEnusk  = jnubeam_flux_info->Enusk;
        brNuNormsk = jnubeam_flux_info->normsk;
        brNuAnorm  = jnubeam_flux_info->anorm;
        brNuVersion= jnubeam_flux_info->version;
        brNuNtrig  = jnubeam_flux_info->ntrig;
        brNuTuneid = jnubeam_flux_info->tuneid;
        brNuPint   = jnubeam_flux_info->pint;
        brNuRand   = jnubeam_flux_info->rand;
        brNuFileName = new TObjString(jnubeam_flux_info->fluxfilename.c_str()); 
      }//jnubeam_flux_info
#endif
    }//kConvFmt_t2k_rootracker

    //
    // fill in additional info for the numi_rootracker format
    //
    if(gOptOutFileFormat == kConvFmt_numi_rootracker) {
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
     // Copy flux info if this is the numi rootracker variance.
     if(gnumi_flux_info) {
       brNumiFluxRun      = gnumi_flux_info->run;
       brNumiFluxEvtno    = gnumi_flux_info->evtno;
       brNumiFluxNdxdz    = gnumi_flux_info->ndxdz;
       brNumiFluxNdydz    = gnumi_flux_info->ndydz;
       brNumiFluxNpz      = gnumi_flux_info->npz;
       brNumiFluxNenergy  = gnumi_flux_info->nenergy;
       brNumiFluxNdxdznea = gnumi_flux_info->ndxdznea;
       brNumiFluxNdydznea = gnumi_flux_info->ndydznea;
       brNumiFluxNenergyn = gnumi_flux_info->nenergyn;
       brNumiFluxNwtnear  = gnumi_flux_info->nwtnear;
       brNumiFluxNdxdzfar = gnumi_flux_info->ndxdzfar;
       brNumiFluxNdydzfar = gnumi_flux_info->ndydzfar;
       brNumiFluxNenergyf = gnumi_flux_info->nenergyf;
       brNumiFluxNwtfar   = gnumi_flux_info->nwtfar;
       brNumiFluxNorig    = gnumi_flux_info->norig;
       brNumiFluxNdecay   = gnumi_flux_info->ndecay;
       brNumiFluxNtype    = gnumi_flux_info->ntype;
       brNumiFluxVx       = gnumi_flux_info->vx;
       brNumiFluxVy       = gnumi_flux_info->vy;
       brNumiFluxVz       = gnumi_flux_info->vz;
       brNumiFluxPdpx     = gnumi_flux_info->pdpx;
       brNumiFluxPdpy     = gnumi_flux_info->pdpy;
       brNumiFluxPdpz     = gnumi_flux_info->pdpz;
       brNumiFluxPpdxdz   = gnumi_flux_info->ppdxdz;
       brNumiFluxPpdydz   = gnumi_flux_info->ppdydz;
       brNumiFluxPppz     = gnumi_flux_info->pppz;
       brNumiFluxPpenergy = gnumi_flux_info->ppenergy;
       brNumiFluxPpmedium = gnumi_flux_info->ppmedium;
       brNumiFluxPtype    = gnumi_flux_info->ptype;
       brNumiFluxPpvx     = gnumi_flux_info->ppvx;
       brNumiFluxPpvy     = gnumi_flux_info->ppvy;
       brNumiFluxPpvz     = gnumi_flux_info->ppvz;
       brNumiFluxMuparpx  = gnumi_flux_info->muparpx;
       brNumiFluxMuparpy  = gnumi_flux_info->muparpy;
       brNumiFluxMuparpz  = gnumi_flux_info->muparpz;
       brNumiFluxMupare   = gnumi_flux_info->mupare;
       brNumiFluxNecm     = gnumi_flux_info->necm;
       brNumiFluxNimpwt   = gnumi_flux_info->nimpwt;
       brNumiFluxXpoint   = gnumi_flux_info->xpoint;
       brNumiFluxYpoint   = gnumi_flux_info->ypoint;
       brNumiFluxZpoint   = gnumi_flux_info->zpoint;
       brNumiFluxTvx      = gnumi_flux_info->tvx;
       brNumiFluxTvy      = gnumi_flux_info->tvy;
       brNumiFluxTvz      = gnumi_flux_info->tvz;
       brNumiFluxTpx      = gnumi_flux_info->tpx;
       brNumiFluxTpy      = gnumi_flux_info->tpy;
       brNumiFluxTpz      = gnumi_flux_info->tpz;
       brNumiFluxTptype   = gnumi_flux_info->tptype;
       brNumiFluxTgen     = gnumi_flux_info->tgen;
       brNumiFluxTgptype  = gnumi_flux_info->tgptype;
       brNumiFluxTgppx    = gnumi_flux_info->tgppx;
       brNumiFluxTgppy    = gnumi_flux_info->tgppy;
       brNumiFluxTgppz    = gnumi_flux_info->tgppz;
       brNumiFluxTprivx   = gnumi_flux_info->tprivx;
       brNumiFluxTprivy   = gnumi_flux_info->tprivy;
       brNumiFluxTprivz   = gnumi_flux_info->tprivz;
       brNumiFluxBeamx    = gnumi_flux_info->beamx;
       brNumiFluxBeamy    = gnumi_flux_info->beamy;
       brNumiFluxBeamz    = gnumi_flux_info->beamz;
       brNumiFluxBeampx   = gnumi_flux_info->beampx;
       brNumiFluxBeampy   = gnumi_flux_info->beampy;
       brNumiFluxBeampz   = gnumi_flux_info->beampz;
     } // gnumi_flux_info
#endif
    } // kConvFmt_numi_rootracker

    // fill tree
    rootracker_tree->Fill();
    mcrec->Clear();

  } // event loop

  // Copy POT normalization for the generated sample
  double pot = gtree->GetWeight();
  rootracker_tree->SetWeight(pot);

  // Copy MC job metadata (gconfig and genv TFolders)
  if(gOptCopyJobMeta) {
    TFolder * genv    = (TFolder*) fin.Get("genv");
    TFolder * gconfig = (TFolder*) fin.Get("gconfig");    
    fout.cd();
    genv    -> Write("genv");
    gconfig -> Write("gconfig");
  }

  fin.Close();

  fout.Write();
  fout.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//____________________________________________________________________________________
// GENIE GHEP EVENT TREE -> NEUGEN-style format for AGKY studies 
//____________________________________________________________________________________
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
      GHepParticle * particle = event.Particle(id);
      int pdg = particle->Pdg();
      double px = particle->P4()->Px();
      double py = particle->P4()->Py();
      double pz = particle->P4()->Pz();
      double E  = particle->P4()->Energy();
      double m  = particle->P4()->M();
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
//____________________________________________________________________________________
// GENIE GHEP EVENT TREE -> Summary tree for INTRANUKE studies 
//____________________________________________________________________________________
void ConvertToGINuke(void)
{
  //-- output tree branch variables
  //
  int    brIEv        = 0;  // Event number
  int    brProbe      = 0;  // Incident hadron code
  int    brTarget     = 0;  // Nuclear target pdg code (10LZZZAAAI)
  double brKE         = 0;  // Probe kinetic energy
  double brE          = 0;  // Probe energy
  double brP          = 0;  // Probe momentum
  int    brTgtA       = 0;  // Target A (mass   number)
  int    brTgtZ       = 0;  // Target Z (atomic number)
  double brVtxX       = 0;  // "Vertex x" (initial placement of h /in h+A events/ on the nuclear boundary)
  double brVtxY       = 0;  // "Vertex y"
  double brVtxZ       = 0;  // "Vertex z"
  int    brProbeFSI   = 0;  // Rescattering code for incident hadron
  double brDist       = 0;  // Distance travelled by h before interacting (if at all before escaping)
  int    brNh         = 0;  // Number of final state hadrons
  int    brPdgh  [kNPmax];  // Pdg code of i^th final state hadron
  double brEh    [kNPmax];  // Energy   of i^th final state hadron
  double brPh    [kNPmax];  // P        of i^th final state hadron
  double brPxh   [kNPmax];  // Px       of i^th final state hadron
  double brPyh   [kNPmax];  // Py       of i^th final state hadron
  double brPzh   [kNPmax];  // Pz       of i^th final state hadron
  double brCosth [kNPmax];  // Cos(th)  of i^th final state hadron
  double brMh    [kNPmax];  // Mass     of i^th final state hadron
  int    brNp         = 0;  // Number of final state p
  int    brNn         = 0;  // Number of final state n
  int    brNpip       = 0;  // Number of final state pi+
  int    brNpim       = 0;  // Number of final state pi-
  int    brNpi0       = 0;  // Number of final state pi0

  //-- open output file & create output summary tree & create the tree branches
  //
  LOG("gntpc", pNOTICE)
       << "*** Saving summary tree to: " << gOptOutFileName;
  TFile fout(gOptOutFileName.c_str(),"recreate");
   
TTree * tEvtTree = new TTree("ginuke","GENIE INuke Summary Tree");
  assert(tEvtTree);

  //-- create tree branches
  //
  tEvtTree->Branch("iev",       &brIEv,          "iev/I"       );
  tEvtTree->Branch("probe",     &brProbe,        "probe/I"     );
  tEvtTree->Branch("tgt" ,      &brTarget,       "tgt/I"       );
  tEvtTree->Branch("ke",        &brKE,           "ke/D"        );
  tEvtTree->Branch("e",         &brE,            "e/D"         );
  tEvtTree->Branch("p",         &brP,            "p/D"         );
  tEvtTree->Branch("A",         &brTgtA,         "A/I"         );
  tEvtTree->Branch("Z",         &brTgtZ,         "Z/I"         );
  tEvtTree->Branch("vtxx",      &brVtxX,         "vtxx/D"      );
  tEvtTree->Branch("vtxy",      &brVtxY,         "vtxy/D"      );
  tEvtTree->Branch("vtxz",      &brVtxZ,         "vtxz/D"      );
  tEvtTree->Branch("probe_fsi", &brProbeFSI,     "probe_fsi/I" );
  tEvtTree->Branch("dist",      &brDist,         "dist/D"      );
  tEvtTree->Branch("nh",        &brNh,           "nh/I"        );
  tEvtTree->Branch("pdgh",      brPdgh,          "pdgh[nh]/I " );
  tEvtTree->Branch("Eh",        brEh,            "Eh[nh]/D"    );
  tEvtTree->Branch("ph",        brPh,            "ph[nh]/D"    );
  tEvtTree->Branch("pxh",       brPxh,           "pxh[nh]/D"   );
  tEvtTree->Branch("pyh",       brPyh,           "pyh[nh]/D"   );
  tEvtTree->Branch("pzh",       brPzh,           "pzh[nh]/D"   );
  tEvtTree->Branch("cth",       brCosth,         "cth[nh]/D"   );
  tEvtTree->Branch("mh",        brMh,            "mh[nh]/D"    );
  tEvtTree->Branch("np",        &brNp,           "np/I"        );
  tEvtTree->Branch("nn",        &brNn,           "nn/I"        );
  tEvtTree->Branch("npip",      &brNpip,         "npip/I"      );
  tEvtTree->Branch("npim",      &brNpim,         "npim/I"      );
  tEvtTree->Branch("npi0",      &brNpi0,         "npi0/I"      );

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
    brIEv = iev; 
    er_tree->GetEntry(iev);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    // analyze current event and fill the summary ntuple

    // clean-up arrays
    //
    for(int j=0; j<kNPmax; j++) {
       brPdgh[j] = 0;
       brEh  [j] = 0;
       brPxh [j] = 0;
       brPyh [j] = 0;
       brPzh [j] = 0;
       brMh  [j] = 0;
    }

    //
    // convert the current event
    //

    GHepParticle * probe  = event.Particle(0);
    GHepParticle * target = event.Particle(1);
    assert(probe && target);

    brProbe    = probe  -> Pdg();
    brTarget   = target -> Pdg();
    brKE       = probe  -> KinE();
    brE        = probe  -> E();
    brP        = probe  -> P4()->Vect().Mag();
    brTgtA     = pdg::IonPdgCodeToA(target->Pdg()); 
    brTgtZ     = pdg::IonPdgCodeToZ(target->Pdg());
    brVtxX     = probe  -> Vx();
    brVtxY     = probe  -> Vy();
    brVtxZ     = probe  -> Vz();
    brProbeFSI = probe  -> RescatterCode(); 
    GHepParticle * rescattered_hadron  = event.Particle(probe->FirstDaughter());
    assert(rescattered_hadron);
    if(rescattered_hadron->Status() == kIStStableFinalState) {
        brDist = -1; // hadron escaped nucleus before interacting;
    }
    else {
      double x  = rescattered_hadron->Vx();
      double y  = rescattered_hadron->Vy();
      double z  = rescattered_hadron->Vz();
      double d2 = TMath::Power(brVtxX-x,2) +
                  TMath::Power(brVtxY-y,2) +
                  TMath::Power(brVtxZ-z,2);
      brDist = TMath::Sqrt(d2);
    }

    brNp       = 0;
    brNn       = 0;
    brNpip     = 0;
    brNpim     = 0;
    brNpi0     = 0;

    int i=0;
    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
       if(pdg::IsPseudoParticle(p->Pdg())) continue;
       if(p->Status() != kIStStableFinalState) continue;

       brPdgh[i] = p->Pdg();
       brEh  [i] = p->E();
       brPxh [i] = p->Px();
       brPyh [i] = p->Py();
       brPzh [i] = p->Pz();
       brPh  [i] =
         TMath::Sqrt(brPxh[i]*brPxh[i]+brPyh[i]*brPyh[i]
                     +brPzh[i]*brPzh[i]);
       brCosth[i] = brPzh[i]/brPh[i];
       brMh  [i] = p->Mass();

       if ( p->Pdg() == kPdgProton  ) brNp++;
       if ( p->Pdg() == kPdgNeutron ) brNn++;
       if ( p->Pdg() == kPdgPiP     ) brNpip++;
       if ( p->Pdg() == kPdgPiM     ) brNpim++;
       if ( p->Pdg() == kPdgPi0     ) brNpi0++;

       i++;
    }
    brNh = i;
    
    ///////////////Test Code///////////////////////
    int tempProbeFSI = brProbeFSI;
    brProbeFSI = HAProbeFSI(tempProbeFSI, brProbe, brNh, brEh, brPdgh, brNpip, brNpim, brNpi0);
    //////////////End Test///////////////////////// 


    // fill the summary tree
    tEvtTree->Fill();

    mcrec->Clear();

  } // event loop

  fin.Close();

  fout.Write();
  fout.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//____________________________________________________________________________________
// FUNCTIONS FOR PARSING CMD-LINE ARGUMENTS 
//____________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  // Common run options. 
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // get input ROOT file (containing a GENIE GHEP event tree)
  if( parser.OptionExists('i') ) {
    LOG("gntpc", pINFO) << "Reading input filename";
    gOptInpFileName = parser.ArgAsString('i');
  } else {
    LOG("gntpc", pFATAL)
       << "Unspecified input filename - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // check input GENIE ROOT file
  bool inpok = !(gSystem->AccessPathName(gOptInpFileName.c_str()));
  if (!inpok) {
    LOG("gntpc", pFATAL)
        << "The input ROOT file ["
        << gOptInpFileName << "] is not accessible";
    gAbortingInErr = true;
    exit(2);
  }

  // get output file format
  if( parser.OptionExists('f') ) {
    LOG("gntpc", pINFO) << "Reading output file format";
    string fmt = parser.ArgAsString('f');

         if (fmt == "gst")                   { gOptOutFileFormat = kConvFmt_gst;                   }
    else if (fmt == "gxml")                  { gOptOutFileFormat = kConvFmt_gxml;                  }
    else if (fmt == "ghep_mock_data")        { gOptOutFileFormat = kConvFmt_ghep_mock_data;        }
    else if (fmt == "rootracker")            { gOptOutFileFormat = kConvFmt_rootracker;            }
    else if (fmt == "rootracker_mock_data")  { gOptOutFileFormat = kConvFmt_rootracker_mock_data;  }
    else if (fmt == "t2k_rootracker")        { gOptOutFileFormat = kConvFmt_t2k_rootracker;        }
    else if (fmt == "numi_rootracker")       { gOptOutFileFormat = kConvFmt_numi_rootracker;       }
    else if (fmt == "t2k_tracker")           { gOptOutFileFormat = kConvFmt_t2k_tracker;           }
    else if (fmt == "nuance_tracker" )       { gOptOutFileFormat = kConvFmt_nuance_tracker;        }
    else if (fmt == "ghad")                  { gOptOutFileFormat = kConvFmt_ghad;                  }
    else if (fmt == "ginuke")                { gOptOutFileFormat = kConvFmt_ginuke;                }
    else                                     { gOptOutFileFormat = kConvFmt_undef;                 }

    if(gOptOutFileFormat == kConvFmt_undef) {
      LOG("gntpc", pFATAL) << "Unknown output file format (" << fmt << ")";
      gAbortingInErr = true;
      exit(3);
    }

  } else {
    LOG("gntpc", pFATAL) << "Unspecified output file format";
    gAbortingInErr = true;
    exit(4);
  }

  // get output file name 
  if( parser.OptionExists('o') ) {
    LOG("gntpc", pINFO) << "Reading output filename";
    gOptOutFileName = parser.ArgAsString('o');
  } else {
    LOG("gntpc", pINFO)
       << "Unspecified output filename - Using default";
    gOptOutFileName = DefaultOutputFile();
  }

  // get number of events to convert
  if( parser.OptionExists('n') ) {
    LOG("gntpc", pINFO) << "Reading number of events to analyze";
    gOptN = parser.ArgAsInt('n');
  } else {
    LOG("gntpc", pINFO)
       << "Unspecified number of events to analyze - Use all";
    gOptN = -1;
  }

  // get format version number
  if( parser.OptionExists('v') ) {
    LOG("gntpc", pINFO) << "Reading format version number";
    gOptVersion = parser.ArgAsInt('v');
    LOG("gntpc", pINFO)
       << "Using version number: " << gOptVersion;
  } else {
    LOG("gntpc", pINFO)
       << "Unspecified version number - Use latest";
    gOptVersion = LatestFormatVersionNumber();
    LOG("gntpc", pINFO)
       << "Latest version number: " << gOptVersion;
  }

  // check whether to copy MC job metadata (only if output file is in ROOT format)
  gOptCopyJobMeta = parser.OptionExists('c');

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gntpc", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gntpc", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }
 
  LOG("gntpc", pNOTICE) << "Input filename  = " << gOptInpFileName;
  LOG("gntpc", pNOTICE) << "Output filename = " << gOptOutFileName;
  LOG("gntpc", pNOTICE) << "Conversion to format = " << gOptRanSeed 
                        << ", vrs = " << gOptVersion;
  LOG("gntpc", pNOTICE) << "Number of events to be converted = " << gOptN;
  LOG("gntpc", pNOTICE) << "Copy metadata? = " << ((gOptCopyJobMeta) ? "Yes" : "No");
  LOG("gntpc", pNOTICE) << "Random number seed = " << gOptRanSeed;

  LOG("gntpc", pNOTICE) << *RunOpt::Instance();
}
//____________________________________________________________________________________
string DefaultOutputFile(void)
{
  // filename extension - depending on file format
  string ext="";
  if      (gOptOutFileFormat == kConvFmt_gst                  ) { ext = "gst.root";         }
  else if (gOptOutFileFormat == kConvFmt_gxml                 ) { ext = "gxml";             }
  else if (gOptOutFileFormat == kConvFmt_ghep_mock_data       ) { ext = "mockd.ghep.root";  }
  else if (gOptOutFileFormat == kConvFmt_rootracker           ) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_rootracker_mock_data ) { ext = "mockd.gtrac.root"; }
  else if (gOptOutFileFormat == kConvFmt_t2k_rootracker       ) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_numi_rootracker      ) { ext = "gtrac.root";       }
  else if (gOptOutFileFormat == kConvFmt_t2k_tracker          ) { ext = "gtrac.dat";        }
  else if (gOptOutFileFormat == kConvFmt_nuance_tracker       ) { ext = "gtrac_legacy.dat"; }
  else if (gOptOutFileFormat == kConvFmt_ghad                 ) { ext = "ghad.dat";         }
  else if (gOptOutFileFormat == kConvFmt_ginuke               ) { ext = "ginuke.root";      }

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
//____________________________________________________________________________________
int LatestFormatVersionNumber(void)
{
  if      (gOptOutFileFormat == kConvFmt_gst                  ) return 1;
  else if (gOptOutFileFormat == kConvFmt_gxml                 ) return 1;
  else if (gOptOutFileFormat == kConvFmt_ghep_mock_data       ) return 1;
  else if (gOptOutFileFormat == kConvFmt_rootracker           ) return 1;
  else if (gOptOutFileFormat == kConvFmt_rootracker_mock_data ) return 1;
  else if (gOptOutFileFormat == kConvFmt_t2k_rootracker       ) return 1;
  else if (gOptOutFileFormat == kConvFmt_numi_rootracker      ) return 1;
  else if (gOptOutFileFormat == kConvFmt_t2k_tracker          ) return 2;
  else if (gOptOutFileFormat == kConvFmt_nuance_tracker       ) return 1;
  else if (gOptOutFileFormat == kConvFmt_ghad                 ) return 1;
  else if (gOptOutFileFormat == kConvFmt_ginuke               ) return 1;

  return -1;
}
//____________________________________________________________________________________
void PrintSyntax(void)
{
  string basedir  = string( gSystem->Getenv("GENIE") );
  string thisfile = basedir + string("/src/Apps/gNtpConv.cxx");
  string cmd      = "less " + thisfile;

  gSystem->Exec(cmd.c_str());
}
//____________________________________________________________________________________
/* Converting HN probe_fsi to HA probe_fsi */
int HAProbeFSI(int probe_fsi, int probe_pdg, int numh, double E_had[], int pdg_had[], int numpip, int numpim, int numpi0)
{
  int index = -1;
  double energy = 0;

  for(int i=0; i<numh; i++)
    { energy += E_had[i]; }


// Determine fates (as defined in Intranuke/INukeUtils.cxx/ utils::intranuke::FindhAFate())
  if (probe_fsi==3 && numh==1) // Elastic
    { index=3; }
  else if (energy==E_had[0] && numh==1) // No interaction
    { index=1; }
  else if ( pdg::IsPion(probe_pdg) && numpip+numpi0+numpim==0) // Absorption
    { index=5; }
  else if ( (pdg::IsNucleon(probe_pdg) && numpip+numpi0+numpim==0 && numh>2 )
	    || (probe_pdg==kPdgGamma && energy!=E_had[0] && numpip+numpi0+numpim==0)) // Knock-out
    { index=6; }
  else if ( numpip+numpi0+numpim > (pdg::IsPion(probe_pdg) ? 1 : 0) ) // Pion production
    { index=7; }
  else if ( numh>=2 ) // Inelastic or Charge Exchange
    {
      for(int i = 0; i < numh; i++)
	{
	  if ( (pdg::IsPion(probe_pdg) && ( probe_pdg==pdg_had[i] ))
	       || pdg::IsNucleon(probe_pdg) ) 
	    {index=4;}	
	  if(index!=4) 
	    {index=2;}
	}
    }
      else //Double Charge Exchange or Undefined
	{
	  bool undef = true;
	  if ( pdg::IsPion(probe_pdg) )
	    {
	      for (int iter = 0; iter < numh; iter++)
		{
		  if      (probe_pdg==211 && pdg_had[iter]==-211) { index=8; undef=false; }
		  else if (probe_pdg==-211 && pdg_had[iter]==211) { index=8; undef=false; }
		}
	    }
	  if (undef) { index=0; }
	}

  return index;
}
