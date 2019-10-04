//_____________________________________________________________________________________________
/*!

\program gevpick

\brief   Reads a list of GENIE event files (GHEP format), `cherry-picks' events with a given 
         topology and writes them out in a separate file. The output event tree contains two
         additional branches to aid book-keeping by maintaining a link to the source location
         of each cherry-picked event. For each such event we store a) the name of the original
         file and b) its original event number.

         This is the _only_recommended_ way to obtain event files that contain specific final 
         states (by cherry-picking events from files generated running GENIE in a comprehensive 
         mode). We don't recommend you attempt switching off generator-level reaction modes.
         No detector measures generator-level reaction modes like CCQE or NCRES.
         Detectors measure final states / topologies like {1mu-,0pi}, {1mu-,1pi+},
         {0mu-, 1pi0}, {1 track, 1 shower}, {1 mu-like ring} etc depending on granularity, 
         thresholds and PID capabilities.
         No final state / topology is a proxy for any particular reaction mode (and vice versa).
         Intranuclear re-scattering in particular causes significant migration between states
         (see Table 8.1 in the Physics and User manual).
         Examples:
         - {1mu-,0pi} is mostly numuCCQE but this particular final state can also come about 
           by numu resonance production followed by pion absorption.
         - numuCCQE yields mostly {1mu-,0pi} final states but occasionaly can yield {1mu-,1pi} 
           if the recoil nucleon re-interacts.
         - NC1pi0 final states can be caused by all 
           a) NC elastic followed by nucleon rescattering,
           b) NC resonance neutrino-production, 
           c) NC non-resonance background, 
           d) low-W NC DIS
           e) NC coherent scattering. 
           Each such NC1pi0 source contributes differently to the pion momentum distribution.

         Synopsis:
           gevpick -i list_of_input_files -t topology  
                   [-o output_file]
                   [--message-thresholds xmfile]
                   [--event-record-print-level level]

         Options:

           [] denotes an optional argument

           -i 
              Specify input file(s).
              Wildcards accepted, eg `-i "/data/genie/t2k/gntp.*.ghep.root"'
           -t 
              Specify event topology to cherry-pick.
              The input topology can be any of
                - all 
                    all (basically merges all files into one)
                - numu_cc_1pip 
                    numu CC with 1 \pi^{+} (and no other pion) in final state
                - numu_cc_1pi0 
                    numu CC with 1 \pi^{0} (and no other pion) in final state
                - numu_cc_1pim 
                    numu CC with 1 \pi^{-} (and no other pion) in final state
                - numu_nc_1pip    
                    numu NC with 1 \pi^{+} (and no other pion) in final state
                - numu_nc_1pi0    
                    numu NC with 1 \pi^{0} (and no other pion) in final state
                - numu_nc_1pim    
                    numu NC with 1 \pi^{-} (and no other pion) in final state
                - numu_cc_hyperon 
                    numu CC with at least one hyperon 
                    (\Sigma^{+,0,-}, \Lambda^{0}, \Xi^{0,-}, \Omega^{-}) in final state
                - numubar_cc_hyperon 
                    \bar{numu} CC with at least one hyperon 
                    (\Sigma^{+,0,-}, \Lambda^{0}, \Xi^{0,-}, \Omega^{-}) in final state
                - cc_hyperon         
                    any (anti)neutrino CC with at least one hyperon 
                    (\Sigma^{+,0,-}, \Lambda^{0}, \Xi^{0,-}, \Omega^{-}) in final state
                - <can add more / please send request to costas.andreopoulos \at stfc.ac.uk>
           -o 
              Specify output filename.
              (optional, default: gntp.<topology>.ghep.root)
          --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
          --event-record-print-level
              Allows users to set the level of information shown when the event
              record is printed in the screen. See GHepRecord::Print().

         Examples:

           (1)  % gevpick -i "*.ghep.root" -t numu_nc_1pi0

                Will read all events in all *.ghep.root files and will cherry-pick 
                numu NC 1pi0 events. All cherry-picked events will be saved in the 
                output file gntp.numu_nc_1pi0.ghep.root (default name).

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created August 09, 2010

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_____________________________________________________________________________________________

#include <cassert>
#include <string>
#include <sstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TChainElement.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/RunOpt.h"

using std::string;
using std::ostringstream;

using namespace genie;

// func prototypes
void   GetCommandLineArgs (int argc, char ** argv);
void   RunCherryPicker    (void);
bool   AcceptEvent        (const EventRecord & event);
void   PrintSyntax        (void);
string DefaultOutputFile  (void);

// cherry-picked topologies
typedef enum EGPickTopo {
  kPtUndefined = 0,
  kPtAll,
  kPtNumuCC1pip,
  kPtNumuCC1pi0,
  kPtNumuCC1pim,
  kPtNumuNC1pip,
  kPtNumuNC1pi0,
  kPtNumuNC1pim,
  kPtNumuCChyperon,
  kPtNumubarCChyperon,
  kPtCChyperon

} GPickTopo_t;

// input options (from command line arguments):
string      gOptInpFileNames;  ///< input file name
string      gOptOutFileName;   ///< output file name
GPickTopo_t gPickedTopology;   ///< output file format id

//____________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);

  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  RunCherryPicker();

  return 0;
}
//____________________________________________________________________________________
void RunCherryPicker(void)
{
  // Create an NtpWriter for writing out a tree with the cherry-picked events
  // Add 2 additional branches to the output event tree to save the original filename
  // and the event number in the original file (so that all info can be traced back 
  // to its source).

  NtpWriter ntpw(kNFGHEP, 0);
  ntpw.CustomizeFilename(gOptOutFileName);
  ntpw.Initialize();
  TObjString* brOrigFilename = new TObjString;
  Long64_t    brOrigEvtNum;
  ntpw.EventTree()->Branch("orig_filename", "TObjString", &brOrigFilename, 5000,0);
  ntpw.EventTree()->Branch("orig_evtnum", &brOrigEvtNum, "brOrigEvtNum/L");
  Long64_t iev_glob = 0;

  // Load input trees. More than one trees can be loaded here if a wildcard was
  // specified with -f (eg -f /data/myfiles/genie/*.ghep.root)

  TChain gchain;
  gchain.Add(gOptInpFileNames.c_str());

  TObjArray * file_array = gchain.GetListOfFiles();
  int nfiles = file_array->GetEntries();
  LOG("gevpick", pFATAL) 
      << "Processing " << nfiles
      << (nfiles==1 ? " file " : " files ");

  //
  // Loop over input event files
  //

  TIter next_file(file_array);
  TChainElement *chEl=0;

  while (( chEl=(TChainElement*)next_file() )) {

     TFile fin(chEl->GetTitle(),"read");
     TTree * ghep_tree = 
        dynamic_cast <TTree *> ( fin.Get("gtree")  );

     if(!ghep_tree) {
        LOG("gevpick", pWARN) 
           << "No GHEP tree found in " << chEl->GetTitle();
        LOG("gevpick", pWARN) 
           << "Skipping to next file...";
        continue;
     }

     NtpMCEventRecord * mcrec = 0;
     ghep_tree->SetBranchAddress("gmcrec", &mcrec);
     if (!mcrec) {
       LOG("gevpick", pERROR) << "Null MC record";
       return;
     }
     Long64_t nmax = ghep_tree->GetEntries();
     LOG("gevpick", pNOTICE) 
        << "* Analyzing: " << nmax 
        << " events from GHEP tree in file: " << chEl->GetTitle();

     NtpMCTreeHeader * thdr = 
        dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );
     LOG("gevpick", pNOTICE) 
          << "Input tree header: " << *thdr;

     //
     // Loop over events in current file
     //

     for(Long64_t iev = 0; iev < nmax; iev++) {
       ghep_tree->GetEntry(iev);
       NtpMCRecHeader rec_header = mcrec->hdr;
       EventRecord &  event      = *(mcrec->event);
       LOG("gevpick", pDEBUG) << rec_header;
       LOG("gevpick", pDEBUG) << event;
       if(AcceptEvent(event)) {
          brOrigFilename->SetString(chEl->GetTitle());
          brOrigEvtNum = iev;
          EventRecord * event_copy = new EventRecord(event);
          ntpw.AddEventRecord(iev_glob,event_copy);
          iev_glob++;
       }
       mcrec->Clear();

    } // event loop (current file)
  }// file loop

  // save the cherry-picked MC events
  ntpw.Save();
  
  LOG("gevpick", pFATAL) << "Done!";
}
//____________________________________________________________________________________
bool AcceptEvent(const EventRecord & event)
{
  if ( gPickedTopology == kPtAll       ) return true;
  if ( gPickedTopology == kPtUndefined ) return false;

  const Interaction * interaction = event.Summary();

  int  nupdg     = event.Probe()->Pdg();
  bool isnumu    = (nupdg == kPdgNuMu);
  bool isnumubar = (nupdg == kPdgAntiNuMu);
  bool iscc      = interaction->ProcInfo().IsWeakCC();
  bool isnc      = interaction->ProcInfo().IsWeakNC();

  int NfP        = 0; // number of protons         in final state
  int NfPbar     = 0; // number of anti-protons    in final state
  int NfN        = 0; // number of neutrons        in final state
  int NfNbar     = 0; // number of anti-neutrons   in final state
  int NfPip      = 0; // number of \pi^+'s         in final state
  int NfPim      = 0; // number of \pi^-'s         in final state
  int NfPi0      = 0; // number of \pi^0's         in final state
  int NfKp       = 0; // number of \K^+'s          in final state
  int NfKm       = 0; // number of \K^-'s          in final state
  int NfK0       = 0; // number of \K^0's          in final state
  int NfK0bar    = 0; // number of \bar{\K^0}'s    in final state
  int NfSigmap   = 0; // number of \Sigma^+'s      in final state
  int NfSigma0   = 0; // number of \Sigma^0's      in final state
  int NfSigmam   = 0; // number of \Sigma^-'s      in final state
  int NfLambda0  = 0; // number of \Lambda^0's     in final state
  int NfXi0      = 0; // number of \Xi^0's         in final state
  int NfXim      = 0; // number of \Xi^-'s         in final state
  int NfOmegam   = 0; // number of \Omega^-'s      in final state
  int NfOther    = 0; // number of other particles in final state

  TObjArrayIter piter(&event);
  GHepParticle * p = 0;
  int ip=-1;
  while( (p = (GHepParticle *) piter.Next())) {
    ip++;
    int pdgc = p->Pdg();
    int ist  = p->Status();
    // only final state particles
    if(ist!=kIStStableFinalState) continue;
    // don't count final state lepton as part of the hadronic system
    if(event.Particle(ip)->FirstMother()==0) continue;
    // skip pseudo-particles
    if(pdg::IsPseudoParticle(pdgc)) continue;
    // count ...
    if      (pdgc == kPdgProton     ) NfP++;
    else if (pdgc == kPdgAntiProton ) NfPbar++;
    else if (pdgc == kPdgNeutron    ) NfN++;
    else if (pdgc == kPdgAntiNeutron) NfNbar++;
    else if (pdgc == kPdgPiP        ) NfPip++;
    else if (pdgc == kPdgPiM        ) NfPim++;
    else if (pdgc == kPdgPi0        ) NfPi0++;
    else if (pdgc == kPdgKP         ) NfKp++;
    else if (pdgc == kPdgKM         ) NfKm++;
    else if (pdgc == kPdgK0         ) NfK0++;
    else if (pdgc == kPdgAntiK0     ) NfK0bar++;
    else if (pdgc == kPdgSigmaP     ) NfSigmap++;
    else if (pdgc == kPdgSigma0     ) NfSigma0++;
    else if (pdgc == kPdgSigmaM     ) NfSigmam++;
    else if (pdgc == kPdgLambda     ) NfLambda0++;
    else if (pdgc == kPdgXi0        ) NfXi0++;
    else if (pdgc == kPdgXiM        ) NfXim++;
    else if (pdgc == kPdgOmegaM     ) NfOmegam++;
    else                              NfOther++;
  }

  bool is1pipX  = (NfPip==1 && NfPi0==0 && NfPim==0);
  bool is1pi0X  = (NfPip==0 && NfPi0==1 && NfPim==0);
  bool is1pimX  = (NfPip==0 && NfPi0==0 && NfPim==1);
  bool has_hype = (NfSigmap+NfSigma0+NfSigmam+NfLambda0+NfXi0+NfXim+NfOmegam > 0);

  if ( gPickedTopology == kPtNumuCC1pip ) {
    if(isnumu && iscc && is1pipX) return true;
  }
  else
  if ( gPickedTopology == kPtNumuCC1pi0 ) {
    if(isnumu && iscc && is1pi0X) return true;
  }
  else
  if ( gPickedTopology == kPtNumuCC1pim ) {
    if(isnumu && iscc && is1pimX) return true;
  }
  else
  if ( gPickedTopology == kPtNumuNC1pip ) {
    if(isnumu && isnc && is1pipX) return true;
  }
  else
  if ( gPickedTopology == kPtNumuNC1pi0 ) {
    if(isnumu && isnc && is1pi0X) return true;
  }
  else
  if ( gPickedTopology == kPtNumuNC1pim ) {
    if(isnumu && isnc && is1pimX) return true;
  }
  else
  if ( gPickedTopology == kPtNumuCChyperon ) {
    if(isnumu && iscc && has_hype) return true;
  }
  else
  if ( gPickedTopology == kPtNumubarCChyperon ) {
    if(isnumubar && iscc && has_hype) return true;
  }
  else
  if ( gPickedTopology == kPtCChyperon ) {
    if(iscc && has_hype) return true;
  }

  return false;
}
//____________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  // Common run options. Set defaults and read.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // get input ROOT file (containing a GENIE GHEP event tree)
  if( parser.OptionExists('i') ) {
    gOptInpFileNames = parser.ArgAsString('i');
  } else {
    LOG("gevpick", pFATAL)
       << "Unspecified input filename - Exiting";
    PrintSyntax();
    gAbortingInErr = true;
    exit(1);
  }

  // requested topology
  string topo = "";
  if( parser.OptionExists('t') ) {
    topo = parser.ArgAsString('t');
    if      ( topo == "all"                ) { gPickedTopology = kPtAll;              }
    else if ( topo == "numu_cc_1pip"       ) { gPickedTopology = kPtNumuCC1pip;       }
    else if ( topo == "numu_cc_1pi0"       ) { gPickedTopology = kPtNumuCC1pi0;       }
    else if ( topo == "numu_cc_1pim"       ) { gPickedTopology = kPtNumuCC1pim;       }
    else if ( topo == "numu_nc_1pip"       ) { gPickedTopology = kPtNumuNC1pip;       }
    else if ( topo == "numu_nc_1pi0"       ) { gPickedTopology = kPtNumuNC1pi0;       }
    else if ( topo == "numu_nc_1pim"       ) { gPickedTopology = kPtNumuNC1pim;       }
    else if ( topo == "numu_cc_hyperon"    ) { gPickedTopology = kPtNumuCChyperon;    }
    else if ( topo == "numubar_cc_hyperon" ) { gPickedTopology = kPtNumubarCChyperon; }
    else if ( topo == "cc_hyperon"         ) { gPickedTopology = kPtCChyperon;        }
    else                                     { gPickedTopology = kPtUndefined;        }

    if(gPickedTopology == kPtUndefined) {
      LOG("gevpick", pFATAL) << "Unknown topology (" << topo << ")";
      gAbortingInErr = true;
      exit(1);
    }

  } else {
    LOG("gevpick", pFATAL) << "Unspecified event topology";
    gAbortingInErr = true;
    exit(1);
  }

  // get output file name 
  if( parser.OptionExists('o') ) {
    gOptOutFileName = parser.ArgAsString('o');
  } else {
    LOG("gevpick", pINFO)
       << "Unspecified output filename - Using default";
    gOptOutFileName = DefaultOutputFile();
  }

  // Summarize
  LOG("gevpick", pNOTICE) 
    << "\n\n gevpick job info: "
    << "\n - input file(s)          : " << gOptInpFileNames
    << "\n - output file            : " << gOptOutFileName
    << "\n - cherry-picked topology : " << topo
    << "\n";
}
//____________________________________________________________________________________
string DefaultOutputFile(void)
{
  string tp = "";

  if      (gPickedTopology == kPtAll              ) { tp = "all";                }
  else if (gPickedTopology == kPtNumuCC1pip       ) { tp = "numu_cc_1pip";       }
  else if (gPickedTopology == kPtNumuCC1pi0       ) { tp = "numu_cc_1pi0";       }
  else if (gPickedTopology == kPtNumuCC1pim       ) { tp = "numu_cc_1pim";       }
  else if (gPickedTopology == kPtNumuNC1pip       ) { tp = "numu_nc_1pip";       }
  else if (gPickedTopology == kPtNumuNC1pi0       ) { tp = "numu_nc_1pi0";       }
  else if (gPickedTopology == kPtNumuNC1pim       ) { tp = "numu_nc_1pim";       }
  else if (gPickedTopology == kPtNumuCChyperon    ) { tp = "numu_cc_hyperon";    }
  else if (gPickedTopology == kPtNumubarCChyperon ) { tp = "numubar_cc_hyperon"; }
  else if (gPickedTopology == kPtCChyperon        ) { tp = "cc_hyperon";         } 

  ostringstream fnm;
  fnm << "gntp." << tp << ".ghep.root";

  return fnm.str();
}
//____________________________________________________________________________________
void PrintSyntax(void)
{
  string basedir  = string( gSystem->Getenv("GENIE") );
  string thisfile = basedir + string("/src/Apps/gEvPick.cxx");
  string cmd      = "less " + thisfile;

  gSystem->Exec(cmd.c_str());
}
//____________________________________________________________________________________
