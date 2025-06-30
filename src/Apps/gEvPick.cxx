//_____________________________________________________________________________________________
/*!

\program gevpick

\brief   Reads a list of GENIE event files (GHEP format), `cherry-picks' events with a given 
         final state topology (or true interaction mode) and writes them out in a separate file. 

         The output event tree contains 2 additional branches to aid book-keeping by maintaining 
         a link to the source location of each cherry-picked event. For each such event we store 
         a) the name of the original file, and 
         b) its original event number.

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
           gevpick -i list_of_input_files 
                   -t type
                  [-o output_file]
                  [--message-thresholds xmfile]
                  [--event-record-print-level level]

         Options:

           [] denotes an optional argument

           -i 
              Specify input file(s).
              Wildcards accepted, eg `-i "/data/genie/t2k/gntp.*.ghep.root"'

           -t 
              Specify the type of events to cherry-pick.
              The following options are suported currently:

                - all 
                    all (basically merges all files into one)

              .................................................................................
              Event types based on final-state topology.
              Include all reaction modes contributing to the selected topology.
              .................................................................................

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

              .................................................................................
              Event types based on true reaction mode.
              Includes all events with the specified mode, regardless of the resulting f/s.
              .................................................................................

                - cc_qe 
                - numu_cc_qe 
                    Genuine CCQE
                - cc_mec 
                - numu_cc_mec 
                    Genuine CCMEC
                - cc_qe_mec
                - numu_cc_qe_mec
                    Genuine CCQE or genuine CCMEC
                - not_cc_qe_mec
                - not_numu_cc_qe_mec
                    Anything other than genuine CCQE or genuine CCMEC
                - cc_not_qe_mec
                    CC but not genuine CCQE or genuine CCMEC
                - numu_not_cc_qe_mec
                    numu but not genuine CCQE or genuine CCMEC

              .................................................................................

                <can add more / please send request to c.andreopoulos \at cern.ch>

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

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created August 09, 2010

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
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

// cherry-picked event types
typedef enum EGPickType {
  kPtUndefined = 0,  
  kPtAll,
  //
  // Types based on final state topology, regardless of reaction mode
  //
  kPtTopoNumuCC1pip,
  kPtTopoNumuCC1pi0,
  kPtTopoNumuCC1pim,
  kPtTopoNumuNC1pip,
  kPtTopoNumuNC1pi0,
  kPtTopoNumuNC1pim,
  kPtTopoNumuCChyperon,
  kPtTopoNumubarCChyperon,
  kPtTopoCChyperon,
  //
  // Types based on true reaction mode, regardless of final state topology
  //
  kPtReacModeCCQE,
  kPtReacModeNumuCCQE,
  kPtReacModeCCMEC,
  kPtReacModeNumuCCMEC,
  kPtReacModeCCQEMEC,
  kPtReacModeNumuCCQEMEC,
  kPtReacModeNotCCQEMEC,
  kPtReacModeCCNotQEMEC,
  kPtReacModeNotNumuCCQEMEC,
  kPtReacModeNumuNotCCQEMEC

} GPickType_t;

// input options (from command line arguments):
string      gOptInpFileNames;  ///< input file name
string      gOptOutFileName;   ///< output file name
string      gPickedTypeStr;    ///< output file name
GPickType_t gPickedType;       ///< output file format id

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
  LOG("gevpick", pNOTICE) 
      << "Processing " << nfiles
      << (nfiles==1 ? " file " : " files ");

  //
  // Loop over input event files
  //

  TIter next_file(file_array);
  TChainElement *chEl=0;

  unsigned int total_events  = 0;
  unsigned int picked_events = 0;

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
       total_events++;
       ghep_tree->GetEntry(iev);
       NtpMCRecHeader rec_header = mcrec->hdr;
       EventRecord &  event      = *(mcrec->event);
       LOG("gevpick", pDEBUG) << rec_header;
       LOG("gevpick", pDEBUG) << event;
       if(AcceptEvent(event)) {
         picked_events++;
          brOrigFilename->SetString(chEl->GetTitle());
          brOrigEvtNum = iev;
          ntpw.AddEventRecord( iev_glob, &event );
          iev_glob++;
       }
       mcrec->Clear();

    } // event loop (current file)
  }// file loop

  // save the cherry-picked MC events
  ntpw.Save();
  
  LOG("gevpick", pNOTICE) << "Picked " << picked_events << " / " << total_events << " events of type " << gPickedTypeStr;
  LOG("gevpick", pNOTICE) << "Done!";
}
//____________________________________________________________________________________
bool AcceptEvent(const EventRecord & event)
{
  if ( gPickedType == kPtAll       ) return true;
  if ( gPickedType == kPtUndefined ) return false;

  const Interaction * interaction = event.Summary();

  int  nupdg     = event.Probe()->Pdg();
  bool isnumu    = (nupdg == kPdgNuMu);
  bool isnumubar = (nupdg == kPdgAntiNuMu);
  bool iscc      = interaction->ProcInfo().IsWeakCC();
  bool isnc      = interaction->ProcInfo().IsWeakNC();
  bool isqe      = interaction->ProcInfo().IsQuasiElastic();
  bool ismec     = interaction->ProcInfo().IsMEC();
  bool isstr     = interaction->ExclTag().IsStrangeEvent();
  bool ischm     = interaction->ExclTag().IsCharmEvent();

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

  if ( gPickedType == kPtTopoNumuCC1pip ) {
    if(isnumu && iscc && is1pipX) return true;
  }
  else
  if ( gPickedType == kPtTopoNumuCC1pi0 ) {
    if(isnumu && iscc && is1pi0X) return true;
  }
  else
  if ( gPickedType == kPtTopoNumuCC1pim ) {
    if(isnumu && iscc && is1pimX) return true;
  }
  else
  if ( gPickedType == kPtTopoNumuNC1pip ) {
    if(isnumu && isnc && is1pipX) return true;
  }
  else
  if ( gPickedType == kPtTopoNumuNC1pi0 ) {
    if(isnumu && isnc && is1pi0X) return true;
  }
  else
  if ( gPickedType == kPtTopoNumuNC1pim ) {
    if(isnumu && isnc && is1pimX) return true;
  }
  else
  if ( gPickedType == kPtTopoNumuCChyperon ) {
    if(isnumu && iscc && has_hype) return true;
  }
  else
  if ( gPickedType == kPtTopoNumubarCChyperon ) {
    if(isnumubar && iscc && has_hype) return true;
  }
  else
  if ( gPickedType == kPtTopoCChyperon ) {
    if(iscc && has_hype) return true;
  }
  else 
  if ( gPickedType == kPtReacModeCCQE ) {
    if(isstr || ischm) return false;
    if(iscc && isqe) return true;
  }
  else 
  if ( gPickedType == kPtReacModeNumuCCQE ) {
    if(isstr || ischm) return false;
    if(isnumu && iscc && isqe) return true;
  }
  else 
  if ( gPickedType == kPtReacModeCCMEC ) {
    if(isstr || ischm) return false;
    if(iscc && ismec) return true;
  }
  else 
  if ( gPickedType == kPtReacModeNumuCCMEC ) {
    if(isstr || ischm) return false;
    if(isnumu && iscc && ismec) return true;
  }
  else 
  if ( gPickedType == kPtReacModeCCQEMEC ) {
    if(isstr || ischm) return false;
    if(iscc && (isqe || ismec)) return true;
  }
  else 
  if ( gPickedType == kPtReacModeNumuCCQEMEC ) {
    if(isstr || ischm) return false;
    if(isnumu && iscc && (isqe || ismec)) return true;
  }
  else 
  if ( gPickedType == kPtReacModeNotCCQEMEC ) {
    if(isstr || ischm) return false;
    if(!(iscc && (isqe || ismec))) return true;
  }
  else 
  if ( gPickedType == kPtReacModeNotNumuCCQEMEC ) {
    if(isstr || ischm) return false;
    if(!(isnumu && iscc && (isqe || ismec))) return true;
  }
  else 
  if ( gPickedType == kPtReacModeCCNotQEMEC ) {
    if(isstr || ischm) return false;
    if(iscc && !(isqe || ismec)) return true;
  }
  else 
  if ( gPickedType == kPtReacModeNumuNotCCQEMEC ) {
    if(isstr || ischm) return false;
    if(isnumu && !(iscc && (isqe || ismec))) return true;
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

  // requested event type
  string evtype = "";
  if( parser.OptionExists('t') ) {
    evtype = parser.ArgAsString('t');
    if      ( evtype == "all"                ) { gPickedType = kPtAll;                    }

    else if ( evtype == "numu_cc_1pip"       ) { gPickedType = kPtTopoNumuCC1pip;         }
    else if ( evtype == "numu_cc_1pi0"       ) { gPickedType = kPtTopoNumuCC1pi0;         }
    else if ( evtype == "numu_cc_1pim"       ) { gPickedType = kPtTopoNumuCC1pim;         }
    else if ( evtype == "numu_nc_1pip"       ) { gPickedType = kPtTopoNumuNC1pip;         }
    else if ( evtype == "numu_nc_1pi0"       ) { gPickedType = kPtTopoNumuNC1pi0;         }
    else if ( evtype == "numu_nc_1pim"       ) { gPickedType = kPtTopoNumuNC1pim;         }
    else if ( evtype == "numu_cc_hyperon"    ) { gPickedType = kPtTopoNumuCChyperon;      }
    else if ( evtype == "numubar_cc_hyperon" ) { gPickedType = kPtTopoNumubarCChyperon;   }
    else if ( evtype == "cc_hyperon"         ) { gPickedType = kPtTopoCChyperon;          }

    else if ( evtype == "cc_qe"              ) { gPickedType = kPtReacModeCCQE;           }
    else if ( evtype == "numu_cc_qe"         ) { gPickedType = kPtReacModeNumuCCQE;       }
    else if ( evtype == "cc_mec"             ) { gPickedType = kPtReacModeCCMEC;          }
    else if ( evtype == "numu_cc_mec"        ) { gPickedType = kPtReacModeNumuCCMEC;      }
    else if ( evtype == "cc_qe_mec"          ) { gPickedType = kPtReacModeCCQEMEC;        }
    else if ( evtype == "numu_cc_qe_mec"     ) { gPickedType = kPtReacModeNumuCCQEMEC;    }
    else if ( evtype == "not_cc_qe_mec"      ) { gPickedType = kPtReacModeNotCCQEMEC;     }
    else if ( evtype == "cc_not_qe_mec"      ) { gPickedType = kPtReacModeCCNotQEMEC;     }
    else if ( evtype == "not_numu_cc_qe_mec" ) { gPickedType = kPtReacModeNotNumuCCQEMEC; }
    else if ( evtype == "numu_not_cc_qe_mec" ) { gPickedType = kPtReacModeNumuNotCCQEMEC; }

    else                                       { gPickedType = kPtUndefined;              }

    if(gPickedType == kPtUndefined) {
      LOG("gevpick", pFATAL) << "Unknown event type (" << evtype << ")";
      gAbortingInErr = true;
      exit(1);
    }
    gPickedTypeStr = evtype;

  } else {
    LOG("gevpick", pFATAL) << "Unspecified event type";
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
    << "\n - cherry-picked topology : " << evtype
    << "\n";
}
//____________________________________________________________________________________
string DefaultOutputFile(void)
{
  string tp = "";

  if      (gPickedType == kPtAll                    ) { tp = "all";                }

  else if (gPickedType == kPtTopoNumuCC1pip         ) { tp = "numu_cc_1pip";       }
  else if (gPickedType == kPtTopoNumuCC1pi0         ) { tp = "numu_cc_1pi0";       }
  else if (gPickedType == kPtTopoNumuCC1pim         ) { tp = "numu_cc_1pim";       }
  else if (gPickedType == kPtTopoNumuNC1pip         ) { tp = "numu_nc_1pip";       }
  else if (gPickedType == kPtTopoNumuNC1pi0         ) { tp = "numu_nc_1pi0";       }
  else if (gPickedType == kPtTopoNumuNC1pim         ) { tp = "numu_nc_1pim";       }
  else if (gPickedType == kPtTopoNumuCChyperon      ) { tp = "numu_cc_hyperon";    }
  else if (gPickedType == kPtTopoNumubarCChyperon   ) { tp = "numubar_cc_hyperon"; }
  else if (gPickedType == kPtTopoCChyperon          ) { tp = "cc_hyperon";         } 

  else if (gPickedType == kPtReacModeCCQE           ) { tp = "cc_qe";              }
  else if (gPickedType == kPtReacModeNumuCCQE       ) { tp = "numu_cc_qe";         }
  else if (gPickedType == kPtReacModeCCMEC          ) { tp = "cc_mec";             }
  else if (gPickedType == kPtReacModeNumuCCMEC      ) { tp = "numu_cc_mec";        }
  else if (gPickedType == kPtReacModeCCQEMEC        ) { tp = "cc_qe_mec";          }
  else if (gPickedType == kPtReacModeNumuCCQEMEC    ) { tp = "numu_cc_qe_mec";     }
  else if (gPickedType == kPtReacModeNotCCQEMEC     ) { tp = "not_cc_qe_mec";      }
  else if (gPickedType == kPtReacModeCCNotQEMEC     ) { tp = "cc_not_qe_mec";      }
  else if (gPickedType == kPtReacModeNotNumuCCQEMEC ) { tp = "not_numu_cc_qe_mec"; }
  else if (gPickedType == kPtReacModeNumuNotCCQEMEC ) { tp = "numu_not_cc_qe_mec"; }

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
