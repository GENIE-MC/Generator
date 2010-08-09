//_____________________________________________________________________________________________
/*!

\program gevpick

\brief   Reads a list of GENIE event files, cherry-picks events with a given topology and 
         writes them out in a separate file (GHEP event tree format).
         Additional branches are added in the output tree so that for each cherry-picked event
         the name of its original event file and its original event number are saved 

         The program can be used to obtain event files with given detector topologies by
         cherry-picking event files produced by running GENIE in its default/comprehensive mode.

         Note that this is the __only_recommended__ way to obtain such samples.

         No detector measures generator-level reaction modes like CCQE or NCRES.
         Detectors measure final states / topologies like {1mu-,0pi}, {1mu-,1pi+},
         {0mu-, 1pi0}, {1 track, 1 shower}, {1 mu-like ring} etc depending on granularity, 
         thresholds and PID capabilities.
         No final state / topology is a proxy for any particular reaction mode (and vice versa).
         Intranuclear re-scattering in particular causes significant migration between states
         (see Table 8.1 in the Physics and User manual).
         Examples:
         - {1mu-,0pi} is mostly numuCCQE, but can also be caused by numu resonance production 
           followed by pion absorption.
         - numuCCQE gives mostly {1mu-,0pi} but occasionaly can give {1mu-,1pi} if the recoil 
           nucleon re-interacts.
         - NC1pi0 final states can be caused by all a) NC elastic followed by nucleon rescattering,
           b) NC resonance neutrino-production, c) NC non-resonance background, d) low-W NC DIS
           e) NC coherent scattering. 
           Each such NC1pi0 source contributes differently in the pion momentum distribution.


         Syntax:
           gevpick -i list_of_input_files -t topology [-o output_file]

         Options :

           [] denotes an optional argument

           -i Specify input file(s).
              Wildcards accepted, eg `-i "/data/genie/t2k/gntp.*.ghep.root"'

           -t Specify event topology to cherry-pick.
              Accepted values:
                - all             :
                - numu_cc_1pip    :
                - numu_cc_1pi0    :
                - numu_cc_1pim    :
                - numu_nc_1pip    :
                - numu_nc_1pi0    :
                - numu_nc_1pim    :
                - numu_cc_hyperon :

           -o Specify output filename.
              (optional, default: gntp.<topology>.ghep.root)

         Examples:

           (1)  % gevpick -i "*.ghep.root" -t numu_nc_1pi0

                Will read all events in all *.ghep.root files and will cherry-pick 
                numu NC 1pi0 events. All cherry-picked events will be saved in the 
                output file gntp.numu_nc_1pi0.ghep.root (default name).

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created August 09, 2010

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
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

#include "Conventions/GBuild.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepUtils.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Ntuple/NtpWriter.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/SystemUtils.h"

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
  kPtNumuCChyperon
} GPickTopo_t;

// input options (from command line arguments):
string      gOptInpFileNames;  ///< input file name
string      gOptOutFileName;   ///< output file name
GPickTopo_t gPickedTopology;   ///< output file format id

//____________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);
  RunCherryPicker();
  return 0;
}
//____________________________________________________________________________________
void RunCherryPicker(void)
{
  // Create an NtpWriter for writing out a tree with the cherry-picked events
  // Add 2 additional branches to the output event tree to save XXX

  NtpWriter ntpw(kNFGHEP, 0);
  ntpw.CustomizeFilename(gOptOutFileName);
  ntpw.Initialize();
  TObjString brOrigFilename;
  Long64_t   brOrigEvtNum;
//  ntpw.EventTree()->Branch("orig_filename", "TObjString", &brOrigFilename, 1024,1);
  ntpw.EventTree()->Branch("orig_evtnum", &brOrigEvtNum, "brOrigEvtNum/I");
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
//          brOrigFilename = chEl->GetTitle();
          brOrigEvtNum   = iev;
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

  int  nupdg  = event.Probe()->Pdg();
  bool isnumu = (nupdg == kPdgNuMu);
  bool iscc   = interaction->ProcInfo().IsWeakCC();
  bool isnc   = interaction->ProcInfo().IsWeakNC();

  int NfP     = 0;
  int NfN     = 0;
  int NfPip   = 0;
  int NfPim   = 0;
  int NfPi0   = 0;
  int NfKp    = 0;
  int NfKm    = 0;
  int NfK0    = 0;
  int NfOther = 0;

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
    if      (pdgc == kPdgProton  || pdgc == kPdgAntiProton)   NfP++;
    else if (pdgc == kPdgNeutron || pdgc == kPdgAntiNeutron)  NfN++;
    else if (pdgc == kPdgPiP) NfPip++;
    else if (pdgc == kPdgPiM) NfPim++;
    else if (pdgc == kPdgPi0) NfPi0++;
    else if (pdgc == kPdgKP)  NfKp++;
    else if (pdgc == kPdgKM)  NfKm++;
    else if (pdgc == kPdgK0 || pdgc == kPdgAntiK0)  NfK0++;
    else NfOther++;
  }

  bool is1pip = (NfPip==1 && NfPi0==0 && NfPim==0 && NfKp==0 && NfKm==0 && NfK0==0);
  bool is1pi0 = (NfPip==0 && NfPi0==1 && NfPim==0 && NfKp==0 && NfKm==0 && NfK0==0);
  bool is1pim = (NfPip==0 && NfPi0==0 && NfPim==1 && NfKp==0 && NfKm==0 && NfK0==0);

  if ( gPickedTopology == kPtNumuCC1pip ) {
    if(isnumu && iscc && is1pip) return true;
    return false;
  }
  if ( gPickedTopology == kPtNumuCC1pi0 ) {
    if(isnumu && iscc && is1pi0) return true;
    return false;
  }
  if ( gPickedTopology == kPtNumuCC1pim ) {
    if(isnumu && iscc && is1pim) return true;
    return false;
  }
  if ( gPickedTopology == kPtNumuNC1pip ) {
    if(isnumu && isnc && is1pip) return true;
    return false;
  }
  if ( gPickedTopology == kPtNumuNC1pi0 ) {
    if(isnumu && isnc && is1pi0) return true;
    return false;
  }
  if ( gPickedTopology == kPtNumuNC1pim ) {
    if(isnumu && isnc && is1pim) return true;
    return false;
  }

  return false;
}
//____________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
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
    if      ( topo == "all"             ) { gPickedTopology = kPtAll;           }
    else if ( topo == "numu_cc_1pip"    ) { gPickedTopology = kPtNumuCC1pip;    }
    else if ( topo == "numu_cc_1pi0"    ) { gPickedTopology = kPtNumuCC1pi0;    }
    else if ( topo == "numu_cc_1pim"    ) { gPickedTopology = kPtNumuCC1pim;    }
    else if ( topo == "numu_nc_1pip"    ) { gPickedTopology = kPtNumuNC1pip;    }
    else if ( topo == "numu_nc_1pi0"    ) { gPickedTopology = kPtNumuNC1pi0;    }
    else if ( topo == "numu_nc_1pim"    ) { gPickedTopology = kPtNumuNC1pim;    }
    else if ( topo == "numu_cc_hyperon" ) { gPickedTopology = kPtNumuCChyperon; }
    else                                  { gPickedTopology = kPtUndefined;     }
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

  if      (gPickedTopology == kPtAll          ) { tp  = "all";              }
  else if (gPickedTopology == kPtNumuCC1pip   ) { tp  = "numu_cc_1pip";     }
  else if (gPickedTopology == kPtNumuCC1pi0   ) { tp  = "numu_cc_1pi0";     }
  else if (gPickedTopology == kPtNumuCC1pim   ) { tp  = "numu_cc_1pim";     }
  else if (gPickedTopology == kPtNumuNC1pip   ) { tp  = "numu_nc_1pip";     }
  else if (gPickedTopology == kPtNumuNC1pi0   ) { tp  = "numu_nc_1pi0";     }
  else if (gPickedTopology == kPtNumuNC1pim   ) { tp  = "numu_nc_1pim";     }
  else if (gPickedTopology == kPtNumuCChyperon) { tp  = "numu_cc_hyperon";  }

  ostringstream fnm;
  fnm << "gntp." << tp << ".ghep.root";

  return fnm.str();
}
//____________________________________________________________________________________
void PrintSyntax(void)
{
  string basedir  = string( gSystem->Getenv("GENIE") );
  string thisfile = basedir + string("/src/stdapp/gEvPick.cxx");
  string cmd      = "less " + thisfile;

  gSystem->Exec(cmd.c_str());
}
//____________________________________________________________________________________
