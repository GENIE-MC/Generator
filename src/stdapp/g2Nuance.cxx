//____________________________________________________________________________
/*!

\program g2nuance

\brief   Reads in a GENIE ROOT Tree and transfors the GHEP event record to
         NUANCE-style text-format event records.

         Syntax:
                g2nuance -i input_filename [-o output_filename]

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 23, 2005
*/
//____________________________________________________________________________

#include <string>
#include <fstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using std::ofstream;
using std::endl;
using std::ios;
using namespace genie;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
void ConvertToNuance    (ofstream & out, EventRecord & event);
int  GHepToNuanceIst    (GHepParticle * p);
int  GHep2NuancePDGC    (GHepParticle * p);

//input options (from command line arguments):
string gOptInpFile;
string gOptOutFile;

//defaults for non-mandatory input options:
string kDefOptOutFile = "genie-events.nuance-text";

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);

  LOG("g2nuance", pNOTICE)
    << "\n\nNeutrino Event Converter:\n"
    << "  GENIE ER Trees in ROOT files -> NUANCE-style text files\n";

  LOG("g2nuance", pNOTICE)
              << "Input  GENIE/ROOT file: " << gOptInpFile;
  LOG("g2nuance", pNOTICE)
              << "Output NUANCE/txt file: " << gOptOutFile;

  //-- open the ROOT file and get the TTree & its header

  TFile file(gOptInpFile.c_str(),"READ");

  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("g2nuance", pINFO) << "Input tree header: " << *thdr;

  //-- Figure out the TTree format (GENIE supports multiple formats)
  //   This program only translates ER ntupes to nuance ascii files.
  //   Assert that the user has input a correct ntuple ytpe
  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord);

  //-- The ER ntuple contains a single TBranch with NtpMCEventRecord
  //   objects in its leaves. Set the branch address for reading it.
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- Open an ofstream for saving Nuance-style records
  ofstream nuance_output(gOptOutFile.c_str(), ios::out);

  //-- loop over TTree NtpMC records, get the events & translate them
  //   write them out into NUANCE-style ascii-based file
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("g2nuance", pINFO) << rec_header;
    LOG("g2nuance", pINFO) << event;

    ConvertToNuance(nuance_output, event);
  }

  LOG("g2nuance", pINFO)
        << "\nDone translating the GENIE's ER ntuple"
                                 << " to NUANCE-style ascii files!";
  return 0;
}
//___________________________________________________________________
void ConvertToNuance(ofstream & nuance_output, EventRecord & event)
{
  TIter event_iter(&event);

  // Nuance begin tag
  nuance_output << "$ begin" << endl;

  // Nuance event type
  int nuance_event_type = 0;
  nuance_output << "$ nuance " << nuance_event_type << endl;

  // Nuance vertex info
  double vtxx = 0, vtxy = 0, vtxz = 0, vtxt = 0;
  nuance_output << "$ vertex " << vtxx << " " << vtxy
                << " " << vtxz << " " << vtxt << " " << endl;

  // Nuance 'tracks', GENIE's equivalent of GHepParticle
  GHepParticle * p = 0;
  bool info_added  = false;

  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

     // Convert GENIE's GHEP pdgc & status to NUANCE's equivalent
     int pdgc = GHep2NuancePDGC(p);
     int ist  = GHepToNuanceIst(p);

     // Get particle's energy & momentum
     TLorentzVector * p4 = p->P4();
     double E  = p4->Energy();
     double Px = p4->Px();
     double Py = p4->Py();
     double Pz = p4->Pz();
     double P  = p4->P();
     // Compute direction cosines
     double dcosx = (P>0) ? Px/P : -999;
     double dcosy = (P>0) ? Py/P : -999;
     double dcosz = (P>0) ? Pz/P : -999;

     GHepStatus_t gist = (GHepStatus_t) p->Status();
     bool is_init =
             (gist == kIStInitialState || gist == kIstNucleonTarget);

     if(!is_init && !info_added) {
       // Add nuance obsolete and flux info (not filled in by
       // GENIE here). Add it once after the initial state particles
       nuance_output << "$ info 2 949000 0.0000E+00" << endl;
       info_added = true;
     }

     // Add track
     nuance_output << "$ track " << pdgc << " " << E << " "
                   << dcosx << " " << dcosy << " " << dcosz << " "
                   << ist << endl;
  }
  // Nuance end tag
  nuance_output << "$ end" << endl;
}
//___________________________________________________________________
int GHepToNuanceIst(GHepParticle * p)
{
// translate GENIE's GHEP Status to Nuance status

  GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();

  int nuance_ist;

  switch (ghep_ist) {
   case kIStInitialState:             nuance_ist = -1;   break;
   case kIStStableFinalState:         nuance_ist =  0;   break;
   case kIstIntermediateState:        nuance_ist = -2;   break;
   case kIstDecayedState:             nuance_ist = -2;   break;
   case kIstNucleonTarget:            nuance_ist = -1;   break;
   case kIstDISPreFragmHadronicState: nuance_ist = -2;   break;
   case kIstPreDecayResonantState:    nuance_ist = -2;   break;
   case kIstHadronInTheNucleus:       nuance_ist = -2;   break;
   case kIStUndefined:                nuance_ist = -999; break;
   default:                           nuance_ist = -999; break;
  }
  return nuance_ist;
}
//___________________________________________________________________
int GHep2NuancePDGC(GHepParticle * p)
{
// For most particles both generators use the standard PDG codes.
// For nuclei GHEP PDGC follows the MINOS-convention: 1AAAZZZ000
// NUANCE is using: ZZZAAA

  int ghep_pdgc   = p->PdgCode();
  int nuance_pdgc = ghep_pdgc;

  if ( p->IsNucleus() ) {
      int Z = pdg::IonPdgCodeToZ(ghep_pdgc);
      int A = pdg::IonPdgCodeToA(ghep_pdgc);

      nuance_pdgc = 1000*Z + A;
  }
  return nuance_pdgc;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  //get input ROOT file (containing a GENIE ER ntuple)
  try {
    LOG("g2nuance", pINFO) << "Reading input filename";
    gOptInpFile = utils::clap::CmdLineArgAsString(argc,argv,'i');
  } catch(utils::clap::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("g2nuance", pFATAL)
                      << "Unspecified input filename - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  //get name for output NUANCE-style data file
  try {
    LOG("g2nuance", pINFO) << "Reading output filename";
    gOptOutFile = utils::clap::CmdLineArgAsString(argc,argv,'o');
  } catch(genie::utils::clap::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("g2nuance", pNOTICE)
                << "Unspecified output filename - Using default";
      gOptOutFile = kDefOptOutFile;
    }
  }

  // check input GENIE ROOT file
  bool ok = !(gSystem->AccessPathName(gOptInpFile.c_str()));
  if (!ok) {
    LOG("g2nuance", pFATAL)
               << "The input ROOT file ["
                            << gOptInpFile << "] is not accessible";
    exit(2);
  }
}
//___________________________________________________________________
void PrintSyntax(void)
{
  LOG("g2nuance", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   g2nuance -i input_filename [-o output_filename]\n";
}
//___________________________________________________________________
