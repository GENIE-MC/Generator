//____________________________________________________________________________
/*!

\program gNtpConv

\brief   Converts a GENIE neutrino events (from ER-type GENIE ROOT Trees) to a 
         variety of textual formats (typically used in legacy systems) or in 
         XML format.

         Syntax:
              gntpc -i input_filename [-o output_filename] [-f fmt]

         Options :
           [] denotes an optional argument

           -f specifies the output text file format. Default is 0.
		0 : GENIE-style XML file
		1 : GENIE's GHEP in tabular text format
		2 : GMINOS-style text file
		3 : NUANCE-style text file

           -o specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
                0 -> *.gxml 
                1 -> *.gtab
                2 -> *.gminos
                3 -> *.nuance
		
\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 23, 2005
*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

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
using std::ostringstream;
using std::ofstream;
using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

using namespace genie;

void   ConvertToGXML      (ofstream & out, EventRecord & event);
void   AddGXMLHeader      (ofstream & out);
void   AddGXMLFooter      (ofstream & out);
void   ConvertToGTab      (ofstream & out, EventRecord & event);
void   ConvertToGMinos    (ofstream & out, EventRecord & event);
void   ConvertToNuance    (ofstream & out, EventRecord & event);
int    GHepToNuanceIst    (GHepParticle * p);
int    GHep2NuancePDGC    (GHepParticle * p);
void   GetCommandLineArgs (int argc, char ** argv);
void   PrintSyntax        (void);
string DefaultOutputFile  (void);

//input options (from command line arguments):
string gOptInpFileName;
string gOptOutFileName;
int    gOptOutFileFormat;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);

  string out_format="";
  if      (gOptOutFileFormat==0) out_format = "GXML";
  else if (gOptOutFileFormat==1) out_format = "GTAB";
  else if (gOptOutFileFormat==2) out_format = "GMINOS";
  else if (gOptOutFileFormat==3) out_format = "NUANCE";

  LOG("gntpc", pNOTICE)
   << "\n\n Converting:...\n"
   << "  GENIE ER ROOT Tree -> " << out_format << "-format file \n";
   
  LOG("gntpc", pNOTICE)
       << "Input  GENIE/ROOT file: " << gOptInpFileName;
  LOG("gntpc", pNOTICE)
                  << "Output file: " << gOptOutFileName;

  //-- open the ROOT file and get the TTree & its header

  TFile file(gOptInpFileName.c_str(),"READ");

  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  //-- Figure out the TTree format (GENIE supports multiple formats)
  //   This program only translates ER ntupes.
  //   Assert that the user has input a correct ntuple type
  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord);

  //-- The ER ntuple contains a single TBranch with NtpMCEventRecord
  //   objects in its leaves. Set the branch address for reading it.
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- Create an ofstream for saving the converted GHEP records
  ofstream output(gOptOutFileName.c_str(), ios::out);

  //-- Add a file header when it is needed
  if(gOptOutFileFormat==0) AddGXMLHeader(output);

  //-- loop over TTree NtpMC records, get each events & call the right
  //   converter depending on the requested format
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    if      (gOptOutFileFormat==0) ConvertToGXML  (output, event);
    else if (gOptOutFileFormat==1) ConvertToGTab  (output, event);
    else if (gOptOutFileFormat==2) ConvertToGMinos(output, event);
    else if (gOptOutFileFormat==3) ConvertToNuance(output, event);
  }

  //-- Add a file footer when it is needed
  if(gOptOutFileFormat==0) AddGXMLHeader(output);

  LOG("gntpc", pINFO) << "\nDone converting GENIE's ER ntuple";

  return 0;
}
//___________________________________________________________________
// FUNCTIONS FOR CONVERTING:
//      ***** GENIE ER ROOT TREE -> GENIE-STYLE XML FILE ****
//___________________________________________________________________
void AddGXMLHeader(ofstream & output)
{
  output << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  output << endl << endl;
  output << "<!-- generated by GENIE gntpc utility -->";
  output << endl << endl;

  output << "<genie_event_list version=\"1.00\">" << endl;
}
//___________________________________________________________________
void ConvertToGXML(ofstream & output, EventRecord & event)
{
  TIter event_iter(&event);

  output << endl << endl;
  output << "  <!-- GENIE GHEP event -->" << endl;
  output << "  <ghep np=\"" << event.GetEntries() 
         << "\" unphysical=\"" 
         << (event.IsUnphysical() ? "true" : "false") << "\">" << endl;

  output << setiosflags(ios::scientific);

  //-- write-out the event-wide properties

  output << "   ";
  output << "  <!-- event weight   -->";
  output << " <wgt> " << event.GetWeight()   << " </wgt>";
  output << endl;

  output << "   ";
  output << "  <!-- cross sections -->";
  output << " <xsec_evnt> " << event.GetXSec()     << " </xsec_evnt>";
  output << " <xsec_kine> " << event.GetDiffXSec() << " </xsec_kine>";
  output << endl;

  output << "   ";
  output << "  <!-- event vertex   -->";
  output << " <vx> " << event.Vertex()->X() << " </vx>";
  output << " <vy> " << event.Vertex()->Y() << " </vy>";
  output << " <vz> " << event.Vertex()->Z() << " </vz>";
  output << " <vt> " << event.Vertex()->T() << " </vt>";
  output << endl;

  //-- write-out the generated particle list
  output << "     <!-- particle list  -->" << endl;
  unsigned int i=0;
  GHepParticle * p = 0;
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

    string type = "U";
    if(p->IsFake()) type = "F";
    else {
      if(p->IsParticle()) type = "P";
      if(p->IsNucleus() ) type = "N";
    }
    output << "     <entry idx=\"" << i << "\" type=\"" 
                                         << type << "\">" << endl;
    output << "        ";
    output << " <pdg> " << p->PdgCode()       << " </pdg>";
    output << " <ist> " << p->Status()        << " </ist>";
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

    output << "     </entry>" << endl;
    i++;
  }
  output << "  </ghep>" << endl;
}
//___________________________________________________________________
void AddGXMLFooter(ofstream & output)
{
  output << endl << endl;
  output << "<genie_event_list version=\"1.00\">";
}
//___________________________________________________________________
// FUNCTIONS FOR CONVERTING:
//     ***** GENIE ER ROOT TREE -> GHEPs in TABULAR TXT FMT *****
//___________________________________________________________________
void ConvertToGTab(ofstream & output, EventRecord & event)
{
  TIter event_iter(&event);

  output << endl << endl;
  output << "# event " << endl;
  output << "$np     : " << event.GetEntries() << endl;
  output << "$unphys : "
         << (event.IsUnphysical() ? "true" : "false") << endl;

  output << setiosflags(ios::scientific);

  //-- write-out the event-wide properties
  output << "$weight : " << event.GetWeight()   << endl;
  output << "$xsec-ev: " << event.GetXSec()     << endl;
  output << "$xsec-kn: " << event.GetDiffXSec() << endl;
  output << "$ev-vtx : "; 
  output << event.Vertex()->X();
  output << setfill(' ') << setw(16) << event.Vertex()->Y();
  output << setfill(' ') << setw(16) << event.Vertex()->Z();
  output << setfill(' ') << setw(16) << event.Vertex()->T() << endl;

  //-- write-out the generated particle list
  unsigned int i=0;
  GHepParticle * p = 0;
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

    string type = "U";
    if(p->IsFake()) type = "F";
    else {
      if(p->IsParticle()) type = "P";
      if(p->IsNucleus() ) type = "N";
    }
    output << "$entry  : "; 
    output << setfill(' ') << setw(4)  << i;
    output << setfill(' ') << setw(12) << p->PdgCode();
    output << setfill(' ') << setw(3)  << type;
    output << setfill(' ') << setw(4)  << p->Status();
    output << setfill(' ') << setw(4)  << p->FirstMother();
    output << setfill(' ') << setw(4)  << p->LastMother();
    output << setfill(' ') << setw(4)  << p->FirstDaughter();
    output << setfill(' ') << setw(4)  << p->LastDaughter();
    output << endl;

    output << setfill(' ') << setw(18) << p->Px();
    output << setfill(' ') << setw(18) << p->Py();
    output << setfill(' ') << setw(18) << p->Pz();
    output << setfill(' ') << setw(18) << p->E();
    output << endl;

    output << setfill(' ') << setw(18) << p->Vx();
    output << setfill(' ') << setw(18) << p->Vy();
    output << setfill(' ') << setw(18) << p->Vz();
    output << setfill(' ') << setw(18) << p->Vt();
    output << endl;

    if(p->PolzIsSet()) {
      output << setfill(' ') << setw(18) 
             << setprecision(5) << p->PolzPolarAngle();
      output << setfill(' ') << setw(18) 
             << setprecision(5) << p->PolzAzimuthAngle();
      output << endl;
    } 
    i++;
  }
}
//___________________________________________________________________
// FUNCTIONS FOR CONVERTING:
//    ***** GENIE ER ROOT TREE -> GMINOS-STYLE TEXT FILE ****
//___________________________________________________________________
void ConvertToGMinos(ofstream & output, EventRecord & event)
{
  LOG("gntpc", pERROR) << "GMINOS format is not supported yet";
  exit(4);
}
//___________________________________________________________________
// FUNCTIONS FOR CONVERTING:
//    ***** GENIE ER ROOT TREE -> NUANCE-STYLE TEXT FILE ****
//___________________________________________________________________
void ConvertToNuance(ofstream & output, EventRecord & event)
{
  TIter event_iter(&event);

  // Nuance begin tag
  output << "$ begin" << endl;

  // Nuance event type
  int nuance_event_type = 0;
  output << "$ nuance " << nuance_event_type << endl;

  // Nuance vertex info
  double vtxx = 0, vtxy = 0, vtxz = 0, vtxt = 0;
  output << "$ vertex " << vtxx << " " << vtxy
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
       output << "$ info 2 949000 0.0000E+00" << endl;
       info_added = true;
     }

     // Add track
     output << "$ track " << pdgc << " " << E << " "
            << dcosx << " " << dcosy << " " << dcosz << " "
            << ist << endl;
  }
  // Nuance end tag
  output << "$ end" << endl;
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
// FUNCTIONS FOR PARSING CMD-LINE ARGUMENTS & PRINTING SYNTAX ON ERR
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

  // check input GENIE ROOT file
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
      LOG("gntpc", pINFO)
               << "Unspecified output file format - Using default";
      gOptOutFileFormat = 0;
    }
  }

  // check output file format
  bool fmtok = (gOptOutFileFormat>=0 && gOptOutFileFormat<=3);
  if (!fmtok) {
    LOG("gntpc", pFATAL)
        << "Invalid output format [" << gOptOutFileFormat << "]";
    PrintSyntax();
    exit(3);
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
}
//___________________________________________________________________
string DefaultOutputFile(void)
{
  // filename extension - depending on file format
  string ext="";
  if      (gOptOutFileFormat==0) ext = "gxml";
  else if (gOptOutFileFormat==1) ext = "gtab";
  else if (gOptOutFileFormat==2) ext = "gminos";
  else if (gOptOutFileFormat==3) ext = "nuance";

  string inpname = gOptInpFileName;
  unsigned int L = inpname.length();

  // if the last 4 characters are "root" (ROOT file extension) then
  // remove them
  if(inpname.substr(L-4, L).find("root") != string::npos) {
    inpname.erase(L-4, L);
  }

  ostringstream name;
  name << inpname << ext;

  return gSystem->BaseName(name.str().c_str());
}
//___________________________________________________________________
void PrintSyntax(void)
{
  LOG("gntpc", pNOTICE)
    << "\n\n" 
    << "Syntax: \n"
    << "  gntpc -i input_filename [-o output_filename] [-f fmt] \n"
    << "\n"
    << "\n"
    << "  Options : \n"
    << "       [] denotes an optional argument \n"
    << "       -f specifies the output text file format. Default is 0.\n"
    << "	        0 : GENIE-style XML file \n"
    << "		1 : GENIE's GHEP in tabular text format \n"
    << "		2 : GMINOS-style text file \n"
    << "		3 : NUANCE-style text file \n\n"
    << "       -o specifies the output filename. \n"
    << "          If not specified a the default filename is constructed by the  \n"
    << "          input base name and an extension depending on the file format: \n"
    << "                0 -> *.gxml \n"
    << "                1 -> *.gtab \n"
    << "                2 -> *.gminos \n"
    << "                3 -> *.nuance \n" << ENDL;
}
//___________________________________________________________________
