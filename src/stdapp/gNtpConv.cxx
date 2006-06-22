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
		4 : NEUGEN/GENIE cross generator test-style text file

           -o specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
                0 -> *.gxml 
                1 -> *.gtab
                2 -> *.gminos
                3 -> *.nuance
                4 -> *.gneugen
		
\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 23, 2005

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
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
void   ConvertToGNeugen   (ofstream & out, EventRecord & event);
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
  else if (gOptOutFileFormat==4) out_format = "GNEUGEN";

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

    if      (gOptOutFileFormat==0) ConvertToGXML   (output, event);
    else if (gOptOutFileFormat==1) ConvertToGTab   (output, event);
    else if (gOptOutFileFormat==2) ConvertToGMinos (output, event);
    else if (gOptOutFileFormat==3) ConvertToNuance (output, event);
    else if (gOptOutFileFormat==4) ConvertToGNeugen(output, event);

    mcrec->Clear();
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
             (gist == kIStInitialState || gist == kIStNucleonTarget);

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
   case kIStIntermediateState:        nuance_ist = -2;   break;
   case kIStDecayedState:             nuance_ist = -2;   break;
   case kIStNucleonTarget:            nuance_ist = -1;   break;
   case kIStDISPreFragmHadronicState: nuance_ist = -2;   break;
   case kIStPreDecayResonantState:    nuance_ist = -2;   break;
   case kIStHadronInTheNucleus:       nuance_ist = -2;   break;
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
// FUNCTIONS FOR CONVERTING:
// ***** GENIE ER ROOT TREE -> NEUGEN/GENIE TEST-STYLE TEXT FILE ****
//___________________________________________________________________
void ConvertToGNeugen(ofstream & output, EventRecord & event)
{
  // format --
  // CCNC PROC NU NUC TGT q2 W x y p4nu p4nuc p4fsl p4had
  // 
  // notes --
  // CCNC    : 1->CC, 2->NC
  // PROC    : 1->QEL
  // NU,NUC  : std PDG codes for neutrino and hit nucleon
  // TGT     : MINOS PDG code for target (1AAAZZZ000)
  // q2      : selected momentum transfer (q2<0, in GeV^2)
  // W       : selected hadronic invariant mass (in GeV)
  // x       : selected Bjorken scaling var
  // y       : selected inelasticity
  // q2,W,x,y: are the kinematic variables selected by the kinematic
  //           generators using *off-shell* kinematics, not the ones
  //           that can be computed by the generated particle 4-vectors.
  //           To avoid the ambiguity caused by the value of the 
  //           (off-shell) nucleon mass when translating between 
  //           kinematic variables, only the ones actually selected
  //           would be non-zero (QEL:q2, RES:q2,W, DIS:x,y)
  // p4      : all 4-momentum vector in this order = (px,py,pz,E)

  Interaction * interaction = event.GetInteraction();

  const ProcessInfo & proc_info = interaction->GetProcessInfo();

  int ccnc=-1;
  if      (proc_info.IsWeakCC()) ccnc=1;
  else if (proc_info.IsWeakNC()) ccnc=2;

  int proc=-1;
  if (proc_info.IsQuasiElastic()) proc=1;
  if (proc_info.IsResonant())     proc=2;
  if (proc_info.IsCoherent())     proc=4;

  //-- get kinematics as selected by the corresponding kinematics
  //   generator, ignore other variables to avoid possible definition
  //   incosnsistencies due to the off-shell kinematics

  const Kinematics & kinematics = interaction->GetKinematics();
  double q2=0, W=0, x=0, y=0;

  if(proc==1) {
    q2 = kinematics.q2(true);
  }
  if(proc==2) {
    q2 = kinematics.q2(true);
    W  = kinematics.W(true);
  }
  if(proc==4) {
    x  = kinematics.x(true);
    y  = kinematics.y(true);
    q2 = kinematics.q2(true);
  }

  GHepParticle * neutrino = event.Probe();
  GHepParticle * nucleon  = event.StruckNucleon();
  GHepParticle * target   = event.TargetNucleus();
  GHepParticle * lepton   = event.FinalStatePrimaryLepton();
  GHepParticle * hadsyst  = event.FinalStateHadronicSystem();

  if(proc==1) hadsyst=event.Particle(nucleon->FirstDaughter());

  int neu_pdg = neutrino -> PdgCode();
  int nuc_pdg = (nucleon) ? nucleon -> PdgCode() : 0;
  int tgt_pdg = (target)  ? target  -> PdgCode() : 0;

  output << ccnc << " " << proc << " "
         << neu_pdg << " " 
         << nuc_pdg << " " 
         << tgt_pdg << " "
         << q2 << " " 
         << W  << " " 
         << x  << " " 
         << y  << " "; 

  TLorentzVector * p4neu = neutrino->P4();
  output << p4neu->Px() << " "  << p4neu->Py() << " " 
         << p4neu->Pz() << " "  << p4neu->E()  << " ";

  if(nucleon) {
    TLorentzVector * p4nuc = nucleon->P4();
    output << p4nuc->Px() << " "  << p4nuc->Py() << " " 
           << p4nuc->Pz() << "\n" << p4nuc->E()  << " ";
  } else {
    output << "0.0 0.0 0.0 \n 0.0 ";
  }

  TLorentzVector * p4lep = lepton->P4();
  output << p4lep->Px() << " "  << p4lep->Py() << " " 
         << p4lep->Pz() << " "  << p4lep->E()  << " ";

  if(hadsyst) {
     TLorentzVector * p4had = hadsyst->P4();
     output << p4had->Px()  << "\n" << p4had->Py()  << " " 
            << p4had->Pz()  << " "  << p4had->E()   << " ";
  } else {
    output << "0.0 \n 0.0  0.0  0.0 ";
  }
  output << endl;
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
  bool fmtok = (gOptOutFileFormat>=0 && gOptOutFileFormat<=4);
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
  else if (gOptOutFileFormat==4) ext = "gneugen";

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
    << "		3 : NUANCE-style text file \n"
    << "                4:  NEUGEN/GENIE cross generator test-style text file\n\n"
    << "       -o specifies the output filename. \n"
    << "          If not specified a the default filename is constructed by the  \n"
    << "          input base name and an extension depending on the file format: \n"
    << "                0 -> *.gxml \n"
    << "                1 -> *.gtab \n"
    << "                2 -> *.gminos \n"
    << "                3 -> *.nuance \n" 
    << "                4 -> *.gneugen \n" 
    << ENDL;
}
//___________________________________________________________________
