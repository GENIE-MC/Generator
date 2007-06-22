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
		2 : NUANCE-style text file
		3 : ROOT Tree used for T2K cross-generator comparisons

           -o specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
                0 -> *.gxml 
                1 -> *.gdat
                2 -> *.gnuance
                3 -> *.gt2k.root
		
\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 23, 2005

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
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

void   ConvertToTextFormat();
void   ConvertToT2KMCComparisonsRootFormat();

void   ConvertToGXML      (ofstream & out, EventRecord & event);
void   AddGXMLHeader      (ofstream & out);
void   AddGXMLFooter      (ofstream & out);

void   ConvertToGTab      (ofstream & out, EventRecord & event);

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

const int kNPmax = 100;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);

  if(gOptOutFileFormat==0 || 
     gOptOutFileFormat==1 || 
     gOptOutFileFormat==2) ConvertToTextFormat();

  if(gOptOutFileFormat==3) ConvertToT2KMCComparisonsRootFormat();

  return 0;
}
//___________________________________________________________________
void ConvertToTextFormat()
{
  //-- open the ROOT file and get the TTree & its header
  TFile file(gOptInpFileName.c_str(),"READ");
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord);

  //-- get mc record
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- open the output stream
  ofstream output(gOptOutFileName.c_str(), ios::out);

  //-- add a file header when it is needed
  if(gOptOutFileFormat==0) AddGXMLHeader(output);

  //-- event loop
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    if      (gOptOutFileFormat==0) ConvertToGXML   (output, event);
    else if (gOptOutFileFormat==1) ConvertToGTab   (output, event);
    else if (gOptOutFileFormat==2) ConvertToNuance (output, event);

    mcrec->Clear();
  }

  //-- add a file footer when it is needed
  if(gOptOutFileFormat==0) AddGXMLHeader(output);

  LOG("gntpc", pINFO) << "\nDone converting GENIE's ER ntuple";
}
//___________________________________________________________________
//        *** code for converting to text-based formats ***
//___________________________________________________________________
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
  output << "$weight : " << event.Weight()   << endl;
  output << "$xsec-ev: " << event.XSec()     << endl;
  output << "$xsec-kn: " << event.DiffXSec() << endl;
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
    output << setfill(' ') << setw(12) << p->Pdg();
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

  int ghep_pdgc   = p->Pdg();
  int nuance_pdgc = ghep_pdgc;

  if ( p->IsNucleus() ) {
      int Z = pdg::IonPdgCodeToZ(ghep_pdgc);
      int A = pdg::IonPdgCodeToA(ghep_pdgc);

      nuance_pdgc = 1000*Z + A;
  }
  return nuance_pdgc;
}
//___________________________________________________________________
// ***** GENIE ER ROOT TREE -> STD NT FOR T2K CROSS-GENERATOR STUDIES ****
//___________________________________________________________________
void ConvertToT2KMCComparisonsRootFormat()
{
  //-- define branch variables
  //
  int    brIev        = 0;      // Event number (to be used as index for friend trees with detector response variables)
  int    brNeutrino   = 0;      // Neutrino pdg code
  int    brTarget     = 0;      // Nuclear target pdg code (1aaazzz000)
  int    brHitNuc     = 0;      // Hit nucleon pdg code
  int    brHitQrk     = 0;      // Hit quark pdg code
  bool   brIsQel      = false;  // Is QEL?
  bool   brIsRes      = false;  // Is RES?
  bool   brIsDis      = false;  // Is DIS?
  bool   brIsCoh      = false;  // Is COH?
  bool   brIsCC       = false;  // Is CC?
  bool   brIsNC       = false;  // Is NC?
  double brWeight     = 0;      // Event weight
  double brKineXs     = 0;      // Bjorken x (selected)
  double brKineYs     = 0;      // Inelasticity y (selected)
  double brKineTs     = 0;      // Energy transfer to nucleus at COH events (selected)
  double brKineQ2s    = 0;      // Momentum transfer Q^2 (selected)
  double brKineWs     = 0;      // Hadronic invariant mass W (selected)
  double brKineX      = 0;      // Bjorken x  (computed from the event record)
  double brKineY      = 0;      // Inelasticity y (computed from the event record)
  double brKineT      = 0;      // Energy transfer to nucleus at COH events (computed from the event record)
  double brKineQ2     = 0;      // Momentum transfer Q^2 (computed from the event record)
  double brKineW      = 0;      // Hadronic invariant mass W (computed from the event record)
  double brEv         = 0;      // Neutrino energy (neutrino assumed in +z direction)
  double brEn         = 0;      // Initial state hit nucleon energy
  double brPxn        = 0;      // Initial state hit nucleon px
  double brPyn        = 0;      // Initial state hit nucleon py
  double brPzn        = 0;      // Initial state hit nucleon pz
  double brEl         = 0;      // Final state primary lepton energy
  double brPxl        = 0;      // Final state primary lepton px
  double brPyl        = 0;      // Final state primary lepton py
  double brPzl        = 0;      // Final state primary lepton pz 
  int    brNfP        = 0;      // Number of final state p's
  int    brNfN        = 0;      // Number of final state n's
  int    brNfPip      = 0;      // Number of final state pi+'s
  int    brNfPim      = 0;      // Number of final state pi-'s
  int    brNfPi0      = 0;      // Number of final state pi0's
  int    brNfKp       = 0;      // Number of final state K+'s
  int    brNfKm       = 0;      // Number of final state K-'s
  int    brNfK0       = 0;      // Number of final state K0's
  int    brNiP        = 0;      // Number of 'primary' p's   (before intranuclear rescattering)
  int    brNiN        = 0;      // Number of 'primary' n's   (before intranuclear rescattering)
  int    brNiPip      = 0;      // Number of 'primary' pi+'s (before intranuclear rescattering)
  int    brNiPim      = 0;      // Number of 'primary' pi-'s (before intranuclear rescattering)
  int    brNiPi0      = 0;      // Number of 'primary' pi0's (before intranuclear rescattering)
  int    brNiKp       = 0;      // Number of 'primary' K+'s  (before intranuclear rescattering)
  int    brNiKm       = 0;      // Number of 'primary' K-'s  (before intranuclear rescattering)
  int    brNiK0       = 0;      // Number of 'primary' K0's  (before intranuclear rescattering) 
  int    brNf         = 0;      // Number of final state particles in hadronic system
  int    brPdgf[kNPmax];        // Pdg code of i^th final state particle in hadronic system
  double brEf  [kNPmax];        // Energy   of i^th final state particle in hadronic system
  double brPxf [kNPmax];        // Px       of i^th final state particle in hadronic system
  double brPyf [kNPmax];        // Py       of i^th final state particle in hadronic system
  double brPzf [kNPmax];        // Pz       of i^th final state particle in hadronic system
  int    brNi         = 0;      // Number of particles in 'primary' hadronic system (before intranuclear rescattering)
  int    brPdgi[kNPmax];        // Pdg code of i^th particle in 'primary' hadronic system 
  double brEi  [kNPmax];        // Energy   of i^th particle in 'primary' hadronic system 
  double brPxi [kNPmax];        // Px       of i^th particle in 'primary' hadronic system 
  double brPyi [kNPmax];        // Py       of i^th particle in 'primary' hadronic system 
  double brPzi [kNPmax];        // Pz       of i^th particle in 'primary' hadronic system 

  //-- open output file
  //
  TFile fout(gOptOutFileName.c_str(),"recreate");

  //-- create output tree
  //
  TTree * tEvtTree = new TTree("tEvtTree","event tree summary");

  //-- create tree branches
  //
  tEvtTree->Branch("iev",       &brIev,          "iev/I"       );
  tEvtTree->Branch("neu",       &brNeutrino,     "neu/I"       );
  tEvtTree->Branch("tgt" ,      &brTarget,       "tgt/I"       );
  tEvtTree->Branch("hitnuc",    &brHitNuc,       "hitnuc/I"    );
  tEvtTree->Branch("hitqrk",    &brHitQrk,       "hitqrk/I"    );
  tEvtTree->Branch("qel",       &brIsQel,        "qel/O"       );
  tEvtTree->Branch("res",       &brIsRes,        "res/O"       );
  tEvtTree->Branch("dis",       &brIsDis,        "dis/O"       );
  tEvtTree->Branch("coh",       &brIsCoh,        "coh/O"       );
  tEvtTree->Branch("cc",        &brIsCC,         "cc/O"        );
  tEvtTree->Branch("nc",        &brIsNC,         "nc/O"        );
  tEvtTree->Branch("wght",      &brWeight,       "wght/D"      );
  tEvtTree->Branch("xs",        &brKineXs,       "xs/D"        );
  tEvtTree->Branch("ys",        &brKineYs,       "ys/D"        );
  tEvtTree->Branch("ts",        &brKineTs,       "ts/D"        );
  tEvtTree->Branch("Q2s",       &brKineQ2s ,     "Q2s/D"       );
  tEvtTree->Branch("Ws",        &brKineWs,       "Ws/D"        );
  tEvtTree->Branch("x",         &brKineX,        "x/D"         );
  tEvtTree->Branch("y",         &brKineY,        "y/D"         );
  tEvtTree->Branch("t",         &brKineT,        "t/D"         );
  tEvtTree->Branch("Q2",        &brKineQ2,       "Q2/D"        );
  tEvtTree->Branch("W",         &brKineW,        "W/D"         );
  tEvtTree->Branch("Ev",        &brEv,           "Ev/D"        );
  tEvtTree->Branch("En",        &brEn,           "En/D"        );
  tEvtTree->Branch("pxn",       &brPxn,          "pxn/D"       );
  tEvtTree->Branch("pyn",       &brPyn,          "pyn/D"       );
  tEvtTree->Branch("pzn",       &brPzn,          "pzn/D"       );
  tEvtTree->Branch("El",        &brEl,           "El/D"        );
  tEvtTree->Branch("pxl",       &brPxl,          "pxl/D"       );
  tEvtTree->Branch("pyl",       &brPyl,          "pyl/D"       );
  tEvtTree->Branch("pzl",       &brPzl,          "pzl/D"       );
  tEvtTree->Branch("nfp",       &brNfP,          "nfp/I"       );
  tEvtTree->Branch("nfn",       &brNfN,          "nfn/I"       );
  tEvtTree->Branch("nfpip",     &brNfPip,        "nfpip/I"     );
  tEvtTree->Branch("nfpim",     &brNfPim,        "nfpim/I"     );
  tEvtTree->Branch("nfpi0",     &brNfPi0,        "nfpi0/I"     );
  tEvtTree->Branch("nfkp",      &brNfKp,         "nfkp/I"      );
  tEvtTree->Branch("nfkm",      &brNfKm,         "nfkm/I"      );
  tEvtTree->Branch("nfk0",      &brNfK0,         "nfk0/I"      );
  tEvtTree->Branch("nip",       &brNiP,          "np/I"        );
  tEvtTree->Branch("nin",       &brNiN,          "nn/I"        );
  tEvtTree->Branch("nipip",     &brNiPip,        "npip/I"      );
  tEvtTree->Branch("nipim",     &brNiPim,        "npim/I"      );
  tEvtTree->Branch("nipi0",     &brNiPi0,        "npi0/I"      );
  tEvtTree->Branch("niKp",      &brNiKp,         "nKp/I"       );
  tEvtTree->Branch("niKm",      &brNiKm,         "nKm/I"       );
  tEvtTree->Branch("niK0",      &brNiK0,         "nK0/I"       );
  tEvtTree->Branch("ni",        &brNi,           "ni/I"        );
  tEvtTree->Branch("pdgi",      brPdgi,          "pdgi[ni]/I " );
  tEvtTree->Branch("Ei",        brEi,            "Ei[ni]/D"    );
  tEvtTree->Branch("pxi",       brPxi,           "pxi[ni]/D"   );
  tEvtTree->Branch("pyi",       brPyi,           "pyi[ni]/D"   );
  tEvtTree->Branch("pzi",       brPzi,           "pzi[ni]/D"   );
  tEvtTree->Branch("nf",        &brNf,           "nf/I"        );
  tEvtTree->Branch("pdgf",      brPdgf,          "pdgf[nf]/I " );
  tEvtTree->Branch("Ef",        brEf,            "Ef[nf]/D"    );
  tEvtTree->Branch("pxf",       brPxf,           "pxf[nf]/D"   );
  tEvtTree->Branch("pyf",       brPyf,           "pyf[nf]/D"   );
  tEvtTree->Branch("pzf",       brPzf,           "pzf[nf]/D"   );

  //-- open the ROOT file and get the TTree & its header
  TFile fin(gOptInpFileName.c_str(),"READ");
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  tree = dynamic_cast <TTree *>           ( fin.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( fin.Get("header") );

  LOG("gntpc", pINFO) << "Input tree header: " << *thdr;

  NtpMCFormat_t format = thdr->format;
  assert(format == kNFEventRecord);

  //-- get the mc record
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  //-- event loop
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);
    EventRecord &  event = *(mcrec->event);

    LOG("gntpc", pINFO) << event;

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

    // start filling up branches
    //
    brIev        = i;      
    brNeutrino   = 0;      
    brTarget     = 0;      
    brHitNuc     = 0;      
    brHitQrk     = 0;      
    brIsQel      = false;
    brIsRes      = false;
    brIsDis      = false;  
    brIsCoh      = false;  
    brIsCC       = false;  
    brIsNC       = false;  
    brWeight     = 0.;      
    brKineXs     = 0.;      
    brKineYs     = 0.;      
    brKineTs     = 0.;      
    brKineQ2s    = 0.;            
    brKineWs     = 0.;      
    brKineX      = 0.;      
    brKineY      = 0.;      
    brKineT      = 0.;      
    brKineQ2     = 0.;      
    brKineW      = 0.;      
    brEv         = 0.;      
    brEn         = 0.;      
    brPxn        = 0.;      
    brPyn        = 0.;      
    brPzn        = 0.;            
    brEl         = 0.;      
    brPxl        = 0.;      
    brPyl        = 0.;      
    brPzl        = 0.;      
    brNfP        = 0;      
    brNfN        = 0;      
    brNfPip      = 0;      
    brNfPim      = 0;      
    brNfPi0      = 0;      
    brNfKp       = 0;      
    brNfKm       = 0;      
    brNfK0       = 0;      
    brNiP        = 0;      
    brNiN        = 0;      
    brNiPip      = 0;      
    brNiPim      = 0;      
    brNiPi0      = 0;      
    brNiKp       = 0;      
    brNiKm       = 0;      
    brNiK0       = 0;      

    brNi = 0;      
    for(int j=0; j<brNi; j++) {
      brPdgi[j] = 0;     
      brEi  [j] = 0;     
      brPxi [j] = 0;     
      brPyi [j] = 0;     
      brPzi [j] = 0;     
    }

    brNf = 0;      
    for(int j=0; j<brNf; j++) {
      brPdgf[j] = 0;     
      brEf  [j] = 0;     
      brPxf [j] = 0;     
      brPyf [j] = 0;     
      brPzf [j] = 0;     
    }

    tEvtTree->Fill();

    mcrec->Clear();
  }
  fin.Close();

  tEvtTree->Write("numcnt");
  fout.Close();
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
      gOptOutFileFormat = 3;
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
  else if (gOptOutFileFormat==2) ext = "gnuance";
  else if (gOptOutFileFormat==3) ext = "gt2k.root";

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
    << "      [] denotes an optional argument \n"
    << "      -f specifies the output text file format. Default is 3.  \n"
    << "          0 : GENIE-style XML file  \n"
    << "          1 : GENIE's GHEP in tabular text format  \n"
    << "          2 : NUANCE-style text file  \n"
    << "          3 : ROOT Tree used for T2K cross-generator comparisons  \n"
    << "  \n"
    << "    -o specifies the output filename.  \n"
    << "         If not specified a the default filename is constructed by the  \n"
    << "         input base name and an extension depending on the file format:  \n"
    << "           0 -> *.gxml  \n"
    << "           1 -> *.gtab  \n"
    << "           2 -> *.gnuance  \n"
    << "           3 -> *.gt2k.root \n"
    << ENDL;
}
//___________________________________________________________________
