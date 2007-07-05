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
		0 : GENIE XML format event file
		1 : NUANCE-style tracker text format 
	        2 : NEUGEN-style test format for hadronization model studies
	       10 : ROOT Tree used for T2K cross-generator comparisons

           -o specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
                0 -> *.gxml 
                1 -> *.gtrac
                2 -> *.ghad
               10 -> *.gt2k.root
		
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

#include "Conventions/Constants.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
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
using namespace genie::constants;

//func prototypes
void   ConvertToTextFormat (void);
void   AddHeader           (ofstream & out);
void   AddFooter           (ofstream & out);
void   ConvertToT2KMCComparisonsRootFormat();

void   ConvertToGXML       (ofstream & out, EventRecord & event);
void   ConvertToGTrac      (ofstream & out, EventRecord & event);
void   ConvertToGHad       (ofstream & out, EventRecord & event);
void   GetCommandLineArgs (int argc, char ** argv);
void   PrintSyntax        (void);
string DefaultOutputFile  (void);

//input options (from command line arguments):
string gOptInpFileName;
string gOptOutFileName;
int    gOptOutFileFormat;

// glob
int gIEv=0;

//consts
const int kNPmax = 100;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  //-- get the command line arguments
  GetCommandLineArgs(argc, argv);

  if(gOptOutFileFormat==0 || 
     gOptOutFileFormat==1 || 
     gOptOutFileFormat==2)  ConvertToTextFormat();

  else
  if(gOptOutFileFormat==10) ConvertToT2KMCComparisonsRootFormat();

  else {
    LOG("gntpc", pFATAL)
        << "Invalid output format [" << gOptOutFileFormat << "]";
    PrintSyntax();
    exit(3);
  }

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

  //-- add required header
  AddHeader(output);

  //-- event loop
  for(gIEv = 0; gIEv< tree->GetEntries(); gIEv++) {
    tree->GetEntry(gIEv);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    if      (gOptOutFileFormat==0) ConvertToGXML   (output, event);
    else if (gOptOutFileFormat==1) ConvertToGTrac  (output, event);
    else if (gOptOutFileFormat==2) ConvertToGHad   (output, event);

    mcrec->Clear();
  }

  //-- add required footer
  AddFooter(output);

  LOG("gntpc", pINFO) << "\nDone converting GENIE's ER ntuple";
}
//___________________________________________________________________
void AddHeader(ofstream & output)
{
  if(gOptOutFileFormat==0) {
     output << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
     output << endl << endl;
     output << "<!-- generated by GENIE gntpc utility -->";
     output << endl << endl;
     output << "<genie_event_list version=\"1.00\">" << endl;
  }
}
//___________________________________________________________________
void AddFooter(ofstream & output)
{
  if(gOptOutFileFormat==0) {
     output << endl << endl;
     output << "<genie_event_list version=\"1.00\">";
  }
}
//___________________________________________________________________
//        *** code for converting to text-based formats ***
//___________________________________________________________________
//      ***** GENIE ER ROOT TREE -> GENIE XML EVENT FILE ****
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
//  ***** GENIE ER ROOT TREE -> NUANCE-STYLE TRACKER TEXT FILE ****
//___________________________________________________________________
void ConvertToGTrac(ofstream & output, EventRecord & event)
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
    // For nuclei GHEP PDGC follows the MINOS-convention: 1AAAZZZ000
    // NUANCE is using: ZZZAAA
    int ghep_pdgc = p->Pdg();
    int pdgc = ghep_pdgc;
    if ( p->IsNucleus() ) {
      int Z = pdg::IonPdgCodeToZ(ghep_pdgc);
      int A = pdg::IonPdgCodeToA(ghep_pdgc);
      pdgc = 1000*Z + A;
    }

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
// *** GENIE ER ROOT TREE -> NEUGEN-style format for AGKY studies ***
//___________________________________________________________________
void ConvertToGHad(ofstream & output, EventRecord & event)
{
// Neugen-style text format for AGKY hadrinization model studies
// (blank line) 
// event number, neutrino particle code, CCNC, IM, A, Z
// int_type, x, y, w, ihadmod 
// neutrino particle code, 5 vec
// lepton particle code, 5-vec
// outgoing hadronic system, 5-vec
// number of stable daughters of hadronic system
// ... then for each stable daughter
// particle id, 5 vec 

  const Interaction * interaction = event.Summary();

  const ProcessInfo & proc_info = interaction->ProcInfo();
  //bool is_res = proc_info.IsResonant();
  bool is_dis = proc_info.IsDeepInelastic();
  if(!is_dis) return; 

  GHepParticle * neutrino = event.Probe();
  assert(neutrino);
  GHepParticle * target = event.Particle(1);
  assert(target);
  GHepParticle * fsl = event.FinalStatePrimaryLepton();
  assert(fsl);
  GHepParticle * hitnucl = event.HitNucleon();
  assert(hitnucl);
  GHepParticle * hadsyst = event.FinalStateHadronicSystem();


  const TLorentzVector & k1 = *(neutrino->P4());  // v 4-p (k1)
  const TLorentzVector & k2 = *(fsl->P4());       // l 4-p (k2)
  const TLorentzVector & p1 = *(hitnucl->P4());   // N 4-p (p1)      
  const TLorentzVector & ph = *(hadsyst->P4());   // had-syst 4-p 
     
  double M  = kNucleonMass;
  TLorentzVector q  = k1-k2;    // q=k1-k2, 4-p transfer
  double v  = (q*p1)/M;         // v (E transfer in hit nucleon rest frame)
  double Q2 = -1 * q.M2();      // momemtum transfer
  double x  = 0.5*Q2/(M*v);     // Bjorken x
  double y  = v*M/(k1*p1);      // Inelasticity, y = q*P1/k1*P1
  double W2 = M*M + 2*M*v - Q2; // Hadronic Invariant mass ^ 2
  double W  = TMath::Sqrt(W2); 

  int     nupdg  = neutrino->Pdg();
  int     fslpdg = fsl->Pdg();
  int     ccnc   = proc_info.IsWeakCC() ? 1 : 0;;
  int     im     = -1;
  int     A      = target->A();
  int     Z      = target->Z();;
  int     inttyp = 2;
  int     hadmod = 0;

  output << endl;
  output << gIEv   << "\t"  
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
  output << ph.Px()     << "\t" << ph.Py() << "\t" << ph.Pz() << "\t"
         << ph.Energy() << "\t" << ph.M()  << endl;

  int id1 = hadsyst->FirstDaughter();
  int id2 = hadsyst->LastDaughter();
  for(int i=id1; i<=id2; i++) {
    GHepParticle * p = event.Particle(i);
    GHepStatus_t ist = p->Status();
    if(ist == kIStStableFinalState) {
      int pdg = p->Pdg();
      double px = p->P4()->Px();
      double py = p->P4()->Py();
      double pz = p->P4()->Pz();
      double E  = p->P4()->Energy();
      double m  = p->P4()->M();
      output << pdg << "\t" 
             << px  << "\t" << py << "\t" << pz << "\t"
             << E   << "\t" << m  << endl;
    }
  }
}
//___________________________________________________________________
// ** GENIE ER ROOT TREE -> STD NT FOR T2K CROSS-GENERATOR STUDIES **
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
  bool   brFromSea    = false;  // Hit quark is from sea
  bool   brIsQel      = false;  // Is QEL?
  bool   brIsRes      = false;  // Is RES?
  bool   brIsDis      = false;  // Is DIS?
  bool   brIsCoh      = false;  // Is COH?
  bool   brIsImd      = false;  // Is IMD?
  bool   brIsNuEL     = false;  // Is ve elastic?
  bool   brIsCC       = false;  // Is CC?
  bool   brIsNC       = false;  // Is NC?
  bool   brIsCharmPro = false;  // Produces charm?
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

// changes:
// sea
// imd
// nuel
// charm

  //-- create tree branches
  //
  tEvtTree->Branch("iev",       &brIev,          "iev/I"       );
  tEvtTree->Branch("neu",       &brNeutrino,     "neu/I"       );
  tEvtTree->Branch("tgt" ,      &brTarget,       "tgt/I"       );
  tEvtTree->Branch("hitnuc",    &brHitNuc,       "hitnuc/I"    );
  tEvtTree->Branch("hitqrk",    &brHitQrk,       "hitqrk/I"    );
  tEvtTree->Branch("sea",       &brFromSea,      "sea/O"       );
  tEvtTree->Branch("qel",       &brIsQel,        "qel/O"       );
  tEvtTree->Branch("res",       &brIsRes,        "res/O"       );
  tEvtTree->Branch("dis",       &brIsDis,        "dis/O"       );
  tEvtTree->Branch("coh",       &brIsCoh,        "coh/O"       );
  tEvtTree->Branch("imd",       &brIsImd,        "imd/O"       );
  tEvtTree->Branch("nuel",      &brIsNuEL,       "nuel/O"      );
  tEvtTree->Branch("cc",        &brIsCC,         "cc/O"        );
  tEvtTree->Branch("nc",        &brIsNC,         "nc/O"        );
  tEvtTree->Branch("charm",     &brIsCharmPro,   "charm/O"     );
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

  TLorentzVector pdummy(0,0,0,0);

  //-- event loop
  for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);
    EventRecord &  event = *(mcrec->event);

    LOG("gntpc", pINFO) << event;

    // go further only if the event is physical
    bool is_unphysical = event.IsUnphysical();
    if(is_unphysical) {
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

    // computing event characteristics
    // This section is largely copied from gMCSampleTest.cxx which creates GENIE's
    // native 'reduced' event ntuple and analyzes the generated event sample.
    // Refer to that file for more comments

    //input particles
    GHepParticle * neutrino = event.Probe();
    assert(neutrino);
    GHepParticle * target = event.Particle(1);
    assert(target);
    GHepParticle * fsl = event.FinalStatePrimaryLepton();
    assert(fsl);
    GHepParticle * hitnucl = event.HitNucleon();
  
    //summary info
    const Interaction * interaction = event.Summary();
    const InitialState & init_state = interaction->InitState();
    const ProcessInfo &  proc_info  = interaction->ProcInfo();
    const Kinematics &   kine       = interaction->Kine();
    const XclsTag &      xcls       = interaction->ExclTag();
    const Target &       tgt        = init_state.Tgt();

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
    int  qrk  = 0;     
    bool seaq = false; 
    if(is_dis) {
      assert(tgt.HitQrkIsSet());
      qrk  = tgt.HitQrkPdg();
      seaq = tgt.HitSeaQrk(); 
    } else {
      assert(!tgt.HitQrkIsSet());
    }

    // chram production?
    bool charm = xcls.IsCharmEvent();

    //weight
    double weight = event.Weight();

    //input 4-momenta    
    const TLorentzVector & k1 = *(neutrino->P4());                     // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());                          // l 4-p (k2)
    const TLorentzVector & p1 = (hitnucl) ? *(hitnucl->P4()) : pdummy; // N 4-p (p1)      
     
    // compute kinematical params
    double M  = kNucleonMass;
    TLorentzVector q  = k1-k2;                     // q=k1-k2, 4-p transfer
    double v  = (hitnucl) ? (q*p1)/M : -1;         // v (E transfer in hit nucleon rest frame)
    double Q2 = -1 * q.M2();                       // momemtum transfer
    double x  = (hitnucl) ? 0.5*Q2/(M*v)     : -1; // Bjorken x
    double y  = (hitnucl) ? v*M/(k1*p1)      : -1; // Inelasticity, y = q*P1/k1*P1
    double W2 = (hitnucl) ? M*M + 2*M*v - Q2 : -1; // Hadronic Invariant mass ^ 2
    double W  = (hitnucl) ? TMath::Sqrt(W2)  : -1; 
    double t  = 0;

    // also, access kinematical params _exactly_ as they were selected
    // (possibly using off-shell kinematics)
    bool get_selected = true;
    double xs  = kine.x (get_selected);
    double ys  = kine.y (get_selected);
    double ts  = (is_coh) ? kine.t (get_selected) : -1;
    double Q2s = kine.Q2(get_selected);
    double Ws  = kine.W (get_selected);

    //
    TObjArrayIter piter(&event);
    GHepParticle * p = 0;
    int ip=-1;

    vector<int> final_had_syst;
    vector<int> prim_had_syst;

    // number of final state hadrons 
    int np    = 0;
    int nn    = 0;
    int npip  = 0;
    int npim  = 0;
    int npi0  = 0;
    int nKp   = 0;
    int nKm   = 0;
    int nK0   = 0;
    int ngpi0 = 0.;
    while( (p = (GHepParticle *) piter.Next()) )
    {
      ip++;
      int pdgc = p->Pdg();
      int ist  = p->Status();
      if(ist!=kIStStableFinalState) continue;
      if (pdgc == kPdgProton)  np++;
      if (pdgc == kPdgNeutron) nn++;
      if (pdgc == kPdgPiP)     npip++;
      if (pdgc == kPdgPiM)     npim++;
      if (pdgc == kPdgPi0)     npi0++;
      if (pdgc == kPdgKP)      nKp++;
      if (pdgc == kPdgKM)      nKm++;
      if (pdgc == kPdgK0)      nK0++;
      if (pdgc == kPdgGamma)  {
          int igmom = p->FirstMother();
          if(igmom!=-1) {
            if( event.Particle(igmom)->Pdg() == kPdgPi0) ngpi0++;
          }
      }
    }
    npi0 += (ngpi0/2);

    // reset iterator so as too loop again
    piter.Reset();

    // number of primary hadrons (before intranuclear rescattering)
    int np_prim    = 0;
    int nn_prim    = 0;
    int npip_prim  = 0;
    int npim_prim  = 0;
    int npi0_prim  = 0;
    int nKp_prim   = 0;
    int nKm_prim   = 0;
    int nK0_prim   = 0;
    ip=-1;
    while( (p = (GHepParticle *) piter.Next()) )
    {
      ip++;
      int pdgc = p->Pdg();
      int ist  = p->Status();
      if(ist!=kIStHadronInTheNucleus) continue;
      if (pdgc == kPdgProton)  np_prim++;
      if (pdgc == kPdgNeutron) nn_prim++;
      if (pdgc == kPdgPiP)     npip_prim++;
      if (pdgc == kPdgPiM)     npim_prim++;
      if (pdgc == kPdgPi0)     npi0_prim++;
      if (pdgc == kPdgKP)      nKp_prim++;
      if (pdgc == kPdgKM)      nKm_prim++;
      if (pdgc == kPdgK0)      nK0_prim++;
    }

    // start filling up branches
    //
    brIev        = i;      
    brNeutrino   = neutrino->Pdg();      
    brTarget     = target->Pdg();      
    brHitNuc     = (hitnucl) ? hitnucl->Pdg() : 0;      
    brHitQrk     = qrk;     
    brFromSea    = seaq;  
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
    brEn         = p1.Energy();      
    brPxn        = p1.Px();      
    brPyn        = p1.Py();      
    brPzn        = p1.Pz();            
    brEl         = k2.Energy();      
    brPxl        = k2.Px();      
    brPyl        = k2.Py();      
    brPzl        = k2.Pz();      
    brNfP        = np;      
    brNfN        = nn;      
    brNfPip      = npip;      
    brNfPim      = npim;      
    brNfPi0      = npi0;      
    brNfKp       = nKp;      
    brNfKm       = nKm;      
    brNfK0       = nK0;      
    brNiP        = np_prim;      
    brNiN        = nn_prim;      
    brNiPip      = npip_prim;      
    brNiPim      = npim_prim;      
    brNiPi0      = npi0_prim;      
    brNiKp       = nKp_prim;      
    brNiKm       = nKm_prim;      
    brNiK0       = nK0_prim;      

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
  fout.Write();
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
      gOptOutFileFormat = 10;
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
}
//___________________________________________________________________
string DefaultOutputFile(void)
{
  // filename extension - depending on file format
  string ext="";
  if      (gOptOutFileFormat==0)  ext = "gxml";
  else if (gOptOutFileFormat==1)  ext = "gtrac";
  else if (gOptOutFileFormat==2)  ext = "ghad";
  else if (gOptOutFileFormat==10) ext = "gt2k.root";

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
    << "          0 : GENIE XML format event file  \n"
    << "          1 : NUANCE-style tracker text format  \n"
    << "          2 : NEUGEN-style test format for hadronization model studies  \n"
    << "         10 : ROOT Tree used for T2K cross-generator comparisons  \n"
    << "  \n"
    << "    -o specifies the output filename.  \n"
    << "         If not specified a the default filename is constructed by the  \n"
    << "         input base name and an extension depending on the file format:  \n"
    << "           0 -> *.gxml   \n"
    << "           1 -> *.gtrac  \n"
    << "           2 -> *.ghad   \n"
    << "          10 -> *.gt2k.root \n"
    << ENDL;
}
//___________________________________________________________________
