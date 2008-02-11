//____________________________________________________________________________
/*!

\program gNtpConv

\brief   Converts a GENIE neutrino events (from GHEP GENIE ROOT Trees) to a 
         variety of textual formats (typically used in legacy systems), in XML
         format or in summary ROOT ntuples.

         Syntax:
              gntpc -i input_filename [-o output_filename] -f fmt

         Options :
           [] denotes an optional argument
           -f specifies the output file format. 
	      < T2K/GENIE formats >
   		    1 : NUANCE-style tracker text-based format 
		    2 : A slightly tweaked NUANCE-style tracker text-based 
                        format - A fast & dirty way for getting GENIE event 
                        samples into the T2K (nd280/2km/SuperK) Monte Carlo.
   	            3 : A standardized bare-ROOT event-tree for getting GENIE
                        neutrino  & pass-through JPARC flux info into the
                        T2K (nd280/2km/SuperK) Monte Carlo
	      < Generic GENIE XML / tabular or bare-ROOT formats >
                  100 : GENIE XML format 
	      < GENIE test / cross-generator comparisons >
	          901 : NEUGEN-style text-based format for hadronization 
                        model studies
           -o specifies the output filename. 
              If not specified a the default filename is constructed by the 
              input base name and an extension depending on the file format: 
                1 -> *.gtrac
                2 -> *.gtrac2
                3 -> *.gt2k.root
              100 -> *.gxml 
              901 -> *.ghad
		
\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created September 23, 2005

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
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

#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
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
void   ConvertToGTrac       (void);
void   ConvertToT2KRootTree (void);
void   ConvertToGXML        (void);
void   ConvertToGHad        (void);
void   GetCommandLineArgs   (int argc, char ** argv);
void   PrintSyntax          (void);
string DefaultOutputFile    (void);

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

  //-- call the appropriate conversion function
  switch(gOptOutFileFormat) {
   case (1)  :  
   case (2)  :  
	ConvertToGTrac();        
	break;
   case (3) :  
	ConvertToT2KRootTree(); 
	break;
   case (100) :  
	ConvertToGXML();         
	break;
   case (901) :  
	ConvertToGHad();         
	break;
   default:
     LOG("gntpc", pFATAL)
          << "Invalid output format [" << gOptOutFileFormat << "]";
     PrintSyntax();
     exit(3);
  }
  return 0;
}
//___________________________________________________________________
//    **** GENIE GHEP EVENT TREE -> GENIE XML EVENT FILE ****
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

  //-- event loop
  for(gIEv = 0; gIEv< tree->GetEntries(); gIEv++) {
    tree->GetEntry(gIEv);
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
// **** GENIE GHEP EVENT TREE -> NUANCE-STYLE TRACKER TEXT FILE ****
//___________________________________________________________________
void ConvertToGTrac(void)
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

//#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
//  flux::GJPARCNuFluxPassThroughInfo * flux_info = 0;
//  tree->SetBranchAddress("flux", &flux_info);
//#endif

  //-- open the output stream
  ofstream output(gOptOutFileName.c_str(), ios::out);

  //-- event loop
  for(gIEv = 0; gIEv< tree->GetEntries(); gIEv++) {
    tree->GetEntry(gIEv);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);
    Interaction * interaction = event.Summary();

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    //
    // convert the current event
    //

    TIter event_iter(&event);

    // Nuance begin tag
    output << "$ begin" << endl;

    if(gOptOutFileFormat==1) {
    	// Nuance event type
    	int nuance_event_type = 0;
    	output << "$ nuance " << nuance_event_type << endl;
    } else {
    	output << "$ genie " << interaction->AsString() << endl;
    }

    // Nuance vertex info
    double vtxx = 0, vtxy = 0, vtxz = 0, vtxt = 0;
    output << "$ vertex " << vtxx << " " << vtxy
           << " " << vtxz << " " << vtxt << " " << endl;

    // Nuance 'tracks', GENIE's equivalent of GHepParticle
    GHepParticle * p = 0;
    bool info_added  = false;

    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

       // Neglect all GENIE pseudo-particles
       //
       if(p->IsFake()) continue;

       // Convert GENIE's GHEP pdgc & status to NUANCE's equivalent
       //
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
       // For nuclei GENIE follows the PDG-convention: 10LZZZAAAI
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
       double E  = p4->Energy() / units::MeV;
       double Px = p4->Px()     / units::MeV;
       double Py = p4->Py()     / units::MeV;
       double Pz = p4->Pz()     / units::MeV;
       double P  = p4->P()      / units::MeV;
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

    mcrec->Clear();
  } // event loop

  // Nuance end-of-file tag
  output << "$ stop" << endl;

  output.close();
  fin.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//___________________________________________________________________
//    *** GENIE GHEP EVENT TREE -> T2K BARE-ROOT EVENT TREE ***
//___________________________________________________________________
void ConvertToT2KRootTree(void)
{
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__

  //-- open the output stream
  TFile fout(gOptOutFileName.c_str(), "RECREATE");

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

  //-- event loop
  for(gIEv = 0; gIEv< tree->GetEntries(); gIEv++) {
    tree->GetEntry(gIEv);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);
    Interaction * interaction = event.Summary();

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;
    LOG("gntpc", pINFO) << *interaction;

    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {

    }

    mcrec->Clear();

  } // event loop

  fin.Close();
  fout.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";

#else

  LOG("gntpc", pWARN) 
    << "\n You should enable --with-flux-drivers during the GENIE"
    << " installation step so as to access the JPARC neutrino flux"
    << " pass-through info for the simulated neutrino intarections.";

#endif
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

  //-- event loop
  for(gIEv = 0; gIEv< tree->GetEntries(); gIEv++) {
    tree->GetEntry(gIEv);
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
    if(!is_dis) return; 

    int ccnc   = proc_info.IsWeakCC() ? 1 : 0;
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
//  GHepParticle * hitnucl = event.HitNucleon();
//  assert(hitnucl);
    GHepParticle * hadsyst = event.FinalStateHadronicSystem();

    int nupdg  = neutrino->Pdg();
    int fslpdg = fsl->Pdg();
    int A      = target->A();
    int Z      = target->Z();

    const TLorentzVector & k1 = *(neutrino->P4());  // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());       // l 4-p (k2)
//  const TLorentzVector & p1 = *(hitnucl->P4());   // N 4-p (p1)      
    const TLorentzVector & ph = *(hadsyst->P4());   // had-syst 4-p 
     
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
    brIev = gIEv;   
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
//            FUNCTIONS FOR PARSING CMD-LINE ARGUMENTS 
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
  if      (gOptOutFileFormat==1)    ext = "gtrac";
  else if (gOptOutFileFormat==2)    ext = "gtrac2";
  else if (gOptOutFileFormat==3)    ext = "gt2k.root";
  else if (gOptOutFileFormat==100)  ext = "gxml";
  else if (gOptOutFileFormat==901)  ext = "ghad";

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
  LOG("gntpc", pNOTICE)
    << "\n\n" 
    << "Syntax: \n"
    << "  gntpc -i input_filename [-o output_filename] [-f fmt] \n"
    << "\n"
    << "\n"
    << "  Options : \n"
    << "       [] denotes an optional argument \n"
    << "       -f specifies the output file format. \n"
    << "          < T2K/GENIE formats > \n"
    << "                1 : NUANCE-style tracker text-based format \n"
    << "                2 : A slightly tweaked NUANCE-style tracker text-based \n"
    << "                    format - A fast & dirty way for getting GENIE event \n"
    << "                    samples into the T2K (nd280/2km/SuperK) Monte Carlo. \n"
    << "                3 : A standardized bare-ROOT event-tree for getting GENIE \n"
    << "                    neutrino  & pass-through JPARC flux info into the \n"
    << "                    T2K (nd280/2km/SuperK) Monte Carlo \n"
    << "          < Generic GENIE XML / tabular or bare-ROOT formats > \n"
    << "              100 : GENIE XML format \n"
    << "          < GENIE test / cross-generator comparisons > \n"
    << "              901 : NEUGEN-style text-based format for hadronization \n"
    << "                    model studies \n"
    << "      -o specifies the output filename. \n"
    << "                 If not specified a the default filename is constructed by the \n"
    << "                 input base name and an extension depending on the file format: \n"
    << "                   1 -> *.gtrac \n"
    << "                   2 -> *.gtrac2 \n"
    << "                   3 -> *.gt2k.root \n"
    << "                 100 -> *.gxml \n"
    << "                 901 -> *.ghad \n";
}
//___________________________________________________________________
