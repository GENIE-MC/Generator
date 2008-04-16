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
                3 -> *.gtrac.root
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
#include <TBits.h>
#include <TObjString.h>

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
#include "PDG/PDGLibrary.h"
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
void   ConvertToGT2KTracker    (void);
void   ConvertToGT2KRooTracker (void);
void   ConvertToGXML           (void);
void   ConvertToGHad           (void);
void   GetCommandLineArgs      (int argc, char ** argv);
void   PrintSyntax             (void);
string DefaultOutputFile       (void);

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
	ConvertToGT2KTracker();        
	break;
   case (3) :  
	ConvertToGT2KRooTracker(); 
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
void ConvertToGT2KTracker(void)
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

  //-- event loop
  for(gIEv = 0; gIEv< tree->GetEntries(); gIEv++) {
    tree->GetEntry(gIEv);
    NtpMCRecHeader rec_header = mcrec->hdr;
    EventRecord &  event      = *(mcrec->event);
    Interaction * interaction = event.Summary();

    LOG("gntpc", pINFO) << rec_header;
    LOG("gntpc", pINFO) << event;

    GHepParticle * p = 0;
    TIter event_iter(&event);

    //
    // convert the current event
    //

    // add tracker begin tag
    output << "$ begin" << endl;

    // add event type
    if(gOptOutFileFormat==1) {
    	// nuance event type
    	int evtype = 0;
        /*
        conversion not tested yet
        const ProcessInfo &  proc = interaction->ProcInfo();
        const InitialState & init = interaction->InitState();
        if      (proc.IsQuasiElastic()   && proc.IsWeakCC()) evtype =  1;
        else if (proc.IsQuasiElastic()   && proc.IsWeakNC()) evtype =  2;
        else if (proc.IsDeepInelastic()  && proc.IsWeakCC()) evtype = 91;        
        else if (proc.IsDeepInelastic()  && proc.IsWeakNC()) evtype = 92;        
        else if (proc.IsCoherentPiProd() && proc.IsWeakNC()) evtype = 96;        
        else if (proc.IsCoherentPiProd() && proc.IsWeakCC()) evtype = 97;        
        else if (proc.IsNuElectronElastic())                 evtype = 98;
        else if (proc.IsInverseMuDecay())                    evtype = 99;
        else if (proc.IsResonant()) {
           int nn=0, np=0, npi0=0, npip=0, npim=0;
           bool nuclear_target = init.Tgt().IsNucleus();
           GHepStatus_t matched_ist = (nuclear_target) ? 
                     kIStHadronInTheNucleus : kIStStableFinalState;
           while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) 
           {
               GHepStatus_t ghep_ist = (GHepStatus_t) p->Status();  
               if(ghep_ist != matched_ist) continue;
               int ghep_pdgc = p->Pdg();
               if(ghep_pdgc == kPdgProton ) np++;
               if(ghep_pdgc == kPdgNeutron) nn++;
               if(ghep_pdgc == kPdgPi0)     npi0++;
               if(ghep_pdgc == kPdgPiP)     npip++;
               if(ghep_pdgc == kPdgPiM)     npim++;
           }
           if(proc.IsWeakCC() && init.IsNuP()) {
             // v p -> l- p pi+ 
             if(np==1 && nn==0 && npip==1 && npi0==0 && npim==0) evtype = 3;  
           }
           if(proc.IsWeakCC() && init.IsNuN()) {
             // v n -> l- p pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 4;  
             // v n -> l- n pi+ 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 5;  
           }
           if(proc.IsWeakNC() && init.IsNuP()) {
             // v p -> v p pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 6;  
             // v p -> v n pi+ 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 7;  
           }
           if(proc.IsWeakNC() && init.IsNuN()) {
             // v n -> v n pi0 
             if(np==0 && nn==1 && npip==0 && npi0==1 && npim==0) evtype = 8;  
             // v n -> v p pi- 
             if(np==1 && nn==0 && npip==0 && npi0==0 && npim==1) evtype = 9;  
           }
           if(proc.IsWeakCC() && init.IsNuBarN()) {
             // \bar{v} n -> l+ n pi- 
             if(np==1 && nn==0 && npip==1 && npi0==0 && npim==0) evtype = 10; 
           }
           if(proc.IsWeakCC() && init.IsNuBarP()) {
             // \bar{v} p -> l+ n pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 11; 
             // \bar{v} p -> l+ p pi- 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 12; 
           }
           if(proc.IsWeakNC() && init.IsNuBarP()) {
             // \bar{v} p -> \bar{v} p pi0 
             if(np==1 && nn==0 && npip==0 && npi0==1 && npim==0) evtype = 13;
             // \bar{v} p -> \bar{v} n pi+ 
             if(np==0 && nn==1 && npip==1 && npi0==0 && npim==0) evtype = 14; 
           }
           if(proc.IsWeakNC() && init.IsNuBarN()) {
             // \bar{v} n -> \bar{v} n pi0 
             if(np==0 && nn==1 && npip==0 && npi0==1 && npim==0) evtype = 15;
             // \bar{v} n -> \bar{v} p pi-  
             if(np==1 && nn==0 && npip==0 && npi0==0 && npim==1) evtype = 16;
           }
        }
        */
    	output << "$ nuance " << evtype << endl;
    } else {
    	// genie "event type"
    	output << "$ genie " << interaction->AsString() << endl;
    }

    // add tracker vertex info
    double vtxx = 0, vtxy = 0, vtxz = 0, vtxt = 0;
    output << "$ vertex " << vtxx << " " << vtxy
           << " " << vtxz << " " << vtxt << " " << endl;

    // add 'tracks' (GENIE's equivalent of GHepParticles)
    bool info_added  = false;
    event_iter.Reset();
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) 
    {
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
       output << "$ track " << pdgc << " " << E << " "
              << dcosx << " " << dcosy << " " << dcosz << " "
              << ist << endl;
    }
    //add  tracker end tag
    output << "$ end" << endl;

    mcrec->Clear();
  } // event loop

  // add tracker end-of-file tag
  output << "$ stop" << endl;

  output.close();
  fin.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";
}
//___________________________________________________________________
//  *** GENIE GHEP EVENT TREE -> A T2K ROOT TRACKER TYPE FORMAT ***
//___________________________________________________________________
void ConvertToGT2KRooTracker(void)
{
#ifdef __GENIE_FLUX_DRIVERS_ENABLED__

  //-- get pdglib
  PDGLibrary * pdglib = PDGLibrary::Instance();

  //-- define the output tree branches
  TBits*      brEvtFlags = 0;             // generator-specific event flags
  TObjString* brEvtCode = 0;              // generator-specific string with 'event code'
  int         brEvtNum;                   // event num.
  double      brEvtXSec;                  // cross section for selected event (1E-38 cm2)
  double      brEvtDXSec;                 // cross section for selected event kinematics (1E-38 cm2 /{K^n})
  double      brEvtWght;                  // weight for that event
  double      brEvtProb;                  // probability for that event (given cross section, path lengths, etc)
  double      brEvtVtx[4];                // event vertex position in detector coord syst (in geom units)
  int         brStdHepN;                  // number of particles in particle array 
  // > stdhep-like particle array:
  int         brStdHepPdg   [kNPmax];     // pdg codes (& generator specific codes for pseudoparticles)
  int         brStdHepStatus[kNPmax];     // generator-specific status code
  double      brStdHepX4    [kNPmax][4];  // 4-x (x, y, z, t) of particle in hit nucleus frame (fm)
  double      brStdHepP4    [kNPmax][4];  // 4-p (px,py,pz,E) of particle in LAB frame (GeV)
  double      brStdHepPolz  [kNPmax][3];  // polarization vector
  int         brStdHepFd    [kNPmax];     // first daughter
  int         brStdHepLd    [kNPmax];     // last  daughter 
  int         brStdHepFm    [kNPmax];     // first mother
  int         brStdHepLm    [kNPmax];     // last  mother
  // > neutrino parent info:
  int         brNuParentPdg;              // parent hadron pdg code
  int         brNuParentDecMode;          // parent hadron decay mode
  double      brNuParentDecP4 [4];        // parent hadron 4-momentum at decay
  double      brNuParentDecX4 [4];        // parent hadron 4-position at decay
  double      brNuParentProP4 [4];        // parent hadron 4-momentum at production
  double      brNuParentProX4 [4];        // parent hadron 4-position at production
  int         brNuParentProNVtx;          // parent hadron vtx id

  //-- open the output ROOT file
  TFile fout(gOptOutFileName.c_str(), "RECREATE");

  //-- create the output ROOT tree
  TTree * rootracker_tree = new TTree("gRooTracker","");

  //-- create the output ROOT tree branches
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
  rootracker_tree->Branch("StdHepX4",         brStdHepX4,        "StdHepX4[StdHepN][4]/D"); 
  rootracker_tree->Branch("StdHepP4",         brStdHepP4,        "StdHepP4[StdHepN][4]/D"); 
  rootracker_tree->Branch("StdHepPolz",       brStdHepPolz,      "StdHepPolz[StdHepN][3]/D"); 
  rootracker_tree->Branch("StdHepFd",         brStdHepFd,        "StdHepFd[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepLd",         brStdHepLd,        "StdHepLd[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepFm",         brStdHepFm,        "StdHepFm[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepLm",         brStdHepLm,        "StdHepLm[StdHepN]/I"); 
  rootracker_tree->Branch("NuParentPdg",     &brNuParentPdg,     "NuParentPdg/I");       
  rootracker_tree->Branch("NuParentDecMode", &brNuParentDecMode, "NuParentDecMode/I");   
  rootracker_tree->Branch("NuParentDecP4",    brNuParentDecP4,   "NuParentDecP4[4]/D");     
  rootracker_tree->Branch("NuParentDecX4",    brNuParentDecX4,   "NuParentDecX4[4]/D");     
  rootracker_tree->Branch("NuParentProP4",    brNuParentProP4,   "NuParentProP4[4]/D");     
  rootracker_tree->Branch("NuParentProX4",    brNuParentProX4,   "NuParentProX4[4]/D");     
  rootracker_tree->Branch("NuParentProNVtx", &brNuParentProNVtx, "NuParentProNVtx/I");   

  //-- open the input GENIE ROOT file and get the TTree & its header
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
    LOG("gntpc", pINFO) << *flux_info;

    //
    // clear output tree branches
    //

    brEvtFlags  = 0;
    brEvtCode   = 0;
    brEvtNum    = 0;    
    brEvtXSec   = 0;
    brEvtDXSec  = 0;
    brEvtWght   = 0;
    brEvtProb   = 0;
    for(int k=0; k<4; k++) { 
      brEvtVtx[k] = 0;
    }
    brStdHepN = event.GetEntries(); 
    for(int i=0; i<kNPmax; i++) {
       brStdHepPdg   [i] = 0;  
       brStdHepStatus[i] = 0;  
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

    //
    // copy current event info to output tree
    //

    brEvtFlags  = new TBits(*event.EventFlags());   
    brEvtCode   = new TObjString(event.Summary()->AsString().c_str());   
    brEvtNum    = gIEv;    
    brEvtXSec   = (1E+38/units::cm2) * event.XSec();    
    brEvtDXSec  = (1E+38/units::cm2) * event.DiffXSec();    
    brEvtWght   = event.Weight();    
    brEvtProb   = event.Probability();    
    brEvtVtx[0] = event.Vertex()->X();    
    brEvtVtx[1] = event.Vertex()->Y();    
    brEvtVtx[2] = event.Vertex()->Z();    
    brEvtVtx[3] = event.Vertex()->T();    

    brStdHepN = event.GetEntries(); 

    int iparticle=0;
    GHepParticle * p = 0;
    TIter event_iter(&event);
    while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
        assert(p);
        brStdHepPdg   [iparticle] = p->Pdg(); 
        brStdHepStatus[iparticle] = (int) p->Status(); 
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

    brNuParentPdg       = flux_info->pdg;        
    brNuParentDecMode   = flux_info->decayMode;        

    brNuParentDecP4 [0] = flux_info->decayP * flux_info->decayDirX; // px
    brNuParentDecP4 [1] = flux_info->decayP * flux_info->decayDirY; // py
    brNuParentDecP4 [2] = flux_info->decayP * flux_info->decayDirZ; // px
    brNuParentDecP4 [3] = TMath::Sqrt(
                              TMath::Power(pdglib->Find(flux_info->pdg)->Mass(), 2.)
                            + TMath::Power(flux_info->decayP, 2.)
                           ); // E
    brNuParentDecX4 [0] = flux_info->decayX; // x
    brNuParentDecX4 [1] = flux_info->decayY; // y       
    brNuParentDecX4 [2] = flux_info->decayZ; // x   
    brNuParentDecX4 [3] = 0; // t

    brNuParentProP4 [0] = flux_info->prodP * flux_info->prodDirX; // px
    brNuParentProP4 [1] = flux_info->prodP * flux_info->prodDirY; // py
    brNuParentProP4 [2] = flux_info->prodP * flux_info->prodDirZ; // px
    brNuParentProP4 [3] = TMath::Sqrt(
                              TMath::Power(pdglib->Find(flux_info->pdg)->Mass(), 2.)
                            + TMath::Power(flux_info->prodP, 2.)
                           ); // E
    brNuParentProX4 [0] = flux_info->prodX; // x
    brNuParentProX4 [1] = flux_info->prodY; // y       
    brNuParentProX4 [2] = flux_info->prodZ; // x   
    brNuParentProX4 [3] = 0; // t

    brNuParentProNVtx   = flux_info->prodNVtx;

    rootracker_tree->Fill();
    mcrec->Clear();

  } // event loop

  fin.Close();

  fout.Write();
  fout.Close();

  LOG("gntpc", pINFO) << "\nDone converting GENIE's GHEP ntuple";

#else

  LOG("gntpc", pWARN) 
    << "\n You should enable --with-flux-drivers during the GENIE"
    << " installation step so as to access the JPARC neutrino flux"
    << " pass-through info for the simulated neutrino interactions.";

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
  else if (gOptOutFileFormat==3)    ext = "gtrac.root";
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
  string basedir  = string( gSystem->Getenv("GENIE") );
  string thisfile = basedir + string("/src/stdapp/gNtpConv.cxx");
  string cmd      = "more " + thisfile;

  gSystem->Exec(cmd.c_str());
}
//___________________________________________________________________
