//____________________________________________________________________________
/*!

\program grwght1scan

\brief   Generates weights given an input GHEP event file and for a given
         systematic parameter (supported by the ReWeight package).
         It outputs a ROOT file containing a tree with an entry for every 
         input event. Each such tree entry contains a TArrayF of all computed 
         weights and a TArrayF of all used tweak dial values. 

\syntax  grwght1scan \
           -f input_event_file 
          [-n n1[,n2]] 
           -s systematic 
           -t n_twk_diall_values
          [-p neutrino_codes] 
          [-o output_weights_file]
          [--seed random_number_seed]
          [--message-thresholds xml_file]
          [--event-record-print-level level]

         where 
         [] is an optional argument.

         -f 
            Specifies a GHEP input file.
         -n 
            Specifies an event range.
            Examples:
            - Type `-n 50,2350' to process all 2301 events from 50 up to 2350.
              Note: Both 50 and 2350 are included.
            - Type `-n 1000' to process the first 1000 events;
              from event number 0 up to event number 999.
            This is an optional argument. 
            By default GENIE will process all events.
         -t 
            Specified the number of tweak dial values between -1 and 1 
            (must be odd so as to include al -1, 0 and 1 / if it is an even
             number it will be incremented by 1)
         -s 
            Specifies the systematic param to tweak.
            See $GENIE/src/ReWeight/GSyst.h for a list of parameters and
            their corresponding label, which is what should be input here.
         -p 
            If set, grwght1scan reweights *only* the specified neutrino 
            species. The input is a comma separated list of PDG codes.
            This is an optional argument. 
            By default GENIE will reweight all neutrino species.
         -o 
            Specifies the filename of the output weight file.
            This is an optional argument. 
            By default filename is weights_<name_of_systematic_param>.root.
         --seed
            Random number seed.
         --message-thresholds
            Allows users to customize the message stream thresholds.
            The thresholds are specified using an XML file.
            See $GENIE/config/Messenger.xml for the XML schema.

\author  Jim Dobson
         Imperial College London

         Costas Andreopoulos, Jelena Ilic, Nick Grant
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 10, 2010

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>
#include <sstream>
#include <cassert>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TArrayF.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSyst.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSecNCEL.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"  
#include "ReWeight/GReWeightNuXSecDIS.h"    
#include "Utils/AppInit.h"
#include "Utils/RunOpt.h"
#include "Utils/CmdLnArgParser.h"

using std::string;
using std::ostringstream;

using namespace genie;
using namespace genie::rew;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);
void GetEventRange      (Long64_t nev_in_file, Long64_t & nfirst, Long64_t & nlast);

string      gOptInpFilename; ///< name for input file (contains input event tree)
string      gOptOutFilename; ///< name for output file (contains the output weight tree)
Long64_t    gOptNEvt1;       ///< range of events to process (1st input, if any)
Long64_t    gOptNEvt2;       ///< range of events to process (2nd input, if any)
GSyst_t     gOptSyst;        ///< input systematic param
int         gOptInpNTwk;     ///< # of tweaking dial values between [-1,1]
PDGCodeList gOptNu(false);   ///< neutrinos to consider
long int    gOptRanSeed;    ///< random number seed

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::RandGen(gOptRanSeed);
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());

  // Get the input event sample
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;
  TFile file(gOptInpFilename.c_str(),"READ");
  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );
  LOG("grwght1scan", pNOTICE) << "Input tree header: " << *thdr;
  if(!tree){
    LOG("grwght1scan", pFATAL) 
      << "Can't find a GHEP tree in input file: "<< file.GetName();
    gAbortingInErr = true;
    PrintSyntax();
    exit(1);
  }
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  Long64_t nev_in_file = tree->GetEntries();

  // The tweaking dial takes N values between [-1,1]

  const int   n_points      = gOptInpNTwk;
  const float twk_dial_min  = -1.0;
  const float twk_dial_max  =  1.0;
  const float twk_dial_step = (twk_dial_max - twk_dial_min) / (n_points-1);

  // Work-out the range of events to process
  Long64_t nfirst = 0;
  Long64_t nlast  = 0;
  GetEventRange(nev_in_file, nfirst, nlast);

  Long64_t nev = (nlast - nfirst + 1);

  //
  // Summarize
  //

  LOG("grwght1scan", pNOTICE) 
    << "\n"
    << "\n** grwght1scan: Will start processing events promptly."
    << "\nHere is a summary of inputs: "
    << "\n - Input event file: " << gOptInpFilename 
    << "\n - Processing: " << nev << " events in the range [" << nfirst << ", " << nlast << "]"
    << "\n - Systematic parameter to tweak: " << GSyst::AsString(gOptSyst)
    << "\n - Number of tweak dial values in [-1,1] : " << gOptInpNTwk
    << "\n - Neutrino species to reweight : " << gOptNu
    << "\n - Output weights to be saved in : " << gOptOutFilename 
    << "\n - Specified random number seed : " << gOptRanSeed
    << "\n\n";


  // Declare the weights and twkdial arrays 
  const int n_events = (const int) nev;
  float weights  [n_events][n_points]; 
  float twkdials [n_events][n_points]; 

  // Create a GReWeight object and add to it a set of weight calculators

  GReWeight rw;
  rw.AdoptWghtCalc( "xsec_ncel",       new GReWeightNuXSecNCEL      );
  rw.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );
  rw.AdoptWghtCalc( "xsec_ccqe_vec",   new GReWeightNuXSecCCQEvec   );
  rw.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
  rw.AdoptWghtCalc( "xsec_ncres",      new GReWeightNuXSecNCRES     );
  rw.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
  rw.AdoptWghtCalc( "xsec_coh",        new GReWeightNuXSecCOH       );
  rw.AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );
  rw.AdoptWghtCalc( "nuclear_qe",      new GReWeightFGM             );
  rw.AdoptWghtCalc( "nuclear_dis",     new GReWeightDISNuclMod      );
  rw.AdoptWghtCalc( "hadro_res_decay", new GReWeightResonanceDecay  );
  rw.AdoptWghtCalc( "hadro_fzone",     new GReWeightFZone           );
  rw.AdoptWghtCalc( "hadro_intranuke", new GReWeightINuke           );
  rw.AdoptWghtCalc( "hadro_agky",      new GReWeightAGKY            );

  // Get GSystSet and include the (single) input systematic parameter

  GSystSet & syst = rw.Systematics();
  syst.Init(gOptSyst);

  // Fine-tune weight calculators

  if(gOptSyst == kXSecTwkDial_MaCCQE) {
     // By default GReWeightNuXSecCCQE is in `NormAndMaShape' mode 
     // where Ma affects the shape of dsigma/dQ2 and a different param affects the normalization
     // If the input is MaCCQE, switch the weight calculator to `Ma' mode
     GReWeightNuXSecCCQE * rwccqe =      
        dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccqe"));  
     rwccqe->SetMode(GReWeightNuXSecCCQE::kModeMa);
  }
  if(gOptSyst == kXSecTwkDial_MaCCRES || gOptSyst == kXSecTwkDial_MvCCRES) {
     // As above, but for the GReWeightNuXSecCCRES weight calculator
     GReWeightNuXSecCCRES * rwccres = 
        dynamic_cast<GReWeightNuXSecCCRES *> (rw.WghtCalc("xsec_ccres"));  
     rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);
  }
  if(gOptSyst == kXSecTwkDial_MaNCRES || gOptSyst == kXSecTwkDial_MvNCRES) {
     // As above, but for the GReWeightNuXSecNCRES weight calculator
     GReWeightNuXSecNCRES * rwncres = 
        dynamic_cast<GReWeightNuXSecNCRES *> (rw.WghtCalc("xsec_ncres"));  
     rwncres->SetMode(GReWeightNuXSecNCRES::kModeMaMv);
  }
  if(gOptSyst == kXSecTwkDial_AhtBYshape  || gOptSyst == kXSecTwkDial_BhtBYshape ||
     gOptSyst == kXSecTwkDial_CV1uBYshape || gOptSyst == kXSecTwkDial_CV2uBYshape ) {
     // Similarly for the GReWeightNuXSecDIS weight calculator.
     // There the default behaviour is for the Aht, Bht, CV1u and CV2u Bodek-Yang
     // params to affects both normalization and dsigma/dxdy shape.
     // Switch mode if a shape-only param is specified.
     GReWeightNuXSecDIS * rwdis = 
        dynamic_cast<GReWeightNuXSecDIS *> (rw.WghtCalc("xsec_dis"));  
     rwdis->SetMode(GReWeightNuXSecDIS::kModeABCV12uShape);
  }

  // Twk dial loop
  for(int ith_dial = 0; ith_dial < n_points; ith_dial++){  

     // Set non-default values and re-configure.    
     double twk_dial = twk_dial_min + ith_dial * twk_dial_step;  
     LOG("grwght1scan", pNOTICE) 
       << "\n\nReconfiguring systematic: " << GSyst::AsString(gOptSyst)
       << " - Setting tweaking dial to: " << twk_dial;
     syst.Set(gOptSyst, twk_dial);
     rw.Reconfigure();

     // Event loop
     for(int iev = nfirst; iev <= nlast; iev++) {

          if(iev%100 == 0) {
              LOG("grwght1scan", pNOTICE) 
                 << "***** Currently at event number: "<< iev;
          }

          // Get next event
          tree->GetEntry(iev);
          EventRecord & event = *(mcrec->event);
          LOG("grwght1scan", pINFO) << "Event: " << iev << "\n" << event;

          // Reset arrays
          int idx = iev - nfirst;
          weights  [idx][ith_dial] = -99999.0;
          twkdials [idx][ith_dial] = twk_dial;

          // Reweight this event?
          int nupdg = event.Probe()->Pdg();
          bool do_reweight = gOptNu.ExistsInPDGCodeList(nupdg);

          // Calculate weight
          double wght=1.;
          if(do_reweight) {
	     wght = rw.CalcWeight(event);
          }

          // Print/store
          LOG("grwght1scan", pDEBUG) 
              << "Overall weight = " << wght;
          weights[idx][ith_dial] = wght;

          // Clean-up
          mcrec->Clear();

      } // evt loop
  } // twk_dial loop

  // Close event file
  file.Close();

  //
  // Save weights 
  //

  // Make an output tree for saving the weights. As only considering 
  // varying a single systematic use this for name of tree.
  TFile * wght_file = new TFile(gOptOutFilename.c_str(), "RECREATE");
  TTree * wght_tree = new TTree(GSyst::AsString(gOptSyst).c_str(), "GENIE weights tree");
  int branch_eventnum = 0;
  TArrayF * branch_weight_array   = new TArrayF(n_points);
  TArrayF * branch_twkdials_array = new TArrayF(n_points);  
  wght_tree->Branch("eventnum", &branch_eventnum);
  wght_tree->Branch("weights",  &branch_weight_array);
  wght_tree->Branch("twkdials", &branch_twkdials_array);

  for(int iev = nfirst; iev <= nlast; iev++) {
    int idx = iev - nfirst;
    branch_eventnum = iev; 
    for(int ith_dial = 0; ith_dial < n_points; ith_dial++){  
        LOG("grwght1scan", pDEBUG)
          << "Filling tree with wght = " << weights[idx][ith_dial] 
          << ", twk dial = "<< twkdials[idx][ith_dial];
       branch_weight_array   -> AddAt (weights [idx][ith_dial], ith_dial);
       branch_twkdials_array -> AddAt (twkdials[idx][ith_dial], ith_dial);
    } // twk_dial loop
    wght_tree->Fill();
  } 

  wght_file->cd();
  wght_tree->Write();
  delete wght_tree; 
  wght_tree = 0; 
  wght_file->Close();

  LOG("grwght1scan", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("grwght1scan", pINFO) 
     << "*** Parsing command line arguments";

  LOG("grwght1scan", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if(parser.OptionExists('f')) {
    LOG("grwght1scan", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("grwght1scan", pFATAL) 
        << "Unspecified input filename - Exiting";
    gAbortingInErr = true;
    PrintSyntax();
    exit(1);
  }

  // range of event numbers to process
  if ( parser.OptionExists('n') ) {
    //
    LOG("grwght1scan", pINFO) << "Reading number of events to analyze";
    string nev =  parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n',",");
      if(vecn.size()!=2) {
         LOG("grwght1scan", pFATAL) << "Invalid syntax";
         gAbortingInErr = true;
         PrintSyntax();
         exit(1);
      }
      // User specified a comma-separated set of values n1,n2.
      // Use [n1,n2] as the event range to process.
      gOptNEvt1 = vecn[0];
      gOptNEvt2 = vecn[1];
    } else {
      // User specified a single number n.
      // Use [0,n] as the event range to process.
      gOptNEvt1 = -1;
      gOptNEvt2 = parser.ArgAsLong('n');
    }  
  } else {
    LOG("grwght1scan", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt1 = -1;
    gOptNEvt2 = -1;
  }
  LOG("grwght1scan", pDEBUG) 
    << "Input event range: " << gOptNEvt1 << ", " << gOptNEvt2;

  // get the number of tweak dials to scan
  if(parser.OptionExists('t')) {
    LOG("grwght1scan", pINFO) 
       << "Reading number of tweak dial values";
    gOptInpNTwk = parser.ArgAsInt('t');
    if(gOptInpNTwk % 2 == 0)
    { 
      gOptInpNTwk+=1;
    }
    if(gOptInpNTwk < 3)
    { 
      LOG("grwght1scan", pFATAL)
	 << "Specified number of tweak dial is too low, min value is 3 - Exiting";
      gAbortingInErr = true;
      PrintSyntax();
      exit(1);
    } 
  } else {
     LOG("grwght1scan", pFATAL) 
       << "Unspecified number of tweak dials - Exiting";
     gAbortingInErr = true;
     PrintSyntax();
     exit(1);
  }

  // get the systematics
  if(parser.OptionExists('s')) {
   LOG("grwght1scan", pINFO) 
      << "Reading input systematic parameter";
   string systematic = parser.ArgAsString('s');
   gOptSyst = GSyst::FromString(systematic);
   if(gOptSyst == kNullSystematic) {
      LOG("grwght1scan", pFATAL) << "Unknown systematic: " << systematic;
      gAbortingInErr = true;
      PrintSyntax();
      exit(1);
   }
  } else {
    LOG("grwght1scan", pFATAL) 
       << "You need to specify a systematic param using -s";
    gAbortingInErr = true;
    PrintSyntax();
    exit(1);
  }

  // output weight file
  if(parser.OptionExists('o')) {
    LOG("grwght1scan", pINFO) << "Reading requested output filename";
    gOptOutFilename = parser.ArgAsString('o');
  } else {
    LOG("grwght1scan", pINFO) << "Setting default output filename";
    ostringstream nm;
    nm << "weights_" << GSyst::AsString(gOptSyst) << ".root";
    gOptOutFilename = nm.str();
  }

  // which species to reweight?
  if(parser.OptionExists('p')) {
   LOG("grwght1scan", pINFO) 
      << "Reading input list of neutrino codes";
   vector<int> vecpdg = parser.ArgAsIntTokens('p',",");
   if(vecpdg.size()==0) {
      LOG("grwght1scan", pFATAL) 
         << "Empty list of neutrino codes!?";
      gAbortingInErr = true;
      PrintSyntax();
      exit(1);
   }
   vector<int>::const_iterator it = vecpdg.begin();
   for( ; it!=vecpdg.end(); ++it) {
     gOptNu.push_back(*it);
   }
  } else {
    LOG("grwght1scan", pINFO) 
       << "Considering all neutrino species";
    gOptNu.push_back (kPdgNuE      );
    gOptNu.push_back (kPdgAntiNuE  );
    gOptNu.push_back (kPdgNuMu     );
    gOptNu.push_back (kPdgAntiNuMu );
    gOptNu.push_back (kPdgNuTau    );
    gOptNu.push_back (kPdgAntiNuTau);
  }

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("grwght1scan", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("grwght1scan", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

}
//_________________________________________________________________________________
void GetEventRange(Long64_t nev_in_file, Long64_t & nfirst, Long64_t & nlast)
{
  nfirst = 0;
  nlast  = 0;

  if(gOptNEvt1>=0 && gOptNEvt2>=0) {
    // Input was `-n N1,N2'.
    // Process events [N1,N2].
    // Note: Incuding N1 and N2.
    nfirst = gOptNEvt1;
    nlast  = TMath::Min(nev_in_file-1, gOptNEvt2);
  }
  else 
  if(gOptNEvt1<0 && gOptNEvt2>=0) {
    // Input was `-n N'.
    // Process first N events [0,N). 
    // Note: Event N is not included.
    nfirst = 0;
    nlast  = TMath::Min(nev_in_file-1, gOptNEvt2-1);
  }
  else
  if(gOptNEvt1<0 && gOptNEvt2<0) {
    // No input. Process all events.
    nfirst = 0;
    nlast  = nev_in_file-1;    
  }

  assert(nfirst < nlast && nfirst >= 0 && nlast <= nev_in_file-1);
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("grwght1scan", pFATAL)
     << "\n\n"
     << "grwght1scan                  \n"
     << "     -f input_event_file     \n"
     << "    [-n n1[,n2]]             \n"
     << "     -s systematic           \n"
     << "     -t n_twk_diall_values   \n"
     << "    [-p neutrino_codes]      \n"
     << "    [-o output_weights_file] \n"
     << "    [--seed random_number_seed] \n"
     << "    [--message-thresholds xml_file]\n"
     << "    [--event-record-print-level level]\n\n\n"
     << " See the GENIE Physics and User manual for more details";      
}
//_________________________________________________________________________________

