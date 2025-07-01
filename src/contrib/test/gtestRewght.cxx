//____________________________________________________________________________
/*!

\program gtestRewght

\brief   A simple program to illustrate how to use the GENIE event reweighting.

\syntax  gtestRewght -f filename [-n nev]

         where 
         [] is an optional argument
         -f specifies a GENIE event file (GHEP format)
         -n specifies the number of events to process (default: all)

\author  Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

\created May 19, 2010

\cpright Copyright (c) 2003-2025, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________


#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Tools/ReWeight/GReWeightI.h"
#include "Tools/ReWeight/GSystSet.h"
#include "Tools/ReWeight/GReWeight.h"
#include "Tools/ReWeight/GReWeightNuXSecCCQE.h"
#include "Tools/ReWeight/GReWeightNuXSecCCQEvec.h"
#include "Tools/ReWeight/GReWeightNuXSecCCRES.h"
#include "Tools/ReWeight/GReWeightNuXSecNCRES.h"
#include "Tools/ReWeight/GReWeightNuXSecDIS.h"
#include "Tools/ReWeight/GReWeightNuXSecCOH.h"
#include "Tools/ReWeight/GReWeightNonResonanceBkg.h"
#include "Tools/ReWeight/GReWeightFGM.h"
#include "Tools/ReWeight/GReWeightDISNuclMod.h"
#include "Tools/ReWeight/GReWeightResonanceDecay.h"
#include "Tools/ReWeight/GReWeightFZone.h"
#include "Tools/ReWeight/GReWeightINuke.h"
#include "Tools/ReWeight/GReWeightAGKY.h"
#include "Framework/Utils/CmdLnArgParser.h"

using std::string;

using namespace genie;
using namespace genie::rew;

void GetCommandLineArgs (int argc, char ** argv);

int    gOptNEvt;
string gOptInpFilename;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  // open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

  LOG("test", pNOTICE) << "Input tree header: " << *thdr;

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  LOG("test", pNOTICE) << "Will process " << nev << " events";

  //
  // Create a GReWeight object and add to it a set of 
  // weight calculators
  //

  GReWeight rw;

  rw.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );
  rw.AdoptWghtCalc( "xsec_ccqe_vec",   new GReWeightNuXSecCCQEvec   );
  rw.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
  rw.AdoptWghtCalc( "xsec_ncres",      new GReWeightNuXSecNCRES     );
  rw.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
  rw.AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );
  rw.AdoptWghtCalc( "xsec_coh",        new GReWeightNuXSecCOH       );
  rw.AdoptWghtCalc( "nuclear_qe",      new GReWeightFGM             );
  rw.AdoptWghtCalc( "nuclear_dis",     new GReWeightDISNuclMod      );
  rw.AdoptWghtCalc( "hadro_res_decay", new GReWeightResonanceDecay  );
  rw.AdoptWghtCalc( "hadro_fzone",     new GReWeightFZone           );
  rw.AdoptWghtCalc( "hadro_intranuke", new GReWeightINuke           );
  rw.AdoptWghtCalc( "hadro_agky",      new GReWeightAGKY            );

  //
  // Create a list of systematic params (more to be found at GSyst.h)
  // set non-default values and re-configure.
  // Weight calculators included above must be able to handle the tweaked params.
  // Each tweaking dial t modifies a physics parameter p as:
  // p_{tweaked} = p_{default} ( 1 + t * dp/p )
  // So setting a tweaking dial to +/-1 modifies a physics quantity
  // by +/- 1sigma.
  // Default fractional errors are defined in GSystUncertainty
  // and can be overriden.
  //

  GSystSet & syst = rw.Systematics();

  syst.Set(kXSecTwkDial_NormCCQE,        +1.0);
  syst.Set(kXSecTwkDial_MaCCQEshape,     +1.0);
  syst.Set(kXSecTwkDial_NormCCRES,       -1.0);
  syst.Set(kXSecTwkDial_VecFFCCQEshape,  -1.0);
  syst.Set(kXSecTwkDial_MaCCRESshape,    -1.0);
  syst.Set(kXSecTwkDial_MvCCRESshape,    +0.5);
  syst.Set(kXSecTwkDial_NormNCRES,       +1.0);
  syst.Set(kXSecTwkDial_MaNCRESshape,    -0.7);
  syst.Set(kXSecTwkDial_MvNCRESshape,    +0.3);
  syst.Set(kXSecTwkDial_RvpCC1pi,        +0.5);
  syst.Set(kXSecTwkDial_RvnCC1pi,        +0.5);
  syst.Set(kXSecTwkDial_MaCOHpi,         -0.5);
  syst.Set(kINukeTwkDial_MFP_pi,         +1.0);
  syst.Set(kINukeTwkDial_MFP_N,          -1.0);
  syst.Set(kINukeTwkDial_FrPiProd_pi,    -0.7);
  syst.Set(kHadrAGKYTwkDial_xF1pi,       -1.0);
  syst.Set(kHadrAGKYTwkDial_pT1pi,       +1.0);
  syst.Set(kHadrNuclTwkDial_FormZone,    +1.0);
  syst.Set(kRDcyTwkDial_Theta_Delta2Npi, +1.0);

  rw.Reconfigure();

  //
  // Concrete weight calculators can be retrieved and fine-tuned.
  // For example:

  GReWeightNuXSecCCQE * rwccqe = 
    dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccqe"));
  rwccqe -> RewNue    (false); 
  rwccqe -> RewNuebar (false); 
  rwccqe -> RewNumubar(false); 

  //
  // Event loop
  //

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  for(int i = 0; i < nev; i++) {
    tree->GetEntry(i);

    EventRecord & event = *(mcrec->event);
    LOG("test", pNOTICE) << event;

    double wght = rw.CalcWeight(event);
    LOG("test", pNOTICE) << "Overall weight = " << wght;

    mcrec->Clear();
  }

  file.Close();

  LOG("test", pNOTICE)  << "Done!";
  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("test", pINFO) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {  
    LOG("test", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("test", pFATAL) 
      << "Unspecified input filename - Exiting";
    exit(1);
  }

  // number of events:
  if( parser.OptionExists('n') ) {  
    LOG("test", pINFO) << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("test", pINFO)
       << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }
}
//_________________________________________________________________________________
