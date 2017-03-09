//____________________________________________________________________________
/*!

\program gRwghtZExpDirect

\brief   A simple program to reweight the GENIE z-expansion axial form factor
         from one set of z-expansion parameters directly to another

\syntax  grwghtzexpaxff -f filename -v val1,val2,val3,val4 [-n nev] [-o fileOutName] [-m norm]

         where 
         [] is an optional argument
         -f specifies a GENIE event file (GHEP format)
         -o specifies a GENIE output filename
         -n specifies the number of events to process (default: all)
         -v z-expansion values to reweight to
         -m reweight normalization of form factor (default: 1)

\author  Aaron Meyer <asmeyer2012 \at uchicago.edu>
         University of Chicago, Fermi National Accelerator Laboratory

         based on gtestRewght by

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created Mar 14, 2016

\cpright Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#include <string>
#include <sstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TArrayF.h>

#include "Algorithm/AlgConfigPool.h"
#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"
#include "ReWeight/GReWeightNuXSecDIS.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GSystUncertainty.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"

// number of parameter tweaks which are coded into reweighting
#define MAX_COEF 4

using namespace genie;
using namespace genie::rew;
using std::string;
using std::ostringstream;

void PrintSyntax();
void GetEventRange      (Long64_t nev_in_file, Long64_t & nfirst, Long64_t & nlast);
void GetCommandLineArgs (int argc, char ** argv);
GSyst_t GetZExpSystematic(int ip);

string gOptInpFilename;
string gOptOutFilename;
Long64_t gOptNEvt1;
Long64_t gOptNEvt2;
bool  gOptDoNorm = false;
double gOptNormValue = 1.;
double gOptParameterValue[MAX_COEF] = {0.};

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
  LOG("rwghtzexpaxff", pNOTICE) << "Input tree header: " << *thdr;
  if(!tree){
    LOG("grwghtzexpaxff", pFATAL)
      << "Can't find a GHEP tree in input file: "<< file.GetName();
    gAbortingInErr = true;
    PrintSyntax();
    exit(1);
  }
  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  Long64_t nev_in_file = tree->GetEntries();
  Long64_t nfirst = 0;
  Long64_t nlast  = 0;
  GetEventRange(nev_in_file, nfirst, nlast);
  int nev = int(nlast - nfirst + 1);

  LOG("rwghtzexpaxff", pNOTICE) << "Will process " << nev << " events";

  //
  // Create a GReWeight object and add to it a set of 
  // weight calculators
  //
  // If seg-faulting here, need to change
  // AxialFormFactorModel in UserPhysicsOptions.xml and LwlynSmithFFCC.xml
  //

  GReWeight rw;
  rw.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );

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

  // Create a concrete weight calculator to fine-tune
  GReWeightNuXSecCCQE * rwccqe = 
    dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccqe"));
  rwccqe->SetMode(GReWeightNuXSecCCQE::kModeZExp);
  // In case uncertainties need to be altered
  GSystUncertainty * unc = GSystUncertainty::Instance();

  // Set up algorithm for loading central values
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("genie::LwlynSmithQELCCPXSec","ZExp");

  Algorithm * alg = algf->AdoptAlgorithm(id);
  XSecAlgorithmI* fXSecModel = dynamic_cast<XSecAlgorithmI*>(alg);
  fXSecModel->AdoptSubstructure();

  Registry * fXSecModelConfig = new Registry(fXSecModel->GetConfig());
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // further optional fine-tuning
  //rwccqe -> RewNue    (false); 
  //rwccqe -> RewNuebar (false); 
  //rwccqe -> RewNumubar(false); 

  // Declare the weights and twkvals arrays 
  const int n_events = (const int) nev;
  const int n_params = (const int) MAX_COEF;
  // if segfaulting here, may need to increase MAX_COEF
  float weights [n_events];

  // Initialize
  for (int iev = 0; iev < nev; iev++) { weights[iev] = 1.; }

  // set first values for weighting
  if (gOptDoNorm)
  {
    //LOG("rwghtzexpaxff", pNOTICE) << "Setting z-expansion tweak for norm : "
    //  << (gOptNormValue-1.)
    //  /unc->OneSigmaErr(kXSecTwkDial_ZNormCCQE,TMath::Sign(1,(int)(gOptNormValue-1.)));
    syst.Set(kXSecTwkDial_ZNormCCQE, (gOptNormValue-1.)
     /unc->OneSigmaErr(kXSecTwkDial_ZNormCCQE,TMath::Sign(1,(int)(gOptNormValue-1.))));
  }
  GSyst_t gsyst;
  double cval = 0.; // central value for parameter
  ostringstream zval_name; // string to use to extract central value
  for (int ipr = 0; ipr < n_params; ipr++)
  {
    gsyst = GetZExpSystematic(ipr+1); // get tweak dial
    zval_name.str("");
    zval_name << "QEL-Z_A" << ipr+1;

    cval = fXSecModelConfig->GetDoubleDef(zval_name.str(),gc->GetDouble(zval_name.str()));
    unc->SetUncertainty(gsyst,1.,1.); // easier to deal with
    
    //LOG("rwghtzexpaxff", pNOTICE) << "Setting z-expansion tweak for param " 
    //  <<ipr<<" : " << (gOptParameterValue[ipr]-cval)/cval;
    syst.Set(gsyst, (gOptParameterValue[ipr]-cval)/cval);
  }

  rw.Reconfigure();
  // Event loop
  for(int iev = nfirst; iev <= nlast; iev++) {
    tree->GetEntry(iev);
  
    EventRecord & event = *(mcrec->event);
    LOG("rwghtzexpaxff", pNOTICE) << "Event number   = " << iev;
    LOG("rwghtzexpaxff", pNOTICE) << event;
  
    double wght = rw.CalcWeight(event);

    LOG("rwghtzexpaxff", pNOTICE) << "Overall weight = " << wght;

    // add to arrays
    weights[iev - nfirst] = wght;
  
    mcrec->Clear();
  } // events

  // Close event file
  file.Close();

  //
  // Save weights 
  //

  // Make an output tree for saving the weights.
  TFile * wght_file = new TFile(gOptOutFilename.c_str(), "RECREATE");
  TTree * wght_tree = new TTree("ZExpCCQE","GENIE weights tree");
  // objects to pass elements into tree
  int branch_eventnum = 0;
  float branch_weight = 1.;
  float branch_norm_val = gOptNormValue;
  float branch_zparam_val[MAX_COEF] = {0.};

  wght_tree->Branch("eventnum", &branch_eventnum);
  wght_tree->Branch("weights",  &branch_weight);
  wght_tree->Branch("norm", &branch_norm_val);

  // create and add branches for each z-expansion coefficient
  ostringstream zparam_brnch_name;
  for (int ipr = 0; ipr < n_params; ipr++) {
    zparam_brnch_name.str("");
    zparam_brnch_name << "param_" << ipr+1;
    LOG("rwghtzexpaxff", pINFO) << "Branch name = " << zparam_brnch_name.str();
    branch_zparam_val[ipr] = gOptParameterValue[ipr];
    LOG("rwghtzexpaxff", pINFO) << "Setting parameter value = " << gOptParameterValue[ipr];
    wght_tree->Branch(zparam_brnch_name.str().c_str(), &branch_zparam_val[ipr]);
  }

  ostringstream str_wght;
  for(int iev = nfirst; iev <= nlast; iev++) {
    branch_eventnum = iev;

    // printout
    LOG("grwghtzexpaxff", pNOTICE)
       << "Filling tree with wght = " << weights[iev - nfirst];
    branch_weight = weights[iev - nfirst];
    wght_tree->Fill();
  }

  wght_file->cd();
  wght_tree->Write();
  wght_tree = 0;
  wght_file->Close();

  // free memory
  delete wght_tree;
  delete wght_file;

  LOG("rwghtzexpaxff", pNOTICE)  << "Done!";
  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("rwghtzexpaxff", pINFO) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {  
    LOG("rwghtzexpaxff", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("rwghtzexpaxff", pFATAL) 
      << "Unspecified input filename - Exiting";
    PrintSyntax();
    exit(1);
  }

  // output weight file
  if(parser.OptionExists('o')) {
    LOG("rwghtzexpaxff", pINFO) << "Reading requested output filename";
    gOptOutFilename = parser.ArgAsString('o');
  } else {
    LOG("rwghtzexpaxff", pINFO) << "Setting default output filename";
    gOptOutFilename = "test_rw_zexp_axff.root";
  }

  if ( parser.OptionExists('n') ) {
    //
    LOG("grwghtzexpaxff", pINFO) << "Reading number of events to analyze";
    string nev =  parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n',",");
      if(vecn.size()!=2) {
         LOG("grwghtzexpaxff", pFATAL) << "Invalid syntax";
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
    LOG("grwghtzexpaxff", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt1 = -1;
    gOptNEvt2 = -1;
  }
  LOG("grwghtzexpaxff", pDEBUG)
    << "Input event range: " << gOptNEvt1 << ", " << gOptNEvt2;

  // values to reweight to:
  if( parser.OptionExists('v') ) {
    LOG("rwghtzexpaxff", pINFO) << "Reading requested parameter values";
    string coef = parser.ArgAsString('v');
    
    vector<string> zvals = utils::str::Split(coef, ",");
    // MAX_COEF should be set to the number of z-expansion tweaks which exist
    assert(zvals.size() == (unsigned int) MAX_COEF);
    for (int ik = 0;ik<MAX_COEF;ik++)
    {
      gOptParameterValue[ik] = atof(zvals[ik].c_str());
    }
    for (int ik = 0;ik<MAX_COEF;ik++)
    {
      LOG("rwghtzexpaxff",pINFO)<<"Parameter value "<<ik+1<<": "<<
       gOptParameterValue[ik];
    }
  }

  // norm to reweight to:
  if( parser.OptionExists('m') ) {
    LOG("rwghtzexpaxff", pINFO) << "Reading requested normalization";
    string coef = parser.ArgAsString('m');
    gOptDoNorm = true;
    gOptNormValue = atof(coef.c_str());
    LOG("rwghtzexpaxff",pINFO)<<"Requested normalization: "<< gOptNormValue;
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

  assert(nfirst <= nlast && nfirst >= 0 && nlast <= nev_in_file-1);
}
//_________________________________________________________________________________
GSyst_t GetZExpSystematic(int ip)
{
    switch(ip){
      case 1: return kXSecTwkDial_ZExpA1CCQE; break;
      case 2: return kXSecTwkDial_ZExpA2CCQE; break;
      case 3: return kXSecTwkDial_ZExpA3CCQE; break;
      case 4: return kXSecTwkDial_ZExpA4CCQE; break;
      default:
        LOG("rwghtzexpaxff", pFATAL) 
          << "Cannot find systematic corresponding to parameter " << ip;
        exit(0);
        break;
    }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("grwghtzexpaxff", pFATAL)
     << "\n\n"
     << "grwghtzexpaxff               \n"
     << "     -f input_event_file     \n"
     << "     -v val1,val2,val3,val4  \n"
     << "    [-n nev]                 \n"
     << "    [-o output_weights_file]";
}
//_________________________________________________________________________________
