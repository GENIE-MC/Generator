//____________________________________________________________________________
/*!

\program gRwghtZExpAxFF

\brief   A simple program to illustrate how to use the GENIE event reweighting
         for use with the z-expansion axial form factor

\syntax  grwghtzexpaxff -f filename -t NTwk1,NTwk2,... [-n nev] [-o fileOutName]
         [-s SigmaLo1,SigmaHi1,SigmaLo2,SigmaHi2,...] [-m NTwkN]

         where 
         [] is an optional argument
         -f specifies a GENIE event file (GHEP format)
         -o specifies a GENIE output filename
         -n specifies the number of events to process (default: all)
         -t specify number of tweaks on each z-expansion coefficient
            values are comma separated (# < 2 are ignored)
         -s specify +- one-sigma bounds on all coefficients up to max
            values are comma separated, given as percentages
            requires 2x number of fields from -t option
            default value is 10% on all coefficients
         -m number of tweaks on normalization
            puts reweighting into norm+shape mode

\author  Aaron Meyer <asmeyer2012 \at uchicago.edu>
         University of Chicago, Fermi National Accelerator Laboratory

         based on gtestRewght by

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created Dec 26, 2014

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

// number of coefficient values to vary
#define MAX_COEF 4

using namespace genie;
using namespace genie::rew;
using std::string;
using std::ostringstream;

void PrintSyntax();
void GetEventRange      (Long64_t nev_in_file, Long64_t & nfirst, Long64_t & nlast);
void GetCommandLineArgs (int argc, char ** argv);
int  GetNumberOfWeights (int* ntwk, int kmaxinc, int normtwk, bool donorm);
bool IncrementCoefficients(int* ntwk, int kmaxinc, int normtwk, bool donorm,
                         float* twkvals, GSystSet& syst);
GSyst_t GetZExpSystematic(int ip);

string gOptInpFilename;
string gOptOutFilename;
//int    gOptNEvt;
Long64_t gOptNEvt1;
Long64_t gOptNEvt2;
int    gOptKmaxInc = 0;
int    gOptNormTweaks = 0;
bool   gOptDoNorm = false; // whether to be in norm+shape mode or not
bool   gOptSigmaDefined = false; // handles setting of SigMin, SigMax
int    gOptNTweaks[MAX_COEF] = {0 };
float  gOptSigMin [MAX_COEF] = {0.};
float  gOptSigMax [MAX_COEF] = {0.};

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

  // further optional fine-tuning
  //rwccqe -> RewNue    (false); 
  //rwccqe -> RewNuebar (false); 
  //rwccqe -> RewNumubar(false); 

  // Declare the weights and twkvals arrays 
  const int n_events = (const int) nev;
  const int n_params = (const int) (gOptKmaxInc + 1); // +1 for norm
  const int n_points =
   (const int) GetNumberOfWeights(gOptNTweaks,gOptKmaxInc,gOptNormTweaks,gOptDoNorm);
  // if segfaulting here, may need to increase MAX_COEF
  // copied from gtestRwght... seems inefficient
  // -- couldn't we just load straight to the tree? doing so would prevent segfaults
  float weights  [n_events][n_points];
  float twkvals  [n_points][n_params];

  // Initialize
  for (int ipt = 0; ipt < n_points; ipt++)
  {
    for (int iev = 0; iev < nev; iev++) { weights[iev][ipt] = 1.; }
    twkvals[ipt][0] = (gOptDoNorm && (gOptNormTweaks > 1)) ? -1 : 0;
    for (int ipr = 1; ipr < n_params; ipr++)
    {
      twkvals[ipt][ipr] = (gOptNTweaks[ipr-1] > 1 ? -1 : 0);
    }
  }

  // set first values for weighting
  if (gOptDoNorm)
  {
    syst.Set(kXSecTwkDial_ZNormCCQE, twkvals[0][0]);
    LOG("rwghtzexpaxff", pNOTICE) << "Setting z-expansion tweak for norm : "
      << twkvals[0][0];
  }
  GSyst_t gsyst;
  for (int ipr = 1; ipr < n_params; ipr++)
  {
    gsyst = GetZExpSystematic(ipr);
    syst.Set(gsyst, twkvals[0][ipr]);
    LOG("rwghtzexpaxff", pNOTICE) << "Setting z-expansion tweak for param " 
      <<ipr<<" : " << twkvals[0][ipr];
    if (gOptSigmaDefined)
    {
      unc->SetUncertainty(gsyst,gOptSigMin[ipr-1],gOptSigMax[ipr-1]);
      LOG("rwghtzexpaxff", pNOTICE) << "Setting z-expansion sigma for param " 
        <<ipr<<" : " << gOptSigMin[ipr-1] <<","<< gOptSigMax[ipr-1];
    }
  }

  // point loop (number of parameter combinations)
  for (int ipt = 0; ipt < n_points; ipt++) {
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
      weights[iev - nfirst][ipt] = wght;
  
      mcrec->Clear();
    } // events

    // set the next set of coefficients to match the previous, then increment
    if (ipt < n_points-1) {
      for (int ipr=0;ipr<n_params;ipr++)
      {
        twkvals[ipt+1][ipr] = twkvals[ipt][ipr];
      }
      IncrementCoefficients(gOptNTweaks,n_params,gOptNormTweaks,gOptDoNorm,
        twkvals[ipt+1],syst);
    }
  }   // points

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
  TArrayF *  branch_weight_array   = new TArrayF(n_points);
  TArrayF ** branch_twkdials_array = new TArrayF* [n_params];

  wght_tree->Branch("eventnum", &branch_eventnum);
  wght_tree->Branch("weights",  &branch_weight_array);

  // create and add branches for each z-expansion coefficient
  ostringstream twk_dial_brnch_name;
  for (int ipr = 0; ipr < n_params; ipr++) {
    if (!gOptDoNorm && ipr == 0) { continue; } // skip norm if not requested
    twk_dial_brnch_name.str("");
    if (ipr == 0) { twk_dial_brnch_name << "twk_dial_param_norm";      }
    else          { twk_dial_brnch_name << "twk_dial_param_" << ipr; }
    LOG("rwghtzexpaxff", pWARN) << "Branch name = " << twk_dial_brnch_name.str();
    branch_twkdials_array[ipr] = new TArrayF(n_points);
    wght_tree->Branch(twk_dial_brnch_name.str().c_str(), branch_twkdials_array[ipr]);
  }

  ostringstream str_wght;
  for(int iev = nfirst; iev <= nlast; iev++) {
    branch_eventnum = iev;

    for(int ipt = 0; ipt < n_points; ipt++){

       // printout
       str_wght.str("");
       str_wght << ", tweaked parameter values : ";
       for (int ipr = 0; ipr < n_params; ipr++) {
          if (ipr > 0) str_wght << ", ";
          str_wght << ipr << " -> " << twkvals[ipt][ipr];
       }
       LOG("grwghtzexpaxff", pNOTICE)
          << "Filling tree with wght = " << weights[iev - nfirst][ipt] << str_wght.str();

       // fill tree
       branch_weight_array   -> AddAt (weights [iev - nfirst][ipt], ipt);
       for (int ipr = (gOptDoNorm ? 0 : 1); ipr < n_params; ipr++) // skip norm if not requested
         { branch_twkdials_array[ipr] -> AddAt (twkvals[ipt][ipr], ipt); }

    } // twk_dial loop
    wght_tree->Fill();
  }

  wght_file->cd();
  wght_tree->Write();
  wght_tree = 0;
  wght_file->Close();

  // free memory
  delete wght_tree;
  delete wght_file;
  for (int ipr = 0; ipr < n_params; ipr++) {
    if (!gOptDoNorm && ipr == 0) { continue; }
    delete branch_twkdials_array[ipr];
  }
  delete branch_twkdials_array;
  delete branch_weight_array;

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

  // number of tweaks:
  if( parser.OptionExists('t') ) {
    LOG("rwghtzexpaxff", pINFO) << "Reading number of tweaks";
    string coef = parser.ArgAsString('t');
    
    // split into sections of min,max,inc(rement)
    vector<string> coefrange = utils::str::Split(coef, ",");
    gOptKmaxInc = coefrange.size();
    if (gOptKmaxInc > MAX_COEF)
    {
      LOG("rwghtzexpaxff", pFATAL) 
        << "Too many coefficients: increase MAX_COEF in source code and recompile";
      exit(1);
    }
    //AlgConfigPool * confp = AlgConfigPool::Instance();
    //const Registry * gc = confp->GlobalParameterList();
    //int kDef = fConfig->GetIntDef("QEL-Kmax",gc->GetInt("QEL-Kmax"));
    //else if (gOptKmaxInc > kDef)
    //{
    //  LOG("rwghtzexpaxff", pFATAL) 
    //    << "Too many coefficients:
    //        requested number of coefficients is more than defined in UserPhysicsOptions.xml";
    //  exit(1);
    //}

    LOG("rwghtzexpaxff", pINFO) << "Largest coefficient to tweak : " << gOptKmaxInc;
    for (int ik = 0;ik<gOptKmaxInc;ik++)
    {
      gOptNTweaks[ik] = atof(coefrange[ik].c_str());
    }
    for (int ik = 0;ik<gOptKmaxInc;ik++)
    {
      LOG("rwghtzexpaxff",pINFO)<<"Number of tweaks on coefficient "<<ik+1<<" : "<< gOptNTweaks[ik];
    }
  } else {
    LOG("rwghtzexpaxff", pFATAL) 
      << "Unspecified tweaks for parameters - Exiting";
    PrintSyntax();
    exit(1);
  }

  // lower/upper sigma:
  if( parser.OptionExists('s') ) {
    LOG("rwghtzexpaxff", pINFO) << "Reading specified parameter uncertainties";
    string coef = parser.ArgAsString('s');
    
    // split into sections of min,max
    vector<string> sigrange = utils::str::Split(coef, ",");
    // gOptKmaxInc defined by number of tweaks (-t)
    assert(sigrange.size() == (unsigned int) 2*gOptKmaxInc);
    gOptSigmaDefined = true;
    for (int ik = 0;ik<gOptKmaxInc;ik++)
    {
      gOptSigMin[ik] = atof(sigrange[ik*2  ].c_str());
      gOptSigMax[ik] = atof(sigrange[ik*2+1].c_str());
    }
    for (int ik = 0;ik<gOptKmaxInc;ik++)
    {
      LOG("rwghtzexpaxff",pINFO)<<ik+1<<": "<< gOptSigMin[ik] <<","<< gOptSigMax[ik];
    }
  }

  // number of norm tweaks:
  if( parser.OptionExists('m') ) {
    LOG("rwghtzexpaxff", pINFO) << "Reading number of tweaks on normalization";
    string coef = parser.ArgAsString('m');
    gOptDoNorm = true;
    gOptNormTweaks   = atof(coef.c_str());
    LOG("rwghtzexpaxff",pINFO)<<"Number of tweaks on norm : "<< gOptNormTweaks;
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
bool IncrementCoefficients(int* ntwk, int kmaxinc, int normtwk, bool donorm,
                           float* twkvals, GSystSet& syst) 
{
  if (kmaxinc < 2 && ! donorm)
  {
    LOG("grwghtzexpaxff",pERROR) << "No coefficients to increment";
    return false;
  } else {

  int ip = -1;
  bool stopflag = false;
  GSyst_t gsyst = kXSecTwkDial_ZNormCCQE;
  do
  {
    if (ip > 0 || (ip == 0 && donorm))
    { // fully incremented a coefficient; reset it and continue to the next
      if (ip == 0) { twkvals[0 ] = (normtwk    > 1 ? -1. : 0.); }
      else         { twkvals[ip] = (ntwk[ip-1] > 1 ? -1. : 0.); }

      // set the value manually
      syst.Set(gsyst, twkvals[ip]);
      LOG("rwghtzexpaxff", pNOTICE) << "Setting z-expansion tweak for param " 
        <<ip<<" : " << twkvals[ip];
    }
    stopflag = true;

    ip++;                                 // increment index
    if (ip == kmaxinc) { return false; }  // done with incrementing
    if (ip == 0 && ! donorm)
    {
      stopflag = false;
      continue;  // skip when not doing norm
    }
    // set to next systematic
    if (ip == 0) { gsyst = kXSecTwkDial_ZNormCCQE; }
    else         { gsyst = GetZExpSystematic(ip); }

    // increment systematic
    if (ip == 0)
    {
      if (normtwk > 1) { twkvals[0] += 2./float(normtwk-1); }
      else             { stopflag = false; continue; }
    }
    else
    {
      if (ntwk[ip-1] > 1) { twkvals[ip] += 2./float(ntwk[ip-1]-1); }
      else                { stopflag = false; continue; }
    }

    // set the systematic to the new tweak value
    // actual updating of this will be handled by GReWeight::Reconfigure
    syst.Set(gsyst, twkvals[ip]);
    LOG("rwghtzexpaxff", pNOTICE) << "Setting z-expansion tweak for param " 
      <<ip<<" : " << twkvals[ip];
    if (twkvals[ip] > 1. + controls::kASmallNum) { stopflag=false; } // went over

  } while (! stopflag); // loop

  return true;
  } // if kmaxinc >= 1

  return false;
}
//_________________________________________________________________________________
int GetNumberOfWeights(int* ntwk, int kmaxinc, int normtwk, bool donorm)
{
  int  num_pts = 1;
  for (int i=0;i<kmaxinc;i++)
  {
    if (ntwk[i] > 1) num_pts *= ntwk[i];
  }
  if (donorm)
  {
    if (normtwk > 1) num_pts *= normtwk;
  }
  return num_pts;
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
     << "     -t ntwk1[,ntwk2[,...]]  \n"
     << "    [-n nev]                 \n"
     << "    [-s sigLo1,sigHi1[,...]] \n"
     << "    [-o output_weights_file] \n"
     << "    [-m ntwkNorm]" ;
}
//_________________________________________________________________________________
