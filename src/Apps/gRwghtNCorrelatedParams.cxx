//____________________________________________________________________________
/*!

\program grwghtnp

\brief   Generates weights given an input GHEP event file, a set of (possibly
         correlated) systematic parameters (supported by the ReWeight package),
         and a covariance matrix for the input set of parameters.
         The covariance matrix should be in the form a ROOT file containing
         only a TMatrixD object which is square and symmetric.
         It outputs a ROOT file containing a tree with an entry for every 
         input event. Each such tree entry contains a TArrayD of all computed 
         weights and TArrayD for each requested systematic of all of the
         corresponding randomly generated tweak dial values.

\syntax  grwghtcov \
           -f input_event_file 
           -c input_covariance_file
           -s systematic1[,systematic2[,...]] 
           -v central_value1[,central_value2[,...]]
           -t n_twk_dial_values
          [-n n1[,n2]] 
          [-r run_key]
          [-o output_weights_file]

         where 
         [] is an optional argument.

         -f 
            Specifies a GHEP input file.
         -c
            Specifies a binary ROOT file which contains the covariance matrix
            as a TMatrixD object.
         -s 
            Specifies the systematic parameters to tweak.
            See $GENIE/src/ReWeight/GSyst.h for a list of parameters and
            their corresponding label, which is what should be input here.
         -v
            Central values specified in $GENIE/config/UserPhysicsOptions.xml
            for the reweighted parameters. Assigns uncertainties based on
            covariance matrix diagonal. Should be automated, but cannot for now.
         -t 
            Number of random drawings of tweak values between -1 and 1.
            Values for tweaks respect the covariance of systematics
         -n 
            Specifies an event range.
            Examples:
            - Type `-n 50,2350' to process all 2301 events from 50 up to 2350.
              Note: Both 50 and 2350 are included.
            - Type `-n 1000' to process the first 1000 events;
              from event number 0 up to event number 999.
            This is an optional argument. 
            By default GENIE will process all events.
         -r
            Specifies an integer run key.
            Changes temporary file names so that multiple instances can run
            without overwriting each other's temporary tree files

\author  Aaron Meyer <asmeyer2012 \at uchicago.edu>
         University of Chicago, Fermi National Accelerator Laboratory

\created July 20, 2015

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________


#include <TArrayD.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TRandom.h>

#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Numerical/RandomGen.h"
#include "Messenger/Messenger.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"
#include "ReWeight/GReWeightNuXSecDIS.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GSystUncertainty.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/MathUtils.h"
#include "Utils/StringUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::rew;
using namespace genie::utils::math;
using std::stringstream;

void PrintSyntax();
void GetEventRange       (Long64_t nev_in_file, Long64_t & nfirst, Long64_t & nlast);
void GetCommandLineArgs  (int argc, char ** argv);
void GetCorrelationMatrix(string fname, TMatrixD *& cmat);
void AdoptWeightCalcs    (vector<GSyst_t> lsyst, GReWeight & rw);
bool FindIncompatibleSystematics(vector<GSyst_t> lsyst);

vector<GSyst_t> gOptVSyst;
vector<double>  gOptVCentVal;
string   gOptInpFilename;
string   gOptInpCovariance;
string   gOptOutFilename;
Long64_t gOptNEvt1;
Long64_t gOptNEvt2;
int      gOptRunKey= 0;
int      gOptNSyst = 0;
int      gOptNTwk  = 0;
TRandom *tRnd = new TRandom(); // to access normal distribution

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
  if(!tree){
    LOG("grwghtcov", pFATAL)
      << "Can't find a GHEP tree in input file: "<< file.GetName();
    gAbortingInErr = true;
    PrintSyntax();
    exit(1);
  }
  LOG("grwghtcov", pNOTICE) << "Input tree header: " << *thdr;
  if(!FindIncompatibleSystematics(gOptVSyst))
  {
    LOG("grwghtcov", pFATAL) << "Error: conflicting systematics";
    gAbortingInErr = true;
    exit(1);
  }

  //
  // Preparation for finding correlated vectors
  //
  // Solutions are subject to constraint of an error ellipse with equation:
  //  1 = x^T.Cor^-1.x 
  // x is a vector of tweaks to apply to the systematics
  // Cholesky Decomposition solves for U in equation:
  //  Cor^-1 = U^T.U
  // LU decomposition used to solve for x:
  //  L.U.x = b
  // where L=I and b is a vector of random numbers with b^T.b = 1
  //

  TMatrixD *cmat = NULL;
  // Gets Cor, which is needed in decompositions
  // Assumed errors from covariance are stored in one sigma errors for parameters
  GetCorrelationMatrix(gOptInpCovariance,cmat);
  TMatrixD lTri = CholeskyDecomposition(*cmat);

  //LOG("grwghtcov", pNOTICE) << "Correlation matrix:";
  //cmat->Print();
  //LOG("grwghtcov", pNOTICE) << "Lower triangular matrix:";
  //lTri.Print();

  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  Long64_t nev_in_file = tree->GetEntries();
  Long64_t nfirst = 0;
  Long64_t nlast  = 0;
  GetEventRange(nev_in_file, nfirst, nlast);
  int nev = int(nlast - nfirst + 1);

  LOG("grwghtcov", pNOTICE) << "Will process " << nev << " events";

  //
  // Create a GReWeight object and add to it a set of 
  // weight calculators
  //
  // If seg-faulting here, probably need to change
  // models in UserPhysicsOptions.xml and other config files
  //

  GReWeight rw;
  AdoptWeightCalcs(gOptVSyst, rw);

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

  // Declare the weights, twkvals
  const int n_params = (const int) gOptNSyst;
  const int n_tweaks = (const int) gOptNTwk;
  TVectorD twkvals(n_params);

  // Initialize
  for (int ipr = 0; ipr < n_params; ipr++) { twkvals(ipr) = 0.; }

  // objects to pass elements into tree
  int     branch_eventnum = 0;
  double  branch_weight   = 0.;

  // objects used in processing
  vector<GSyst_t>::iterator it;
  stringstream twk_dial_brnch_name;
  stringstream tmpName;
  int ip;

  //
  // REWEIGHTING LOOP
  // -- do all of reweighting, save to temporary files
  //
  TFile * wght_file = NULL;
  TTree * wght_tree = NULL;
  for (int itk = 0; itk < gOptNTwk; itk++) {
    // Make temporary output trees for saving the weights.
    // This step is necessary because ROOT trees cannot be edited once filled
    // Later consolidate the trees into a single tree with the requested filename
    tmpName.str("");
    tmpName << "_temporary_rwght." <<itk <<"." <<gOptRunKey <<".root";
    LOG("grwghtcov", pINFO) <<"temporary file: " <<tmpName.str();
    wght_file = new TFile(tmpName.str().c_str(),"RECREATE");
    wght_tree = new TTree("covrwt","GENIE weights tree");

    // Create tree branches
    wght_tree->Branch("eventnum", &branch_eventnum);
    wght_tree->Branch("weights",  &branch_weight);

    // Construct multiple branches to streamline loading later
    // Load tweaks into reweighting
    twkvals = CholeskyGenerateCorrelatedParamVariations(lTri);
    // scale the size of the vector to avg length 1
    //twkvals *= 1./(TMath::Sqrt((double)gOptNSyst));
    ip = 0;
    for (it = gOptVSyst.begin();it != gOptVSyst.end(); it++, ip++) {
      twk_dial_brnch_name.str("");
      twk_dial_brnch_name << "twk_" << GSyst::AsString(*it);
      // each array element individually
      wght_tree->Branch(twk_dial_brnch_name.str().c_str(), &twkvals(ip));
      //LOG("grwghtcov", pINFO) << "Setting systematic : "
      //  <<GSyst::AsString(*it) <<", " <<twkvals(ip);
      syst.Set(*it,twkvals(ip));
    }
    rw.Reconfigure();

    stringstream str_wght;
    str_wght.str("");
    str_wght << ", parameter tweaks : ";
    for (int ipr=0; ipr < n_params; ipr++) {
       if (ipr > 0) str_wght << ", ";
       str_wght << ipr << " -> " << twkvals(ipr);
    }

    for(int iev = nfirst; iev <= nlast; iev++) {
      branch_eventnum = iev;
      tree->GetEntry(iev);

      EventRecord & event = *(mcrec->event);
      LOG("rwghtzexpaxff", pNOTICE) << "Event_num  => " << iev;
      //LOG("rwghtzexpaxff", pNOTICE) << event;

      branch_weight = rw.CalcWeight(event);
      mcrec->Clear();
      wght_tree->Fill();

    } // event loop

    // close out temporary file
    wght_file->cd();
    wght_tree->Write();
    wght_file->Close();
    //delete wght_tree; // segfault when deleted
    wght_tree = 0;
    delete wght_file;
  } // tweak loop

  // Close event file
  file.Close();

  // open temporary trees for consolidation
  LOG("rwghtzexpaxff", pNOTICE)
    << "Consolidating temporary files into ROOT file " << gOptOutFilename;
  wght_file = new TFile(gOptOutFilename.c_str(),"RECREATE"); // new file
  wght_tree = new TTree("covrwt","GENIE covariant reweighting tree");
  wght_tree->Branch("n_tweaks", &gOptNTwk);
  wght_tree->Branch("eventnum", &branch_eventnum);
  TFile * file_list[n_tweaks];
  TTree * wght_list[n_tweaks];
  for (int itk=0; itk < n_tweaks; itk++) {
    tmpName.str("");
    tmpName << "_temporary_rwght." <<itk <<"." <<gOptRunKey <<".root";
    file_list[itk] = new TFile(tmpName.str().c_str(),"READ");
    wght_list[itk] = (TTree*)file_list[itk]->Get("covrwt");
  }

  // objects to load data into and fill new tree with
  TArrayD   branch_weights_array(gOptNTwk);
  double  * branch_weights_ptr = branch_weights_array.GetArray();
  TArrayD * branch_twkdials_array[n_params];
  double  * branch_twkdials_ptr  [n_params];

  // set up streamlined weight loading
  wght_tree->Branch("weights",  &branch_weights_array);
  for (int itk = 0; itk < gOptNTwk; itk++) {
    wght_list[itk]->SetBranchAddress("weights",&branch_weights_ptr[itk]);
  }

  ip = 0;
  for (it = gOptVSyst.begin();it != gOptVSyst.end(); it++, ip++) {
    twk_dial_brnch_name.str("");
    twk_dial_brnch_name << "twk_" << GSyst::AsString(*it);

    // access TArrayD memory directly
    branch_twkdials_array[ip] = new TArrayD(gOptNTwk); 
    branch_twkdials_ptr[ip] = branch_twkdials_array[ip]->GetArray();

    // create branch
    wght_tree->Branch(twk_dial_brnch_name.str().c_str(), branch_twkdials_array[ip]);
    LOG("grwghtcov", pINFO) << "Creating tweak branch : " << twk_dial_brnch_name.str();
 
    // set up loading directly into TArrayD
    for (int i=0; i < n_tweaks; i++) { 
      wght_list[i]->SetBranchAddress(twk_dial_brnch_name.str().c_str(),&branch_twkdials_ptr[ip][i]);
      //LOG("grwghtcov", pINFO) << "Loading tweak value : "<<branch_twkdials_array[ip]->fArray[i];
    }
  }
 
  //
  // CONSOLIDATION LOOP
  // -- combine all data from reweighting into single file
  //
  wght_file->cd();
  for(int iev = nfirst; iev <= nlast; iev++) {
    branch_eventnum = iev;
    for (int itk = 0; itk < n_tweaks; itk++) {
      wght_list[itk]->GetEntry(iev);
    } // tweak loop
    wght_tree->Fill();
  } // event loop
  wght_file->cd();
  wght_tree->Write();

  //
  // CLEANUP LOOP
  //

  // delete temporary files
  for (int itk = 0; itk < gOptNTwk; itk++) {
    tmpName.str("");
    tmpName << "_temporary_rwght." <<itk <<"." <<gOptRunKey <<".root";
    if( remove(tmpName.str().c_str()) != 0 )
    { LOG("grwghtcov", pWARN) << "Could not delete temporary file : " << tmpName.str(); }
    //else 
    //{ LOG("grwghtcov", pINFO) << "Deleted temporary file : " << tmpName.str(); }
  }

  // free memory
  wght_file->Close();
  //delete wght_tree;
  wght_tree = 0;
  delete wght_file;
  for (int ipr = 0; ipr < n_params; ipr++) {
    delete branch_twkdials_array[ipr];
  }

  LOG("grwghtcov", pNOTICE)  << "Done!";
  return 0;
}
//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("grwghtcov", pINFO) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {  
    LOG("grwghtcov", pINFO) << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("grwghtcov", pFATAL) 
      << "Unspecified input filename - Exiting";
    PrintSyntax();
    exit(1);
  }

  // get ROOT covariance matrix binary file
  if( parser.OptionExists('c') ) {  
    LOG("grwghtcov", pINFO) << "Reading covariance matrix filename";
    gOptInpCovariance = parser.ArgAsString('c');
  } else {
    LOG("grwghtcov", pFATAL) 
      << "Unspecified covariance filename - Exiting";
    PrintSyntax();
    exit(1);
  }

  // output weight file
  if(parser.OptionExists('o')) {
    LOG("grwghtcov", pINFO) << "Reading requested output filename";
    gOptOutFilename = parser.ArgAsString('o');
  } else {
    LOG("grwghtcov", pINFO) << "Setting default output filename";
    gOptOutFilename = "rwt_cov.root";
  }

  if ( parser.OptionExists('n') ) {
    //
    LOG("grwghtcov", pINFO) << "Reading number of events to analyze";
    string nev =  parser.ArgAsString('n');
    if (nev.find(",") != string::npos) {
      vector<long> vecn = parser.ArgAsLongTokens('n',",");
      if(vecn.size()!=2) {
         LOG("grwghtcov", pFATAL) << "Invalid syntax";
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
    LOG("grwghtcov", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt1 = -1;
    gOptNEvt2 = -1;
  }
  LOG("grwghtcov", pDEBUG)
    << "Input event range: " << gOptNEvt1 << ", " << gOptNEvt2;

  // systematics:
  if( parser.OptionExists('s') ) {
    LOG("grwghtcov", pINFO) << "Reading systematics";
    string insyst = parser.ArgAsString('s');
    vector<string> lsyst = utils::str::Split(insyst, ",");
    vector<string>::iterator it;
    int ik = 0;
    for(it=lsyst.begin();it != lsyst.end();it++,ik++)
    {
      gOptVSyst.push_back(GSyst::FromString(*it));
      LOG("grwghtcov",pINFO)<<"Read systematic "<<ik+1<<" : "<< lsyst[ik];
    }
    // split into strings of systematics
    gOptNSyst = gOptVSyst.size();
    LOG("grwghtcov", pINFO) << "Number of systematics : " << gOptNSyst;
  } else {
    LOG("rwghtcov", pFATAL) 
      << "Unspecified systematics - Exiting";
    PrintSyntax();
    exit(1);
  }

  // systematic central values:
  if( parser.OptionExists('v') ) {
    LOG("grwghtcov", pINFO) << "Reading parameter central values";
    gOptVCentVal  = parser.ArgAsDoubleTokens('v',",");
    // check size
    if (gOptNSyst != (int)gOptVCentVal.size()) {
      LOG("rwghtcov", pFATAL) 
        << "Number of systematic central values does not match number of systematics- Exiting";
      PrintSyntax();
      exit(1);
    }
  } else {
    LOG("rwghtcov", pFATAL) 
      << "Unspecified systematic central values - Exiting";
    PrintSyntax();
    exit(1);
  }

  // number of tweaks:
  if( parser.OptionExists('t') ) {
    LOG("grwghtcov", pINFO) << "Reading number of tweaks";
    gOptNTwk = parser.ArgAsInt('t');
    
    if( gOptNTwk < 1 )
    {
      LOG("grwghtcov", pFATAL) << "Must have at least 1 tweak - Exiting";
      PrintSyntax();
      exit(1);
    }
    LOG("grwghtcov",pINFO)<<"Number of tweaks : "<< gOptNTwk;
  } else {
    LOG("grwghtcov", pFATAL) 
      << "Unspecified tweaks for parameters - Exiting";
    PrintSyntax();
    exit(1);
  }

  // run key:
  if( parser.OptionExists('r') ) {
    LOG("grwghtcov", pINFO) << "Reading run key";
    int gOptRunKey = parser.ArgAsInt('r');
    
    LOG("grwghtcov", pINFO) 
      << "Run key set to " <<gOptRunKey;
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
void GetCorrelationMatrix(string fname, TMatrixD *& cmat)
{
  //
  // Loads a ROOT file containing only a TMatrixD
  // Reads a covariance/correlation matrix and returns correlation matrix
  // 

  // Load covariance matrix diagonal into uncertainties
  GSystUncertainty * unc = GSystUncertainty::Instance();
  vector<GSyst_t>::iterator it;
  vector<double>::iterator itd;

  // open file and prepare for loading
  TFile * fin = new TFile(fname.c_str(),"READ");
  TKey * fkey;
  TMatrixD * inmat = new TMatrixD();
  if(cmat) { delete cmat; }

  // read covariance into cmat and return
  TIter next(fin->GetListOfKeys());
  fkey = (TKey *)next();
  fkey->Read(inmat);

  // checks 
  fin->Close();
  if (!inmat) {
    LOG("grwghtcov", pFATAL) << "Error reading covariance matrix - Exiting";
    gAbortingInErr = true;
    exit(1);
  }
  if (inmat->GetNrows() != inmat->GetNcols()) {
    LOG("grwghtcov", pFATAL) << "Covariance matrix not square - Exiting";
    gAbortingInErr = true;
    exit(1);
  }
  if (inmat->GetNrows() != gOptNSyst ){
    LOG("grwghtcov", pFATAL) << "Number of systematics does not match covariance matrix size- Exiting";
    gAbortingInErr = true;
    exit(1);
  }

  // Convert covariance to correlation
  TMatrixD tmpmat = TMatrixD(*inmat);

  // Diag = Diag(Cov)
  // set off-diagonals to zero and take 1/square root of diagonals
  // is there an easier way?
  int i=0;
  itd = gOptVCentVal.begin();
  for (it = gOptVSyst.begin();it != gOptVSyst.end(); it++, itd++, i++) {
    for(int j=0;j<tmpmat.GetNcols();j++){
      if(i!=j) { tmpmat(i,j) = 0.; }
      else     {
        // convert diagonal entries to uncertainty
        //LOG("grwghtcov", pINFO) <<"Setting uncertainty "<<i<<" to : "<<
        //  TMath::Sqrt(tmpmat(i,i))/(*itd);
        unc->SetUncertainty(*it,TMath::Sqrt(tmpmat(i,i))/(*itd),
          TMath::Sqrt(tmpmat(i,i))/(*itd));
        tmpmat(i,i) = 1./TMath::Sqrt(tmpmat(i,i));
      }
    }
  }
  //LOG("grwghtcov", pWARN) <<"diagonal scaling matrix:";
  //tmpmat.Print();

  // Cor = Diag^(-1/2).Cov.Diag^(-1/2)
  TMatrixD sigmat = TMatrixD(tmpmat);
  tmpmat = TMatrixD(*inmat,TMatrixD::kMult,sigmat);
  cmat = new TMatrixD(sigmat,TMatrixD::kMult,tmpmat);

  delete inmat;
  return;
}
//_________________________________________________________________________________
bool FindIncompatibleSystematics(vector<GSyst_t> lsyst)
{
  //
  // returns false if there are incompatible systematics
  //
  vector<GSyst_t>::iterator it0;
  vector<GSyst_t>::iterator it1;
  // CC QE
  GSyst_t ccqe_ma_shp[1]
    = { kXSecTwkDial_MaCCQE };
  GSyst_t ccqe_ma_shp_norm[2]
    = { kXSecTwkDial_NormCCQE,    kXSecTwkDial_MaCCQEshape };
  GSyst_t ccqe_z_shp[5]
    = { kXSecTwkDial_ZNormCCQE ,  kXSecTwkDial_ZExpA1CCQE,
        kXSecTwkDial_ZExpA2CCQE,  kXSecTwkDial_ZExpA3CCQE,
        kXSecTwkDial_ZExpA4CCQE };
  // CC Res
  GSyst_t ccres_shp[2]
    = { kXSecTwkDial_MaCCRES,     kXSecTwkDial_MvCCRES };
  GSyst_t ccres_shp_norm[3]
    = { kXSecTwkDial_NormCCRES,   kXSecTwkDial_MaCCRESshape,
        kXSecTwkDial_MvCCRESshape };
  // DIS
  GSyst_t dis_shp_norm[4]
    = { kXSecTwkDial_AhtBY,       kXSecTwkDial_BhtBY,
        kXSecTwkDial_CV1uBY,      kXSecTwkDial_CV2uBY };
  GSyst_t dis_shp[4]
    = { kXSecTwkDial_AhtBYshape,  kXSecTwkDial_BhtBYshape,
        kXSecTwkDial_CV1uBYshape, kXSecTwkDial_CV2uBYshape };
  // NC Res
  GSyst_t ncres_shp[2]
    = { kXSecTwkDial_MaNCRES,     kXSecTwkDial_MvNCRES };
  GSyst_t ncres_shp_norm[3]
    = { kXSecTwkDial_NormNCRES,   kXSecTwkDial_MaNCRESshape,
        kXSecTwkDial_MvNCRESshape };

  //
  // Sort systematics by mode (it0)
  // Get the list of systematics for that systematic
  // Look through the rest of systematics the list for others with similar modes (it1)
  // If it1 conflicts with the group for it0, reject whole list
  //
  for(it0=lsyst.begin();it0 != lsyst.end();it0++)
  {
    switch(*it0){
    // CC QE
    case kXSecTwkDial_MaCCQE:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<2;i++){ if(ccqe_ma_shp_norm[i] == *it1) return false; }
        for(int i=0;i<5;i++){ if(ccqe_z_shp[i]       == *it1) return false; }
      }
    case kXSecTwkDial_NormCCQE:
    case kXSecTwkDial_MaCCQEshape:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<1;i++){ if(ccqe_ma_shp[i]      == *it1) return false; }
        for(int i=0;i<5;i++){ if(ccqe_z_shp[i]       == *it1) return false; }
      }
    break;
    case kXSecTwkDial_ZNormCCQE:
    case kXSecTwkDial_ZExpA1CCQE:
    case kXSecTwkDial_ZExpA2CCQE:
    case kXSecTwkDial_ZExpA3CCQE:
    case kXSecTwkDial_ZExpA4CCQE:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<1;i++){ if(ccqe_ma_shp[i]      == *it1) return false; }
        for(int i=0;i<2;i++){ if(ccqe_ma_shp_norm[i] == *it1) return false; }
      }
    break;
    // CC RES
    case kXSecTwkDial_MaCCRES:
    case kXSecTwkDial_MvCCRES:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<3;i++){ if(ccres_shp_norm[i] == *it1) return false; }
      }
    break;
    case kXSecTwkDial_NormCCRES:
    case kXSecTwkDial_MaCCRESshape:
    case kXSecTwkDial_MvCCRESshape:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<2;i++){ if(ccres_shp[i]      == *it1) return false; }
      }
    break;
    // DIS
    case kXSecTwkDial_AhtBYshape:
    case kXSecTwkDial_BhtBYshape:
    case kXSecTwkDial_CV1uBYshape:
    case kXSecTwkDial_CV2uBYshape:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<4;i++){ if(dis_shp_norm[i]   == *it1) return false; }
      }
    break;
    case kXSecTwkDial_AhtBY:
    case kXSecTwkDial_BhtBY:
    case kXSecTwkDial_CV1uBY:
    case kXSecTwkDial_CV2uBY:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<4;i++){ if(dis_shp[i]        == *it1) return false; }
      }
    break;
    // NC RES
    case kXSecTwkDial_MaNCRES:
    case kXSecTwkDial_MvNCRES:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<3;i++){ if(ncres_shp_norm[i] == *it1) return false; }
      }
    break;
    case kXSecTwkDial_NormNCRES:
    case kXSecTwkDial_MaNCRESshape:
    case kXSecTwkDial_MvNCRESshape:
      for(it1=it0;it1 != lsyst.end();it1++){
        for(int i=0;i<2;i++){ if(ncres_shp[i]      == *it1) return false; }
      }
    break;
    default: // no conflicts
    break;
    }
  }
  return true;

}
//_________________________________________________________________________________
void AdoptWeightCalcs (vector<GSyst_t> lsyst, GReWeight & rw)
{
  //
  // Sets of systematics can be incompatible because they request different
  // reweighting modes. If one of these systematics is requested, adopt a
  // reweighting calculator and set it to the appropriate mode.
  //
  vector<GSyst_t>::iterator it;
  for(it=lsyst.begin();it != lsyst.end();it++)
  {
    switch(*it){
    // CC QE
    case kXSecTwkDial_MaCCQE:
    case kXSecTwkDial_NormCCQE:
    case kXSecTwkDial_MaCCQEshape:
      if ( ! rw.WghtCalc("xsec_ccqe") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_ccqe weight calc";
        rw.AdoptWghtCalc( "xsec_ccqe", new GReWeightNuXSecCCQE );
        GReWeightNuXSecCCQE * rwccqe =
          dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccqe"));
        if (*it == kXSecTwkDial_MaCCQE) {
               rwccqe->SetMode(GReWeightNuXSecCCQE::kModeMa); }
        else { rwccqe->SetMode(GReWeightNuXSecCCQE::kModeNormAndMaShape); }
      }
    break;
    case kXSecTwkDial_ZNormCCQE:
    case kXSecTwkDial_ZExpA1CCQE:
    case kXSecTwkDial_ZExpA2CCQE:
    case kXSecTwkDial_ZExpA3CCQE:
    case kXSecTwkDial_ZExpA4CCQE:
      if ( ! rw.WghtCalc("xsec_ccqe") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_ccqe weight calc";
        rw.AdoptWghtCalc( "xsec_ccqe", new GReWeightNuXSecCCQE );
        GReWeightNuXSecCCQE * rwccqe =
          dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccqe"));
        rwccqe->SetMode(GReWeightNuXSecCCQE::kModeZExp);
      }
    break;
    // CC Res
    case kXSecTwkDial_MaCCRES:
    case kXSecTwkDial_MvCCRES:
      if ( ! rw.WghtCalc("xsec_ccres") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_ccres weight calc";
        rw.AdoptWghtCalc( "xsec_ccres", new GReWeightNuXSecCCRES );
        GReWeightNuXSecCCRES * rwccres =
          dynamic_cast<GReWeightNuXSecCCRES *> (rw.WghtCalc("xsec_ccres"));
        rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);
      }
    break;
    case kXSecTwkDial_NormCCRES:
    case kXSecTwkDial_MaCCRESshape:
    case kXSecTwkDial_MvCCRESshape:
      if ( ! rw.WghtCalc("xsec_ccres") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_ccres weight calc";
        rw.AdoptWghtCalc( "xsec_ccres", new GReWeightNuXSecCCQE );
        GReWeightNuXSecCCQE * rwccres =
          dynamic_cast<GReWeightNuXSecCCQE *> (rw.WghtCalc("xsec_ccres"));
        rwccres->SetMode(GReWeightNuXSecCCRES::kModeNormAndMaMvShape);
      }
    break;
    // DIS
    case kXSecTwkDial_AhtBYshape:
    case kXSecTwkDial_BhtBYshape:
    case kXSecTwkDial_CV1uBYshape:
    case kXSecTwkDial_CV2uBYshape:
      if ( ! rw.WghtCalc("xsec_dis") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_dis weight calc";
        rw.AdoptWghtCalc( "xsec_dis", new GReWeightNuXSecDIS );
        GReWeightNuXSecDIS * rwdis =
          dynamic_cast<GReWeightNuXSecDIS *> (rw.WghtCalc("xsec_dis"));
        rwdis->SetMode(GReWeightNuXSecDIS::kModeABCV12uShape);
      }
    break;
    case kXSecTwkDial_AhtBY:
    case kXSecTwkDial_BhtBY:
    case kXSecTwkDial_CV1uBY:
    case kXSecTwkDial_CV2uBY:
      if ( ! rw.WghtCalc("xsec_dis") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_dis weight calc";
        rw.AdoptWghtCalc( "xsec_dis", new GReWeightNuXSecDIS );
        GReWeightNuXSecDIS * rwdis =
          dynamic_cast<GReWeightNuXSecDIS *> (rw.WghtCalc("xsec_dis"));
        rwdis->SetMode(GReWeightNuXSecDIS::kModeABCV12u);
      }
    break;
    // NC Res
    case kXSecTwkDial_MaNCRES:
    case kXSecTwkDial_MvNCRES:
      if ( ! rw.WghtCalc("xsec_ncres") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_ncres weight calc";
        rw.AdoptWghtCalc( "xsec_ncres", new GReWeightNuXSecNCRES );
        GReWeightNuXSecNCRES * rwncres =
          dynamic_cast<GReWeightNuXSecNCRES *> (rw.WghtCalc("xsec_ncres"));
        rwncres->SetMode(GReWeightNuXSecNCRES::kModeMaMv);
      }
    break;
    case kXSecTwkDial_NormNCRES:
    case kXSecTwkDial_MaNCRESshape:
    case kXSecTwkDial_MvNCRESshape:
      if ( ! rw.WghtCalc("xsec_ncres") ){
        LOG("grwghtcov", pNOTICE) << "Adopting xsec_ncres weight calc";
        rw.AdoptWghtCalc( "xsec_ncres", new GReWeightNuXSecNCRES );
        GReWeightNuXSecNCRES * rwncres =
          dynamic_cast<GReWeightNuXSecNCRES *> (rw.WghtCalc("xsec_ncres"));
        rwncres->SetMode(GReWeightNuXSecNCRES::kModeNormAndMaMvShape);
      }
    break;
    default: // no fine-tuning needed
    break;
    }
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("grwghtcov", pFATAL)
     << "\n\n"
     << "grwghtcov                    \n"
     << "     -f input_event_file     \n"
     << "     -c input_covariance_file\n"
     << "     -t num_twk              \n"
     << "     -s syst1[,syst2[,...]]  \n"
     << "     -v cval1[,cval2[,...]]  \n"
     << "    [-n n1[,n2]]             \n"
     << "    [-r run_key]             \n"
     << "    [-o output_weights_file]";
}
//_________________________________________________________________________________
