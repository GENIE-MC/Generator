//______________________________________________________________________________________________
/*!

\program gtune_rijk

\brief   Cross-section tuning utility: 
         Fits parameters controlling the cross section in the transition region.

          - nu_mu N -> mu- X
          - nu_mu p -> mu- p pi+
          - nu_mu n -> mu- p pi0
          - nu_mu n -> mu- n pi+ xsec
          - nu_mu n -> mu- p pi+ pi- xsec
          - nu_mu p -> mu- p pi+ pi0 xsec
          - nu_mu p -> mu- n pi+ pi+ xsec

\syntax  gtune_rijk  [-g genie_inputs] [-d data_archive]

         where 

           -d Location of the neutrino cross-section data archive.
              By default, the program will look-up the one located in:
              $GENIE/data/validation/vA/xsec/integrated/

           -g An XML file with GENIE inputs.
              They are files with calculated cross-sections and event samples
              used for decomposing the inclusive cross-section to various
              exclusive cross-sections.
              For info on the XML file format see the GSimFiles class documentation.
                      
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Hugh Gallagher <gallag \at minos.phy.tufts.edu>
         Tufts University

\created June 06, 2008 

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//______________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TMinuit.h>
#include <TMath.h>

#include "Messenger/Messenger.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"
#include "Utils/GSimFiles.h"

using std::map;
using std::vector;
using std::string;
using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

// datasets included in the fit for each mode
const int kNDataSets = 7;

const char * kDataSets[kNModes] = 
{
/* numu CC inclusive */ 
"ANL_12FT,2;ANL_12FT,4;BEBC,0;BEBC,2;BEBC,5;BEBC,8;BNL_7FT,0;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,0;CHARM,4;FNAL_15FT,1;FNAL_15FT,2;Gargamelle,0;Gargamelle,10;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0",

/* numu p -> mu- p pi+ */ 
"ANL_12FT,0;ANL_12FT,5;ANL_12FT,8;BEBC,4;BEBC,9;BEBC,13;BNL_7FT,5;FNAL_15FT,0;Gargamelle,4;SKAT,4;SKAT,5",

/* v n -> l- p pi0 */ 
"ANL_12FT,6;ANL_12FT,9;BNL_7FT,6;SKAT,6",

/* v n -> l- n pi+ */ 
"ANL_12FT,7;ANL_12FT,10;BNL_7FT,7;SKAT,7",

/* v n -> l- p pi+ pi- */ 
"ANL_12FT,11;BNL_7FT,8",

/*  v p -> l- p pi+ pi0 */ 
"ANL_12FT,12",

/*  v p -> l- n pi+ pi+ */ 
"ANL_12FT,13"
};

// func prototypes
void GetCommandLineArgs (int argc, char ** argv);
void Init               (void);
void DoTheFit           (void);
void FitFunc            (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
void Save               (string filename);
void ReadData           (void);

// command-line arguments
GSimFiles gOptGenieInputs;
string    gOptDataFilename = "";

// default data archive
string kDefDataFile = "data/validation/vA/xsec/integrated/nuXSec.root";  

// default energy range
const double kEmin =  0.1; // GeV
const double kEmax = 20.0; // GeV

// Wcut parameter for the RES/DIS joining algorithm
const double kWcut  =  1.7;  // GeV, 

// simple cross-section functor using the program inputs to 
class NuXSec
{
public:
  NuXSec() 
  {
  }
 ~NuXSec() 
  {
  }

  double operator() (double Ev, int imode)
  { 
  }
};

// globals
TFile * gNuXSecDataFile  = 0;
TTree * gNuXSecDataTree  = 0;
NuXSec  gNuXSecFunc;
vector<TGraphAsymmErrors *> gData[kNModes];
//vector<int>  gEnabledDataSets;

///< fitted parameter (corresponds to `DIS-HMultWgt-vp-CC-m2' in UsersPhysicsConfig.xml)
double gFittedParam_vpCC_m2; 
///< fitted parameter (corresponds to `DIS-HMultWgt-vp-CC-m3' in UsersPhysicsConfig.xml)
double gFittedParam_vpCC_m3; 
///< fitted parameter (corresponds to `DIS-HMultWgt-vn-CC-m2' in UsersPhysicsConfig.xml)
double gFittedParam_vnCC_m2; 
///< fitted parameter (corresponds to `DIS-HMultWgt-vn-CC-m3' in UsersPhysicsConfig.xml)
double gFittedParam_vnCC_m3; 

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  LOG("gtune", pNOTICE) << "Running Rijk tuning program...";

  GetCommandLineArgs (argc,argv);
  Init();
  DoTheFit();
  Save("rijk_tune.root");

  LOG("gtune", pNOTICE) << "Done!";

  return 0;
}
//____________________________________________________________________________
void Init(int argc, char ** argv)
{
  LOG("gtune", pNOTICE) << "Initializing...";

  utils::style::SetDefaultStyle();

  // get TTree with neutrino scattering cross-section data
  if( ! utils::system::FileExists(gOptDataFilename) ) {
      LOG("gtune", pFATAL) 
         << "Can not find file: " << gOptDataFilename;
      gAbortingInErr = true;
      exit(1);
  }
  gNuXSecDataFile = new TFile(gOptDataFilename.c_str(),"read");  
  gNuXSecDataTree = (TTree *) gNuXSecDataFile->Get("nuxsnt");
  if(!gNuXSecDataTree) {
      LOG("gtune", pFATAL) 
         << "Can not find TTree `nuxsnt' in file: " << gOptDataFilename;
      gAbortingInErr = true;
      exit(1);
  }

  // read the specified data-sets from the input tree
  ReadData();

  // configure cross-section functor
  //  gNuXSecIncl.BuildDefaults(gOptGenieInputs);
}
//____________________________________________________________________________
void ReadData(void)
{
  LOG("gtune", pNOTICE) 
    << "Getting specified cross-section data from the data archive";

  const int buffer_size = 100;

  char   dataset  [buffer_size];
  char   citation [buffer_size];
  double E;
  double Emin;
  double Emax;
  double xsec;
  double xsec_err_p;
  double xsec_err_m;

  gNuXSecDataTree->SetBranchAddress ("dataset",    (void*)dataset  );
  gNuXSecDataTree->SetBranchAddress ("citation",   (void*)citation );
  gNuXSecDataTree->SetBranchAddress ("E",          &E              );
  gNuXSecDataTree->SetBranchAddress ("Emin",       &Emin           );
  gNuXSecDataTree->SetBranchAddress ("Emax",       &Emax           );
  gNuXSecDataTree->SetBranchAddress ("xsec",       &xsec           );
  gNuXSecDataTree->SetBranchAddress ("xsec_err_p", &xsec_err_p     );
  gNuXSecDataTree->SetBranchAddress ("xsec_err_m", &xsec_err_m     );

  // loop over modes
  for(int imode = 0; imode < kNModes; imode++) {

    LOG("gtune", pNOTICE) << "Getting datasets for mode ID: " << imode;

    string keys = kDataSets[imode];
    vector<string> keyv = utils::str::Split(keys,";");
    unsigned int ndatasets = keyv.size();

    // loop over datasets for current mode
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
      gNuXSecDataTree->Draw("E", Form("dataset==\"%s\"",keyv[idataset].c_str()), "goff");
      int npoints = gNuXSecDataTree->GetSelectedRows();
      double *  x    = new double[npoints];
      double *  dxl  = new double[npoints];
      double *  dxh  = new double[npoints];
      double *  y    = new double[npoints];
      double *  dyl  = new double[npoints];
      double *  dyh  = new double[npoints];
      string label = "";
      int ipoint=0;
      for(int i = 0; i < gNuXSecDataTree->GetEntries(); i++) {
        gNuXSecDataTree->GetEntry(i);
        if(strcmp(dataset,keyv[idataset].c_str()) == 0) {
          if(ipoint==0) {
            label = Form("%s [%s]", dataset, citation);
          }
          x   [ipoint] = E;
          dxl [ipoint] = (Emin > 0) ? TMath::Max(0., E-Emin) : 0.;
          dxh [ipoint] = (Emin > 0) ? TMath::Max(0., Emax-E) : 0.;
          y   [ipoint] = xsec;
          dyl [ipoint] = xsec_err_m;
          dyh [ipoint] = xsec_err_p;
          ipoint++;
        } 
      }//i
      TGraphAsymmErrors * gr = new TGraphAsymmErrors(npoints,x,y,dxl,dxh,dyl,dyh);
      int sty = kDataPointStyle[idataset];
      int col = kDataPointColor[idataset];
      utils::style::Format(gr,col, kSolid, 1, col, sty, 1.5);
      gr->SetTitle(label.c_str());
      gData[imode].push_back(gr);

      LOG("gtune", pNOTICE) 
        << "Done getting data-points for dataset: " << label;

    }//idataset
  }//imode
}
//____________________________________________________________________________

/*
double GetXSecGENIE(double E, int imode)
{
// get the GENIE cross section, for the current Rijk parameters, for the input
// mode (matching the fitted data ids; see header) at the input energy

  //
  // 0 -> total xsec
  // (still inefficient; integration on the fly - not based on precomputed data)
  //
  if(imode == 0) {
    TLorentzVector p4(0,0,E,E);
    double xsec_tot = 0.5 * (gNuMuPdrv->XSecSum(p4) + gNuMuNdrv->XSecSum(p4));
    xsec_tot /= (1E-38 * units::cm2);

    LOG("RijkFit", pDEBUG) 
       << "xsec(total; E = " << E << " GeV) = " << xsec_tot << " 1E-38 * cm2";

    return xsec_tot;
  }

  //
  // Exclusive reactions
  //
  // single pion channels:
  // 1 -> v p -> l- p pi+ 
  // 2 -> v n -> l- p pi0 
  // 3 -> v n -> l- n pi+ 
  //
  // multi pion channels:
  // 4 -> v n -> l- p pi+ pi- 
  // 5 -> v p -> l- p pi+ pi0 
  // 6 -> v p -> l- n pi+ pi+ 
  //

  // calculate the resonance contribution
  double xsec_res = gSplMapRES[imode]->Evaluate(E);

  // calculate the dis contribution for W <  Wcut (Rijk=1)
  double xsec_dis_b = gSplMapDIS_b[imode]->Evaluate(E); 

  // calculate the dis contribution for W >= Wcut (Rijk=1)
  double xsec_dis_a = gSplMapDIS_a[imode]->Evaluate(E);

  // get the appropriate R factor for the input mode (data set)
  double R = 0;
  if      (imode==1)             { R = gFittedParam_vpCC_m2; }
  else if (imode==2 || imode==3) { R = gFittedParam_vnCC_m2; }
  if      (imode==4)             { R = gFittedParam_vnCC_m3; }
  else if (imode==5 || imode==6) { R = gFittedParam_vpCC_m3; }

  // sum-up
  double xsec = xsec_res + R * xsec_dis_b + xsec_dis_a;

  LOG("RijkFit", pDEBUG) 
   << "xsec(mode: "<< imode<< ", E = "<< E<< " GeV) = "<< xsec<< " 1E-38*cm2";
  LOG("RijkFit", pDEBUG) 
   << " R=" << R << ", xsec_res = " << xsec_res 
   << ", xsec_dis [b,a] = " << xsec_dis_b << ", " << xsec_dis_a;

  return xsec;
}
*/
//____________________________________________________________________________
void DoTheFit(void)
{
  // Initialize MINUIT
  //
  const int np = 4;
  TMinuit * minuit = new TMinuit(np);

  double arglist[10];
  int ierrflg = 0;

  arglist[0] = 1;
  minuit->mnexcm("SET ERR",arglist,1,ierrflg);

  float        value [np] = { 0.10,          1.00,           0.30,           1.00           };
  float        min   [np] = { 0.00,          0.00,           0.00,           0.00           };
  float        max   [np] = { 2.00,          2.00,           2.00,           2.00           };
  float        step  [np] = { 0.05,          0.05,           0.05,           0.05           };
  const char * pname [np] = {"R(vp/CC/m2)",  "R(vp/CC/m3)",  "R(vn/CC/m2)",  "R(vn/CC/m3)"  };

  for(int i=0; i<np; i++) {
    LOG("RijkFit",pDEBUG)
        << "** Setting fit param " << i
          << " (" << pname[i] << ") value = " << value[i]
             << ", range = [" << min[i] << ", " << max[i] <<"]";

    minuit->mnparm(
       i, pname[i], value[i], step[i], min[i], max[i], ierrflg);
  }

  minuit->SetFCN(FitFunc);

  // MINUIT minimization step
  ierrflg    = 0;
  arglist[0] = 500;
  arglist[1] = 1.;
  minuit->mnexcm("MIGRAD",arglist,2,ierrflg);

  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit->mnprin(3,amin);
}
//____________________________________________________________________________
void FitFunc (
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
// the MINUIT fit function

  LOG("RijkFit", pNOTICE) << "Setting: R(vp/CC/m2) = " << par[0];
  LOG("RijkFit", pNOTICE) << "Setting: R(vp/CC/m3) = " << par[1];
  LOG("RijkFit", pNOTICE) << "Setting: R(vn/CC/m2) = " << par[2];
  LOG("RijkFit", pNOTICE) << "Setting: R(vn/CC/m3) = " << par[3];

  gFittedParam_vpCC_m2 = par[0];
  gFittedParam_vpCC_m3 = par[1];
  gFittedParam_vnCC_m2 = par[2];
  gFittedParam_vnCC_m3 = par[3];

  double chisq = 0;

  // loop over all modes included in the fit
  for(int imode = 0; imode < kNModes; imode++) {

    /*    // include current data set?
    vector<int>::const_iterator it =
        find(gEnabledDataSets.begin(), gEnabledDataSets.end(), imode);
    bool skip = (it==gEnabledDataSets.end());
    if(skip) continue;
    */
    // loop over graphs in current data-set (one graph per experiment/publication in this data set)
    unsigned int ndatasets = gData[imode].size();
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {

       TGraphAsymmErrors * data = gData[imode][idataset];
       assert(data);

       LOG("gtune", pNOTICE) 
	 << " Mode: " << imode << ", Dataset : " << idataset+1 << "/" << ndatasets 
         << " [" <<  data->GetTitle() << "]";

       // loop over data-points in current graph
       int np = data->GetN();
       for (int ip = 0; ip < np; ip++) {

          double E = data->GetX()[ip];
          bool in_fit_range = (E >= gOptEmin && E <= gOptEmax);

          if(in_fit_range) {
            double xsec_data     = data->GetY()[ip];    
            double xsec_data_err = data->GetErrorY(ip); 
            double xsec_model    = gNuXSecFunc(E,imode);
            double delta = (xsec_data>0) ? (xsec_data - xsec_model) / xsec_data_err : 0.;
            chisq += delta*delta;

            LOG("gtune", pNOTICE)
               << " > pnt " << ip+1 << "/" << np 
               << " @ E = " << E << " GeV : "
               << "Data = "  << xsec_data << " +/- " << xsec_data_err << " x1E-38 cm^2/GeV/nucleon, "
               << "Model = " << xsec_model << " x1E-38 cm^2/GeV/nucleon "
               << " ==> Running chisq = " << chisq;

	  } else {

            LOG("gtune", pNOTICE)
               << " > pnt " << ip+1 << "/" << np 
               << " @ E = " << E << " GeV : ** not in fit range **";
          }

       } // graph points
    } // graph
  } // data set

  f = chisq;

  // LOG("gtune", pNOTICE) 
  //   << "Chisq (DIS norm = " <<  dis_norm << ") = " << chisq;
}
//____________________________________________________________________________
void Save(string /*filename*/)
{

}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gtune", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // get data archive
  if(parser.OptionExists('d')){
     string filename = parser.ArgAsString('d');
     gOptDataFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataFile;
        gOptDataFilename = filename;
     } else { 
        LOG("gtune", pFATAL) 
          << "\n Please make sure that $GENIE is defined, or use the -d option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }
  // get GENIE inputs
  if( parser.OptionExists('g') ) {
     string inputs = parser.ArgAsString('g');
     bool ok = gOptGenieInputs.LoadFromFile(inputs);
     if(!ok) {
        LOG("gtune", pFATAL) << "Could not read: " << inputs;
        gAbortingInErr = true;
        exit(1);
     }
  } else {
    LOG("gtune", pFATAL) 
       << "Please specify GENIE inputs using the -g option";
    gAbortingInErr = true;
    exit(1);
  }

  gOptEmax = kEmax;
  gOptEmin = kEmin;
}
//____________________________________________________________________________
