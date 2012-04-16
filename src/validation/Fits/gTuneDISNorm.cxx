//____________________________________________________________________________
/*!

\program gtune_dis_norm

\brief   Cross-section tuning utility: 
         Fits DIS scale factor to high energy neutrino and anti-neutrino cross 
         section data.

\syntax  gtune_dis_norm 
             [-g genie_inputs] [-d data_archive] [--emin E1] [--emax E2]

         where

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 06, 2008 

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TMinuit.h>
#include <TMath.h>

#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"
#include "Utils/VldTestInputs.h"

using std::vector;
using std::string;

using namespace genie;
using namespace genie::utils;
using namespace genie::utils::vld;

const int kNModes = 2;

const char * kDataSets[kNModes] = 
{
/* mode 0 : nu_mu+N CC inclusive */ 
"ANL_12FT,2;ANL_12FT,4;BEBC,0;BEBC,2;BEBC,5;BEBC,8;BNL_7FT,0;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,0;CHARM,4;FNAL_15FT,1;FNAL_15FT,2;Gargamelle,0;Gargamelle,10;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0",

/* mode 1 : nu_mu_bar+N CC inclusive */ 
"BEBC,1;BEBC,3;BEBC,6;BEBC,7;BNL_7FT,1;CCFR,3;CHARM,1;CHARM,5;FNAL_15FT,4;FNAL_15FT,5;Gargamelle,1;Gargamelle,11;Gargamelle,13;IHEP_ITEP,1;IHEP_ITEP,3;IHEP_JINR,1"
};

// func prototypes
void GetCommandLineArgs (int argc, char ** argv);
void Init               (void);
void DoTheFit           (void);
void FitFunc            (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
void Save               (string filename);
void ReadData           (void);

// command-line arguments
VldTestInputs  gOptGenieInputs;
string         gOptDataFilename = "";
double         gOptEmax;
double         gOptEmin;

// default data archive
string kDefDataFile = "data/validation/vA/xsec/integrated/nuXSec.root";  

// default energy range
const double kEmin = 10.0; // GeV
const double kEmax = 90.0; // GeV

// simple functor which uses the program inputs to return the GENIE cross-section (inclusive (anti)neutrino 
// cross-section for an isoscalar target) needed for fitting the specified datasets.
class GNuXSec
{
public:
  GNuXSec() 
  {
    fBuiltDefaults = false;
    for(int imode = 0; imode < kNModes; imode++) {
      fXSecIncl [imode] = 0;
      fXSecDIS  [imode] = 0;
    }
  }
 ~GNuXSec() 
  {
    for(int imode = 0; imode < kNModes; imode++) {
      if (fXSecIncl [imode]) delete fXSecIncl [imode];
      if (fXSecDIS  [imode]) delete fXSecDIS  [imode];
    }
  }

  double operator() (double Ev, int imode, double dis_norm)
  { 
    if(imode < 0 || imode >= kNModes) return 0;
    if(!fBuiltDefaults) return 0;
    double xsec_dis_nominal   = fXSecDIS [imode] -> Eval(Ev);
    double xsec_incl_nominal  = fXSecIncl[imode] -> Eval(Ev);
    double xsec_nodis_nominal = xsec_incl_nominal - xsec_dis_nominal;
    double xsec_incl_tweaked  = xsec_nodis_nominal + dis_norm * xsec_dis_nominal;
    return xsec_incl_tweaked;
  }

  void BuildDefaults(const VldTestInputs & inp)
  {
    TFile * genie_xsec_file = inp.XSecFile(0);
    if(!genie_xsec_file) return;

    //TChain * genie_event_tree = inp.EvtChain(0);
    //if(!genie_event_tree) return;

    const int n = 200;
    const double dE = (gOptEmax - gOptEmin)/(n-1);

    double E[n];    
    double xsec_incl[kNModes][n];
    double xsec_dis [kNModes][n];

    for(int i=0; i<n; i++) { E[i] =  gOptEmin + i * dE; }

    for(int imode = 0; imode < kNModes; imode++) {
        for(int i=0; i<n; i++) { 
           xsec_incl[imode][i] = 0;
           xsec_dis [imode][i] = 0;
        }
    }

    const int ngraphs = 4;
    string directory       [ngraphs] = { "nu_mu_n",  "nu_mu_H1", "nu_mu_bar_n", "nu_mu_bar_H1" };
    string graph_xsec_incl [ngraphs] = { "tot_cc_n", "tot_cc_p", "tot_cc_n",    "tot_cc_p"     };
    string graph_xsec_dis  [ngraphs] = { "dis_cc_n", "dis_cc_p", "dis_cc_n",    "dis_cc_p"     };
    int    mode            [ngraphs] = {  0,          0,          1,             1             };        

    for(int igr = 0; igr < ngraphs; igr++) {
      int imode = mode[igr];
      TDirectory * dir = (TDirectory *) genie_xsec_file->Get(directory[igr].c_str());
      if(!dir) return;
      TGraph * gr = 0;
      gr = (TGraph *) dir->Get(graph_xsec_incl[igr].c_str());
      if(!gr) return;
      for(int i=0; i<n; i++) { xsec_incl[imode][i] = 0.5 * gr->Eval(E[i]); }
      gr = (TGraph *) dir->Get(graph_xsec_dis[igr].c_str());
      if(!gr) return;
      for(int i=0; i<n; i++) { xsec_dis[imode][i]  = 0.5 * gr->Eval(E[i]); }
    }

    for(int imode = 0; imode < kNModes; imode++) {
      fXSecIncl [imode] = new TGraph(n,E,xsec_incl[imode]);
      fXSecDIS  [imode] = new TGraph(n,E,xsec_dis [imode]);
    }

    fBuiltDefaults = true;
  }

private:
  bool     fBuiltDefaults;
  TGraph * fXSecIncl [kNModes];
  TGraph * fXSecDIS  [kNModes];
};



// globals
TFile *  gNuXSecDataFile  = 0;
TTree *  gNuXSecDataTree  = 0;
GNuXSec  gNuXSecFunc;
vector<TGraphAsymmErrors *> gData[kNModes];

//vector<int>                 gEnabledDataSets;

const int kNMaxDataSets = 25; // max number of datasets in single plot

const int kDataPointStyle[kNMaxDataSets] = 
{ 
  20,      20,    20,       20,    20,
  21,      21,    21,       21,    21,
  24,      24,    24,       24,    24,
  25,      25,    25,       25,    25,
  29,      29,    29,       29,    29
};
const int kDataPointColor[kNMaxDataSets] = 
{
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kBlack,  kRed,  kGreen+1, kBlue, kMagenta+1
};


//____________________________________________________________________________
int main(int argc, char ** argv)
{
  LOG("gtune", pNOTICE) << "Running DISNorm tuning program...";

  GetCommandLineArgs (argc,argv);
  Init();
  DoTheFit();
  Save("dis_norm_tune.root");

  LOG("gtune", pNOTICE) << "Done!";

  return 0;
}
//____________________________________________________________________________
void Init(void)
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
  gNuXSecFunc.BuildDefaults(gOptGenieInputs);

  /*
  string datasets = GetArgument(argc, argv, "-d");
  vector<string> dsvec = str::Split(datasets,",");
  vector<string>::const_iterator it = dsvec.begin();
  for( ; it != dsvec.end(); ++it) { 
    gEnabledDataSets.push_back( atoi(it->c_str()) );
  }
  */

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
void DoTheFit(void)
{
  // Initialize MINUIT
  //
  const int np = 1;
  TMinuit * minuit = new TMinuit(np);

  double arglist[10];
  int ierrflg = 0;

  arglist[0] = 1;
  minuit->mnexcm("SET ERR",arglist,1,ierrflg);

  float        value [np] = {  1.00     };
  float        min   [np] = {  0.50     };
  float        max   [np] = {  2.00     };
  float        step  [np] = {  0.01     };
  const char * pname [np] = { "DISNorm" };

  for(int i=0; i<np; i++) {
    LOG("gtune", pNOTICE)
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

  double dis_norm = par[0];

  LOG("gtune", pDEBUG) 
      << "DIS cross-section normalization factor set to " << dis_norm;
 
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
            double xsec_model    = gNuXSecFunc(E,imode,dis_norm);
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

  LOG("gtune", pNOTICE) 
     << "Chisq (DIS norm = " <<  dis_norm << ") = " << chisq;
}
//____________________________________________________________________________
void Save(string filename)
{
// write-out fit results

  TFile out(filename.c_str(), "recreate");
  out.cd();

  // save fitted data-sets
  for(int imode=0; imode<kNModes; imode++) {
    unsigned int ndatasets = gData[imode].size();
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gData[imode][idataset];
       assert(data);
       data->Write(Form("dataset_%d_%d",imode,idataset));
    }
  }

  const int n = 100;
  const double dE = (gOptEmax - gOptEmin)/(n-1);

  double E[n];
  for(int i=0; i<n; i++) { E[i] =  gOptEmin + i * dE; }

  double xsec_nominal[kNModes][n];
  double xsec_bestfit[kNModes][n];
  for(int imode=0; imode<kNModes; imode++) {
     for(int i=0; i<n; i++) { 
       xsec_nominal[imode][i] = gNuXSecFunc(E[i],imode,1);
       xsec_bestfit[imode][i] = gNuXSecFunc(E[i],imode,2);
     }
  }
  for(int imode=0; imode<kNModes; imode++) {
    TGraph * gr_xsec_nominal = new TGraph(n,E,xsec_nominal[imode]);
    gr_xsec_nominal->Write(Form("mc_nominal_%d",imode));
    TGraph * gr_xsec_bestfit = new TGraph(n,E,xsec_bestfit[imode]);
    gr_xsec_bestfit->Write(Form("mc_bestfit_%d",imode));
  }


  // write-out fitted params / errors  
  // ...

  out.Close();
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
