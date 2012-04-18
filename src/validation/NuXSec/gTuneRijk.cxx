//______________________________________________________________________________________________
/*!

\program gtune_rijk

\brief   Cross-section tuning utility: 
         Fits GENIE parameters RvpCC1pi, RvnCC1pi, RvpCC2pi and RvnCC2pi 
         controlling the cross-section in the transition region.

         The cross-section model is fit simultaneously to integrated cross-section data on:
          - nu_mu N -> mu- X
          - nu_mu p -> mu- p pi+
          - nu_mu n -> mu- n pi+ 
          - nu_mu n -> mu- p pi0
          - nu_mu p -> mu- n pi+ pi+ 
          - nu_mu p -> mu- p pi+ pi0 
          - nu_mu n -> mu- p pi+ pi- 

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
#include <TCanvas.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TVirtualFitter.h>
#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "Registry/Registry.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"
#include "Utils/GSimFiles.h"
#include "validation/NuXSec/NuXSecData.h"

using std::map;
using std::vector;
using std::string;
using std::ostringstream;

using namespace genie;
using namespace genie::mc_vs_data;

const int kNFitParams = 4;

typedef enum {
  kUndef    = -1,
  kRvpCC1pi =  0,
  kRvnCC1pi =  1,
  kRvpCC2pi =  2,
  kRvnCC2pi =  3
} RijkParam_t;

// default data archive
string kDefDataFile = "data/validation/vA/xsec/integrated/nuXSec.root";  

// energy and W range
const int    kNE   =  50; 
const double kEmin =  0.1;  // GeV
const double kEmax = 10.0;  // GeV
const int    kNW   =  50; 
const double kWmin =  0.1;  // GeV
const double kWmax =  5.0;  // GeV

// Wcut parameter for the RES/DIS joining algorithm
const double kWcut  = 1.7;  // GeV

// Number of modes included in the fit
const int kNModes = 7;
const int kNExclusiveModes = 6;

//
// datasets included in the fit for each mode
//
const char * kDataSets[kNModes] = 
{
// mode 0: nu_mu CC inclusive 
"ANL_12FT,2;ANL_12FT,4;BEBC,0;BEBC,2;BEBC,5;BEBC,8;BNL_7FT,0;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,0;CHARM,4;FNAL_15FT,1;FNAL_15FT,2;Gargamelle,0;Gargamelle,10;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0",

// mode 1: nu_mu p -> mu- p pi+ 
"ANL_12FT,0;ANL_12FT,5;ANL_12FT,8;BEBC,4;BEBC,9;BEBC,13;BNL_7FT,5;FNAL_15FT,0;Gargamelle,4;SKAT,4;SKAT,5",

// mode 2: nu_mu n -> mu- n pi+ 
"ANL_12FT,7;ANL_12FT,10;BNL_7FT,7;SKAT,7",

// mode 3: nu_mu n -> mu- p pi0 
"ANL_12FT,6;ANL_12FT,9;BNL_7FT,6;SKAT,6",

// mode 4: nu_mu p -> mu- n pi+ pi+ 
"ANL_12FT,13",

// mode 5: nu_mu p -> mu- p pi+ pi0 
"ANL_12FT,12",

// mode 6: nu_mu n -> mu- p pi+ pi- 
"ANL_12FT,11;BNL_7FT,8"
};

//
// GENIE cross-section functor used in fit.
//
class XSecFunc
{
public:
   XSecFunc() {}
  ~XSecFunc() {}

  void Init(GSimFiles & genie_inputs, double RNom[kNFitParams])
  {
    for(int iparam = 0; iparam < kNFitParams; iparam++) {
      fRNominal[iparam] = RNom[iparam];
    }

    TFile * genie_xsec_file = genie_inputs.XSecFile(0);
    assert(genie_xsec_file);
    TChain * genie_event_tree = genie_inputs.EvtChain(0);
    assert(genie_event_tree);

    const double dE = (kEmax - kEmin)/(kNE-1);
    double E[kNE];   
    for(int iE=0; iE<kNE; iE++) { E[iE] =  kEmin + iE * dE; } 

    double xsec_numu_n[kNE];
    double xsec_numu_p[kNE];
    double xsec_numu_N[kNE];
    for(int iE=0; iE<kNE; iE++) { 
      xsec_numu_n[iE] = 0.;
      xsec_numu_p[iE] = 0.;
      xsec_numu_N[iE] = 0.;
    } 

    TDirectory * dir = 0;
    TGraph * gr = 0;

    dir = (TDirectory *) genie_xsec_file->Get("nu_mu_n");
    assert(dir);
    gr = (TGraph *) dir->Get("tot_cc_n");
    assert(gr);
    for(int iE = 0; iE < kNE; iE++) { 
      xsec_numu_n[iE] = gr->Eval(E[iE]);
    }

    dir = (TDirectory *) genie_xsec_file->Get("nu_mu_H1");
    assert(dir);
    gr = (TGraph *) dir->Get("tot_cc_p");
    assert(gr);
    for(int iE = 0; iE < kNE; iE++) { 
      xsec_numu_p[iE] = gr->Eval(E[iE]);
    }

    for(int iE = 0; iE < kNE; iE++) { 
      xsec_numu_N[iE] = 0.5 * (xsec_numu_n[iE] + xsec_numu_p[iE]);
    }

    TGraph * gr_xsec_incl_numu_n = new TGraph(kNE, E, xsec_numu_n);
    TGraph * gr_xsec_incl_numu_p = new TGraph(kNE, E, xsec_numu_p);
    TGraph * gr_xsec_incl_numu_N = new TGraph(kNE, E, xsec_numu_N);

    string evt_selection_res [kNExclusiveModes] = { 
        "cc&&res&&neu==14&&Z==1&&A==1&&nfpim==0&&nfpi0==0&&nfpip==1&&nfp==1&&nfn==0", // nu_mu p -> mu- p pi+
        "cc&&res&&neu==14&&Z==0&&A==1&&nfpim==0&&nfpi0==0&&nfpip==1&&nfp==0&&nfn==1", // nu_mu n -> mu- n pi+ 
	"cc&&res&&neu==14&&Z==0&&A==1&&nfpim==0&&nfpi0==1&&nfpip==0&&nfp==1&&nfn==0", // nu_mu n -> mu- p pi0
        "cc&&res&&neu==14&&Z==1&&A==1&&nfpim==0&&nfpi0==0&&nfpip==2&&nfp==0&&nfn==1", // nu_mu p -> mu- n pi+ pi+ 
        "cc&&res&&neu==14&&Z==1&&A==1&&nfpim==0&&nfpi0==1&&nfpip==1&&nfp==1&&nfn==0", // nu_mu p -> mu- p pi+ pi0 
        "cc&&res&&neu==14&&Z==0&&A==1&&nfpim==1&&nfpi0==0&&nfpip==1&&nfp==1&&nfn==0"  // nu_mu n -> mu- p pi+ pi- 
    };
    string evt_selection_nonres [kNExclusiveModes] = { 
        "cc&&!res&&neu==14&&Z==1&&A==1&&nfpim==0&&nfpi0==0&&nfpip==1&&nfp==1&&nfn==0", // nu_mu p -> mu- p pi+
        "cc&&!res&&neu==14&&Z==0&&A==1&&nfpim==0&&nfpi0==0&&nfpip==1&&nfp==0&&nfn==1", // nu_mu n -> mu- n pi+ 
	"cc&&!res&&neu==14&&Z==0&&A==1&&nfpim==0&&nfpi0==1&&nfpip==0&&nfp==1&&nfn==0", // nu_mu n -> mu- p pi0
        "cc&&!res&&neu==14&&Z==1&&A==1&&nfpim==0&&nfpi0==0&&nfpip==2&&nfp==0&&nfn==1", // nu_mu p -> mu- n pi+ pi+ 
        "cc&&!res&&neu==14&&Z==1&&A==1&&nfpim==0&&nfpi0==1&&nfpip==1&&nfp==1&&nfn==0", // nu_mu p -> mu- p pi+ pi0 
        "cc&&!res&&neu==14&&Z==0&&A==1&&nfpim==1&&nfpi0==0&&nfpip==1&&nfp==1&&nfn==0"  // nu_mu n -> mu- p pi+ pi- 
    };
    string evt_selection_incl [kNExclusiveModes] = { 
        "cc&&neu==14&&Z==1&&A==1", 
        "cc&&neu==14&&Z==0&&A==1", 
	"cc&&neu==14&&Z==0&&A==1", 
        "cc&&neu==14&&Z==1&&A==1", 
        "cc&&neu==14&&Z==1&&A==1", 
        "cc&&neu==14&&Z==0&&A==1"  
    };
    TGraph * gr_xsec_incl [kNExclusiveModes] = { 
        gr_xsec_incl_numu_p,
        gr_xsec_incl_numu_n,
        gr_xsec_incl_numu_n,
        gr_xsec_incl_numu_p,
        gr_xsec_incl_numu_p,
        gr_xsec_incl_numu_n
    };

    // Fit modes 1-6 are exclusive CC 1pi and 2pi cross-sections
    for(int imode = 1; imode < kNModes; imode++) {
        LOG("gtune", pNOTICE) << "Building nominal dxsec/dW for mode: " << imode;
        int iexclmode = imode - 1;
        TH2D * hincl       = new TH2D("hincl", "", kNE,kEmin,kEmax,kNW,kWmin,kWmax);
        TH2D * hexclres    = new TH2D(Form("dxsecdW_res_mode%d",   imode),"",kNE,kEmin,kEmax,kNW,kWmin,kWmax);
        TH2D * hexclnonres = new TH2D(Form("dxsecdW_nonres_mode%d",imode),"",kNE,kEmin,kEmax,kNW,kWmin,kWmax);
	genie_event_tree->Draw("Ws:Ev>>hincl",       evt_selection_incl  [iexclmode].c_str(), "goff");
	LOG("gtune", pNOTICE) << "Selection: " <<  evt_selection_incl  [iexclmode] << ", entries = " << genie_event_tree->GetSelectedRows();
        genie_event_tree->Draw("Ws:Ev>>hexclres",    evt_selection_res   [iexclmode].c_str(), "goff");
	LOG("gtune", pNOTICE) << "Selection: " <<  evt_selection_res  [iexclmode] << ", entries = " << genie_event_tree->GetSelectedRows();
        genie_event_tree->Draw("Ws:Ev>>hexclnonres", evt_selection_nonres[iexclmode].c_str(), "goff");
	LOG("gtune", pNOTICE) << "Selection: " <<  evt_selection_nonres  [iexclmode] << ", entries = " << genie_event_tree->GetSelectedRows();
        hexclres    -> Divide (hincl);
        hexclnonres -> Divide (hincl);
        for(int iE = 0; iE<kNE; iE++) { 
	  int iEbin = iE+1;
          double energy    = hincl->GetXaxis()->GetBinCenter(iEbin);
          double xsec_incl = gr_xsec_incl[iexclmode]->Eval(energy);
          for(int iW = 0; iW<kNW; iW++) { 
	     int iWbin = iW+1;
             double evt_frac_res    = hexclres    -> GetBinContent(iEbin,iWbin);
             double evt_frac_nonres = hexclnonres -> GetBinContent(iEbin,iWbin);
             double xsec_res        = xsec_incl * evt_frac_res;
             double xsec_nonres     = xsec_incl * evt_frac_nonres;
             hexclres    -> SetBinContent(iEbin,iWbin,xsec_res);
             hexclnonres -> SetBinContent(iEbin,iWbin,xsec_nonres);
          }//iW
        }//iE
        fXSecRes   [iexclmode] = hexclres;
	fXSecNonRes[iexclmode] = hexclnonres;
    }//imode

  }//init()

  double operator() (int imode, double E, double R[kNFitParams])
  {
    double wght[kNFitParams];
    for(int ip = 0; ip < kNFitParams; ip++) { wght[ip] = R[ip]/fRNominal[ip]; }

    // CC inclusive
    if(imode == 0) {
      double sum_xsec_excl = 0;
      for(int jmode = 1; jmode < kNModes; jmode++) {      
         sum_xsec_excl += (*this)(jmode,E,R);
      }
      return sum_xsec_excl;
    }

    // One of the exclusive CC 1pi and 2pi channels
    if(imode >= 1 && imode < kNModes) {
      int wght_idx[kNExclusiveModes] = { kRvpCC1pi, kRvpCC1pi, kRvnCC1pi, kRvpCC2pi, kRvpCC2pi, kRvnCC2pi };
      int iexclmode = imode - 1;
      // Integrate resonance and non-resonance dxsec/dW(E,W) along W for the input energy.
      // Apply R weighting factors where appropriate
      double xsec = 0;
      int e_bin = fXSecRes[iexclmode]->GetXaxis()->FindBin(E);
      for(int w_bin = 1; w_bin <= fXSecRes[iexclmode]->GetYaxis()->GetNbins(); w_bin++) {
 	 double W = fXSecRes[iexclmode]->GetYaxis()->GetBinCenter(w_bin);
         double dxsec_res    = fXSecRes   [iexclmode]->GetBinContent(e_bin,w_bin);
         double dxsec_nonres = fXSecNonRes[iexclmode]->GetBinContent(e_bin,w_bin);
         double wght_nonres  = (W < kWcut) ? wght_idx[iexclmode] : 1.;
         double dxsec = dxsec_res + wght_nonres * dxsec_nonres;
         xsec += dxsec;
      }//w
      return xsec;
    }

    return 0;
  }

private:

  double fRNominal[kNFitParams]; // nominal non-resonance background params

  TH2D * fXSecRes   [kNExclusiveModes];  // nominal resonance     dxsec/dW = f(Ev,W;mode)
  TH2D * fXSecNonRes[kNExclusiveModes];  // nominal non-resonance dxsec/dW = f(Ev,W;mode)
  TH2D * fXSecOthExcl;                   // nominal dxsec/dW = f(Ev,W) for all other exclusive channels not in fit
  TH2D * fXSecIncl;                      // nominal inclusive dxsec/dW = f(Ev,W)
};

// command-line arguments
string gOptDataFilename  = ""; // -d
string gOptGenieFileList = ""; // -g

// nominal value of fit parameters, as used in GENIE simulation,
// and best-fit value
double gRNominal[kNFitParams];
double gRBestFit[kNFitParams];

// data and xsec function used in fit
vector<TGraphAsymmErrors *> gXSecData[kNModes];
XSecFunc                    gXSecFunc;

//
// func prototypes
//
void GetCommandLineArgs (int argc, char ** argv);
void Init               (void);
void DoTheFit           (void);
void FitFunc            (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
void Save               (string filename);

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
void Init(void)
{
  LOG("gtune", pNOTICE) << "Initializing...";

  utils::style::SetDefaultStyle();

  // Read data archive & retrieve/store specified datasets
  NuXSecData data_reader;
  bool ok = data_reader.OpenArchive(gOptDataFilename);
  if(!ok) {
      LOG("gtune", pFATAL)
         << "Could not open the neutrino cross-section data archive";
      gAbortingInErr = true;
      exit(1);
  }
  for(int imode = 0; imode < kNModes; imode++) {
    string datasets = kDataSets[imode];
    gXSecData[imode] = data_reader.Retrieve(datasets);
  }

  // Read GENIE inputs
  GSimFiles genie_inputs;
  ok = genie_inputs.LoadFromFile(gOptGenieFileList);
  if(!ok) {
     LOG("gtune", pFATAL)
         << "Could not read GENIE inputs specified in XML file: "
         << gOptGenieFileList;
     gAbortingInErr = true;
     exit(1);
  }

  // Get nominal value of fit parameter, as used in GENIE simulation.
  // Obvious caveat here: Assuming that the user runs this program using the exact same
  // configuration used when the GENIE inputs were produced.
  Registry * params = AlgConfigPool::Instance()->GlobalParameterList();
  gRNominal[kRvpCC1pi] = params->GetDouble("DIS-HMultWgt-vp-CC-m2");
  gRNominal[kRvnCC1pi] = params->GetDouble("DIS-HMultWgt-vn-CC-m2");
  gRNominal[kRvpCC2pi] = params->GetDouble("DIS-HMultWgt-vp-CC-m3");
  gRNominal[kRvnCC2pi] = params->GetDouble("DIS-HMultWgt-vn-CC-m3");
  for(int ip = 0; ip < kNFitParams; ip++) {
    gRBestFit[ip] = gRNominal[ip]; // init
  }
  /*
  LOG("gtune", pNOTICE)
    << "\nNominal R(vp;CC;1pi) value used in simulation: " << gRvpCC1piNominal
    << "\nNominal R(vn;CC;1pi) value used in simulation: " << gRvnCC1piNominal
    << "\nNominal R(vp;CC;2pi) value used in simulation: " << gRvpCC2piNominal
    << "\nNominal R(vn;CC;2pi) value used in simulation: " << gRvnCC2piNominal;
  */

  // Configure cross-section functor
  gXSecFunc.Init(genie_inputs, gRNominal);
}
//____________________________________________________________________________
void DoTheFit(void)
{
  float        init [kNFitParams] = { 0.10,       0.30,       1.00,       1.00       };
  float        min  [kNFitParams] = { 0.00,       0.00,       0.00,       0.00       };
  float        max  [kNFitParams] = { 2.00,       2.00,       2.00,       2.00       };
  float        step [kNFitParams] = { 0.05,       0.05,       0.05,       0.05       };
  const char * name [kNFitParams] = { "vpCC1pi",  "RvnCC1pi", "RvpCC2pi", "RvnCC2pi" };

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * fitter = TVirtualFitter::Fitter(0,kNFitParams);
      
  double arglist[100];
  arglist[0] = -1;
  fitter->ExecuteCommand("SET PRINT",arglist,1);

  for(int i = 0; i < kNFitParams; i++) {
    LOG("gtune", pNOTICE)
        << "** Setting fit param " << i
          << " (" << name[i] << ") initial value = " << init[i]
             << ", range = [" << min[i] << ", " << max[i] <<"]";

    fitter->SetParameter(i, name[i], init[i], step[i], min[i], max[i]);
  }

  fitter->SetFCN(FitFunc);

  // MINUIT minimization step
  arglist[0] = 500;   // num of func calls
  arglist[1] = 0.01;  // tolerance
  fitter->ExecuteCommand("MIGRAD",arglist,2);

  // Get fit status code
  double ha(0.);      //minuit dummy vars
  double edm, errdef; //minuit dummy vars
  int nvpar, nparx;   //minuit dummy vars
  int status_code = fitter->GetStats(ha,edm,errdef,nvpar,nparx);
  LOG("gtune", pNOTICE)
    << "Minuit status code = " << status_code;

  // Get best-fit values
  for(int ip = 0; ip < kNFitParams; ip++) {
     gRBestFit[ip] = fitter->GetParameter(ip);
  }

  // Print results
  double amin = 0;
  fitter->PrintResults(3,amin);
}
//____________________________________________________________________________
void FitFunc (
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
// the MINUIT fit function

  double chisq = 0;

  // loop over all modes included in the fit
  for(int imode = 0; imode < kNModes; imode++) {

    // loop over graphs in current data-set (one graph per experiment/publication in this data set)
    unsigned int ndatasets = gXSecData[imode].size();
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {

       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       assert(data);

       LOG("gtune", pNOTICE) 
	 << " Mode: " << imode << ", Dataset : " << idataset+1 << "/" << ndatasets 
         << " [" <<  data->GetTitle() << "]";

       // loop over data-points in current graph
       int np = data->GetN();
       for (int ip = 0; ip < np; ip++) {

          double E = data->GetX()[ip];
          bool in_fit_range = (E >= kEmin && E <= kEmax);

          if(in_fit_range) {
            double xsec_data     = data->GetY()[ip];    
            double xsec_data_err = data->GetErrorY(ip); 
            double xsec_model    = gXSecFunc(E,imode,par);
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
}
//____________________________________________________________________________
void Save(string filename)
{
// write-out fit results

  TFile out(filename.c_str(), "recreate");
  out.cd();

  // generate xsec prediction for nominal and best-fit param values
  TGraph * gr_xsec_bestfit[kNModes];
  TGraph * gr_xsec_nominal[kNModes];
  const int n = 100;
  const double dE = (kEmax - kEmin)/(n-1);
  double E[n];
  for(int i=0; i<n; i++) { E[i] = kEmin + i * dE; }
  double xsec_nominal[kNModes][n];
  double xsec_bestfit[kNModes][n];
  for(int imode=0; imode<kNModes; imode++) {
     for(int i=0; i<n; i++) { 
       xsec_nominal[imode][i] = gXSecFunc(imode,E[i],gRNominal);
       xsec_bestfit[imode][i] = gXSecFunc(imode,E[i],gRBestFit);
     }
  }
  for(int imode=0; imode<kNModes; imode++) {
    gr_xsec_bestfit[imode] = new TGraph(n,E,xsec_bestfit[imode]);
    gr_xsec_bestfit[imode]->SetLineStyle(kSolid);
    gr_xsec_bestfit[imode]->SetLineWidth(2);
    gr_xsec_nominal[imode] = new TGraph(n,E,xsec_nominal[imode]);
    gr_xsec_nominal[imode]->SetLineStyle(kDashed);
    gr_xsec_nominal[imode]->SetLineWidth(2);
  }

  // plot nominal and best-fit MC & datasets
  TCanvas * c[kNModes];
  for(int imode=0; imode<kNModes; imode++) {
    c[imode] = new TCanvas(Form("c_%d",imode), "", 10,10,400,400);
    TGraph * bestfit = gr_xsec_bestfit[imode];
    TGraph * nominal = gr_xsec_nominal[imode];
    bestfit->Draw("al");
    nominal->Draw("l");
    unsigned int ndatasets = gXSecData[imode].size();
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       data->Draw("P");
    }
    c[imode]->Update();
    c[imode]->Write();
  }

  // save fitted data-sets
  for(int imode=0; imode<kNModes; imode++) {
    unsigned int ndatasets = gXSecData[imode].size();
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       assert(data);
       data->Write(Form("dataset_%d_%d",imode,idataset));
    }
  }
  // save nominal and best-fit MC
  for(int imode=0; imode<kNModes; imode++) {
    gr_xsec_bestfit[imode]->Write(Form("mc_bestfit_%d",imode));
    gr_xsec_nominal[imode]->Write(Form("mc_nominal_%d",imode));
  }

  // write-out fitted param values and covariance matrix
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
     gOptGenieFileList = parser.ArgAsString('g');
  } else {
    LOG("gtune", pFATAL)
       << "Please specify GENIE inputs using the -g option";
    gAbortingInErr = true;
    exit(1);
  }
}
//____________________________________________________________________________
