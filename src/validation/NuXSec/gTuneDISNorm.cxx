//____________________________________________________________________________
/*!

\program gtune_dis_norm

\brief   Cross-section tuning utility: 
         Fits DIS scale factor to high-energy neutrino and anti-neutrino 
         cross-section data.

         The nominal cross-section is extracted from the input files and 
         then reweighted. For the DIS norm fit we don't need to use 
         event-by-event reweighting. Only the GENIE cross-section file
         needs to be input.

\syntax  gtune_dis_norm 
             [-g genie_inputs] [-d data_archive] [--emin E1] [--emax E2]

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
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 06, 2008 

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
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
#include <TVirtualFitter.h>
#include <TPostScript.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLatex.h>

#include "Algorithm/AlgConfigPool.h"
#include "Messenger/Messenger.h"
#include "Registry/Registry.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"
#include "Utils/GSimFiles.h"

#include "validation/NuXSec/NuXSecData.h"

using std::vector;
using std::string;

using namespace genie;
using namespace genie::mc_vs_data;

//
// Constants
//

// default data archive
string kDefDataFile = "data/validation/vA/xsec/integrated/nuXSec.root";  

// default energy range
const double kEmin =  20.0; // GeV
const double kEmax = 150.0; // GeV

// Number of modes included in the fit
const int kNModes = 2;

//
// Datasets used in fit for each mode
//
const char * kDataSets[kNModes] = 
{
// mode 0 : nu_mu+N CC inclusive 
"ANL_12FT,4;BEBC,8;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,4;FNAL_15FT,2;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0;MINOS,0",
/*
"ANL_12FT,2;ANL_12FT,4;BEBC,0;BEBC,2;BEBC,5;BEBC,8;BNL_7FT,0;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,0;CHARM,4;FNAL_15FT,1;FNAL_15FT,2;Gargamelle,0;Gargamelle,10;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0;MINOS,0",
*/

// mode 1 : nu_mu_bar+N CC inclusive 
"BEBC,7;BNL_7FT,1;CCFR,3;CHARM,5;FNAL_15FT,5;Gargamelle,13;IHEP_ITEP,3;IHEP_JINR,1;MINOS,1",
/*
"BEBC,1;BEBC,3;BEBC,6;BEBC,7;BNL_7FT,1;CCFR,3;CHARM,1;CHARM,5;FNAL_15FT,4;FNAL_15FT,5;Gargamelle,1;Gargamelle,11;Gargamelle,13;IHEP_ITEP,1;IHEP_ITEP,3;IHEP_JINR,1;MINOS,1",
*/
};

const char * kLabel[kNModes] = 
{
  "#nu_{#mu} CC inclusive",
  "#bar{#nu_{#mu}} CC inclusive"
};

//
// GENIE cross-section functor used in fit.
//
class InclXSecFunc
{
public:
   InclXSecFunc() {}
  ~InclXSecFunc() {}

  void Init(GSimFiles & genie_inputs, double Emin, double Emax, double dis_norm_nominal) 
   {
      const int ngraphs = 4;
      string directory       [ngraphs] = { "nu_mu_n",  "nu_mu_H1", "nu_mu_bar_n", "nu_mu_bar_H1" };
      string graph_xsec_incl [ngraphs] = { "tot_cc_n", "tot_cc_p", "tot_cc_n",    "tot_cc_p"     };
      string graph_xsec_dis  [ngraphs] = { "dis_cc_n", "dis_cc_p", "dis_cc_n",    "dis_cc_p"     };
      int    mode            [ngraphs] = {  0,          0,          1,             1             };        
      TFile * genie_xsec_file = genie_inputs.XSecFile(0);
      assert(genie_xsec_file);
      const int n = 100;
      const double dE = (Emax - Emin)/(n-1);
      double E[n];   
      double xsec_incl[kNModes][n];
      double xsec_dis [kNModes][n];
      for(int i=0; i<n; i++) { E[i] =  Emin + i * dE; }	
      for(int imode = 0; imode < kNModes; imode++) {
         for(int i=0; i<n; i++) { 
            xsec_incl[imode][i] = 0;
            xsec_dis [imode][i] = 0;
         }
      }	
      for(int igr = 0; igr < ngraphs; igr++) {
         int imode = mode[igr];
         TDirectory * dir = (TDirectory *) genie_xsec_file->Get(directory[igr].c_str());
         assert(dir);
         TGraph * gr_incl = (TGraph *) dir->Get(graph_xsec_incl[igr].c_str());
         assert(gr_incl);
         for(int i=0; i<n; i++) { xsec_incl[imode][i] += (0.5 * gr_incl->Eval(E[i])); }
         TGraph * gr_dis  = (TGraph *) dir->Get(graph_xsec_dis[igr].c_str());
         assert(gr_dis);
         for(int i=0; i<n; i++) { xsec_dis[imode][i]  += (0.5 * gr_dis->Eval(E[i]));  }
      }
      for(int imode = 0; imode < kNModes; imode++) {
         fXSecIncl [imode] = new TGraph(n,E,xsec_incl[imode]);
         fXSecDIS  [imode] = new TGraph(n,E,xsec_dis [imode]);
      }

      fDISNormNominal = dis_norm_nominal;
      assert(fDISNormNominal > 0.);
   }

  double operator() (int imode, double E, double dis_norm)
  {
    if(E <= 0) return 0;
    if(imode < 0 || imode >= kNModes) return 0;
    if(!fXSecDIS[imode] || !fXSecIncl[imode]) return 0;
    double wght_dis           = dis_norm / fDISNormNominal;
    double xsec_dis_nominal   = fXSecDIS [imode] -> Eval(E);
    double xsec_incl_nominal  = fXSecIncl[imode] -> Eval(E);
    double xsec_nodis_nominal = xsec_incl_nominal - xsec_dis_nominal;
    double xsec_incl_tweaked  = xsec_nodis_nominal + wght_dis * xsec_dis_nominal;
    LOG("gtune", pDEBUG)
      << "xsec_incl (E = " << E << " GeV, norm = " << dis_norm << ") : " 
      << xsec_incl_nominal << " --> " << xsec_incl_tweaked << " 1E-38 cm^2/GeV/nucleon"; 
    assert(xsec_incl_tweaked > 0);
    xsec_incl_tweaked /= E; // actually return sigma/E
    return xsec_incl_tweaked;
  }    

private:
  double   fDISNormNominal;
  TGraph * fXSecDIS  [kNModes];
  TGraph * fXSecIncl [kNModes];
};

//
// Func prototypes
//
void   GetCommandLineArgs (int argc, char ** argv);
void   Init               (void);
void   DoTheFit           (void);
void   FitFunc            (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
double Chisq              (double dis_norm);
void   Save               (string filename);

//
// Globals
//

// command-line arguments
string gOptDataFilename  = ""; // -d
string gOptGenieFileList = ""; // -g
double gOptEmin = -1;          //
double gOptEmax = -1;          //

// nominal value of fit parameter, as used in GENIE simulation,
// and best-fit value
double gDISNormNominal = -1;
double gDISNormBestFit = -1;

// data and xsec function used in fit
vector<TGraphAsymmErrors *> gXSecData[kNModes];
InclXSecFunc                gXSecFunc;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  LOG("gtune", pNOTICE) << "Running DISNorm tuning program...";

  GetCommandLineArgs (argc,argv);
  Init();
  DoTheFit();
  Save("dis_norm_tune");

  LOG("gtune", pNOTICE) << "Done!";

  return 0;
}
//____________________________________________________________________________
void Init(void)
{
  LOG("gtune", pNOTICE) << "Initializing...";

  // Set GENIE style
  utils::style::SetDefaultStyle();

  // Read data archive & retrieve/store specified datasets
  NuXSecData data_reader;
  bool ok = data_reader.OpenArchive(gOptDataFilename);
  if(!ok) {
      LOG("gvldtest", pFATAL) 
         << "Could not open the neutrino cross-section data archive";
      gAbortingInErr = true;
      exit(1);
  }
  for(int imode = 0; imode < kNModes; imode++) {
    string datasets = kDataSets[imode];
    gXSecData[imode] = data_reader.Retrieve(datasets,gOptEmin,gOptEmax,true);
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
  gDISNormNominal = params->GetDouble("DIS-XSecScale");
  gDISNormBestFit = gDISNormNominal; // init
  LOG("gtune", pNOTICE) 
    << "Nominal DIS norm. value used in simulation: " << gDISNormNominal;

  // Configure cross-section functor
  gXSecFunc.Init(genie_inputs,gOptEmin,gOptEmax,gDISNormNominal);

}
//____________________________________________________________________________
void DoTheFit(void)
{
  const int nparams = 1;

  float        nominal [nparams] = {  1.00     };
  float        min     [nparams] = {  0.50     };
  float        max     [nparams] = {  2.00     };
  float        step    [nparams] = {  0.001    };
  const char * name    [nparams] = { "DISNorm" };

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * fitter = TVirtualFitter::Fitter(0,nparams);

  double arglist[100];
  arglist[0] = -1;
  fitter->ExecuteCommand("SET PRINT",arglist,1);

  for(int i=0; i<nparams; i++) {
    LOG("gtune", pNOTICE)
        << "** Setting fit param " << i
          << " (" << name[i] << ") nominal = " << nominal[i]
             << ", range = [" << min[i] << ", " << max[i] <<"]";

    fitter->SetParameter(i, name[i], nominal[i], step[i], min[i], max[i]);
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

  // Get best-fit value
  gDISNormBestFit = fitter->GetParameter(0);

  // Print results
  double amin = 0;
  fitter->PrintResults(3,amin);
}
//____________________________________________________________________________
void FitFunc (
   Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
// MINUIT fit function with signature expected by TVirtualFitter::SetFCN()

  double dis_norm = par[0];  
  f = Chisq(dis_norm);
}
//____________________________________________________________________________
double Chisq(double dis_norm)
{
  LOG("gtune", pDEBUG) 
      << "DIS cross-section normalization factor set to " << dis_norm;
 
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
          bool in_fit_range = (E >= gOptEmin && E <= gOptEmax);

          if(in_fit_range) {
            double xsec_data     = data->GetY()[ip];    
            double xsec_data_err = data->GetErrorY(ip); 
            double xsec_model    = gXSecFunc(imode,E,dis_norm);
            double delta = (xsec_data>0) ? (xsec_data - xsec_model) / xsec_data_err : 0.;
            chisq += delta*delta;

            LOG("gtune", pINFO)
               << " > pnt " << ip+1 << "/" << np 
               << " @ E = " << E << " GeV : "
               << "Data = "  << xsec_data << " +/- " << xsec_data_err << " x1E-38 cm^2/GeV/nucleon, "
               << "Model = " << xsec_model << " x1E-38 cm^2/GeV/nucleon "
               << " ==> Running chisq = " << chisq;

	  } else {

            LOG("gtune", pINFO)
               << " > pnt " << ip+1 << "/" << np 
               << " @ E = " << E << " GeV : ** not in fit range **";
          }

       } // graph points
    } // graph
  } // data set

  LOG("gtune", pNOTICE) 
     << "Chisq (DIS norm = " <<  dis_norm << ") = " << chisq;

  return chisq;
}
//____________________________________________________________________________
void Save(string filename)
{
// write-out fit results

  // Set GENIE style
  utils::style::SetDefaultStyle();

  // Generate xsec prediction for nominal and best-fit param values
  TGraph * gr_xsec_bestfit[kNModes];
  TGraph * gr_xsec_nominal[kNModes];
  const int n = 100;
  const double dE = (gOptEmax - gOptEmin)/(n-1);
  double E[n];
  for(int i=0; i<n; i++) { E[i] =  gOptEmin + i * dE; }
  double xsec_nominal[kNModes][n];
  double xsec_bestfit[kNModes][n];
  for(int imode=0; imode<kNModes; imode++) {
     for(int i=0; i<n; i++) { 
       xsec_nominal[imode][i] = gXSecFunc(imode,E[i],gDISNormNominal);
       xsec_bestfit[imode][i] = gXSecFunc(imode,E[i],gDISNormBestFit);
     }
  }
  for(int imode=0; imode<kNModes; imode++) {
    gr_xsec_bestfit[imode] = new TGraph(n,E,xsec_bestfit[imode]);
    gr_xsec_nominal[imode] = new TGraph(n,E,xsec_nominal[imode]);
  }
  // Generate chisq vs dis_norm plot
  const int npv = 100;
  const double disnorm_min  = 0.75;
  const double disnorm_max  = 1.25;
  const double disnorm_step = (disnorm_max - disnorm_min)/(npv-1);
  double disnorm[npv];
  double chisq  [npv];
  double chisq_min = 9999999999;
  for(int ipv=0; ipv<npv; ipv++) {
    disnorm[ipv] = disnorm_min + ipv * disnorm_step;
    chisq  [ipv] = Chisq(disnorm[ipv]);
    chisq_min = TMath::Min(chisq_min,chisq[ipv]);
  }
  for(int ipv=0; ipv<npv; ipv++) {
    chisq  [ipv] -= chisq_min;
  }
  TGraph * gr_chisq = new TGraph(npv,disnorm,chisq);

  // Save data, nominal and best-fit MC and chisq plots in a ps file
  TCanvas * c = new TCanvas("c","",20,20,500,650);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  TPostScript * ps = new TPostScript(Form("%s.ps",filename.c_str()), 111);
  ps->NewPage();
  c->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE dis_norm tune");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.Draw();
  c->Update();
  for(int imode=0; imode<kNModes; imode++) {
     ps->NewPage();
     c->Clear();
     c->Divide(2,1);
     c->GetPad(1)->SetPad("mplots_pad","",0.01,0.35,0.99,0.99);
     c->GetPad(2)->SetPad("legend_pad","",0.01,0.01,0.99,0.34);
     c->GetPad(1)->SetFillColor(0);
     c->GetPad(1)->SetBorderMode(0);
     c->GetPad(2)->SetFillColor(0);
     c->GetPad(2)->SetBorderMode(0);
     c->GetPad(1)->cd();
     c->GetPad(1)->SetBorderMode(0);
     TH1F * hframe = 0;
     double xmin =  9999999;
     double xmax = -9999999;
     double ymin =  9999999;
     double ymax = -9999999;
     unsigned int ndatasets = gXSecData[imode].size();
     for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       if(!data) continue;
       xmin  = TMath::Min(xmin, (data->GetX())[TMath::LocMin(data->GetN(),data->GetX())]);
       xmax  = TMath::Max(xmax, (data->GetX())[TMath::LocMax(data->GetN(),data->GetX())]);
       ymin  = TMath::Min(ymin, (data->GetY())[TMath::LocMin(data->GetN(),data->GetY())]);
       ymax  = TMath::Max(ymax, (data->GetY())[TMath::LocMax(data->GetN(),data->GetY())]);
     }
     hframe = (TH1F*) c->GetPad(1)->DrawFrame(0.8*xmin, 0.4*ymin, 1.2*xmax, 1.2*ymax);
     hframe->GetXaxis()->SetTitle("E_{#nu} (GeV)");
     hframe->GetYaxis()->SetTitle("#sigma_{#nu}/E_{#nu} (1E-38 cm^{2}/GeV^{2})");
     hframe->Draw();
     for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       if(!data) continue;
       data->Draw("P");
     }
     TGraph * bestfit = gr_xsec_bestfit[imode];
     TGraph * nominal = gr_xsec_nominal[imode];
     bestfit->SetLineStyle(kSolid);
     bestfit->SetLineWidth(2);
     nominal->SetLineColor(kRed);
     nominal->SetLineStyle(kDashed);
     nominal->SetLineWidth(2);
     bestfit->Draw("l");
     nominal->Draw("l");
     c->GetPad(1)->SetLogx(1);
     c->GetPad(1)->Update();
     c->GetPad(2)->cd();
     TLegend * legend = new TLegend(0.01, 0.01, 0.99, 0.99);
     legend->SetLineStyle(0);
     legend->SetFillColor(0);
     legend->SetTextSize(0.06);
     legend->SetHeader(kLabel[imode]);
     legend->AddEntry(bestfit, "GENIE best-fit", "L");
     legend->AddEntry(nominal, "GENIE nominal",  "L");
     for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       if(!data) continue;
       legend->AddEntry(data, data->GetTitle(), "LP");
     }
     legend->SetTextSize(0.05);
     legend->Draw();
     c->GetPad(2)->Update();
     c->Update();
  }
  {
     ps->NewPage();
     c->Clear();
     TH1F * hframe = 0;
     double xmin =  9999999;
     double xmax = -9999999;
     double ymin =  9999999;
     double ymax = -9999999;
     xmin  = TMath::Min(xmin, (gr_chisq->GetX())[TMath::LocMin(gr_chisq->GetN(),gr_chisq->GetX())]);
     xmax  = TMath::Max(xmax, (gr_chisq->GetX())[TMath::LocMax(gr_chisq->GetN(),gr_chisq->GetX())]);
     ymin  = TMath::Min(ymin, (gr_chisq->GetY())[TMath::LocMin(gr_chisq->GetN(),gr_chisq->GetY())]);
     ymax  = TMath::Max(ymax, (gr_chisq->GetY())[TMath::LocMax(gr_chisq->GetN(),gr_chisq->GetY())]);
     hframe = (TH1F*) c->DrawFrame(0.8*xmin, 0.4*ymin, 1.2*xmax, 1.2*ymax);
     hframe->GetXaxis()->SetTitle("DIS normalization");
     hframe->GetYaxis()->SetTitle("#Delta#chi^{2}");
     hframe->Draw("same");
     TLatex bestfit_info;
     bestfit_info.SetTextColor(kRed);
     bestfit_info.SetTextSize(0.03);
     bestfit_info.DrawLatex(0.9*xmin,1.1*ymax,Form("Best-fit DIS normalization scale = %f",gDISNormBestFit));
     gr_chisq->Draw("l");
     ps->Close();
  }

  // Save data, nominal and best-fit MC and chisq graphs in a root file
  TFile out(Form("%s.root",filename.c_str()), "recreate");
  out.cd();
  for(int imode=0; imode<kNModes; imode++) {
    unsigned int ndatasets = gXSecData[imode].size();
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       assert(data);
       data->Write(Form("dataset_%d_%d",imode,idataset));
    }
  }
  for(int imode=0; imode<kNModes; imode++) {
    gr_xsec_bestfit[imode]->Write(Form("mc_bestfit_%d",imode));
    gr_xsec_nominal[imode]->Write(Form("mc_nominal_%d",imode));
  }
  gr_chisq->Write("chisq");
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

  gOptEmax = kEmax;
  gOptEmin = kEmin;
}
//____________________________________________________________________________

