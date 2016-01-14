//______________________________________________________________________________________________
/*!

\program gtune_rijk

\brief   Cross-section tuning utility: 
         Fits GENIE parameters (MaRES,) RvpCC1pi, RvnCC1pi, RvpCC2pi and RvnCC2pi 
         controlling the cross-section in the Delta and in the transition regions.

         The cross-section model is fit simultaneously to integrated cross-section data on:
          - nu_mu N -> mu- X
          - nu_mu p -> mu- p pi+
          - nu_mu n -> mu- n pi+ 
          - nu_mu n -> mu- p pi0
          - nu_mu p -> mu- n pi+ pi+ 
          - nu_mu p -> mu- p pi+ pi0 
          - nu_mu n -> mu- p pi+ pi- 

\syntax  gtune_rijk  [-g genie_inputs] [-d data_archive] [--fit-params ]

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

         Hugh Gallagher <gallag \at minos.phy.tufts.edu>
         Tufts University

\created June 06, 2008 

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//______________________________________________________________________________________________

#include <cassert>
#include <cstdlib>
#include <sstream>
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

#include "Algorithm/AlgConfigPool.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGCodes.h"
#include "Registry/Registry.h"
#include "ReWeight/GReWeight.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"
#include "Utils/GSimFiles.h"
#include "validation/NuXSec/NuXSecData.h"

using std::vector;
using std::string;
using std::ostringstream;

using namespace genie;
using namespace genie::rew;
using namespace genie::mc_vs_data;

//
// Constants
//

// number and enumeration of fit params
const int kNFitParams = 4;

const int kRvpCC1pi = 0;
const int kRvnCC1pi = 1;
const int kRvpCC2pi = 2;
const int kRvnCC2pi = 3;

// default data archive
string kDefDataFile = "data/validation/vA/xsec/integrated/nuXSec.root";  

// energy and W range
const int    kNE   =  50; 
const double kEmin =  0.4;  // GeV
const double kEmax = 10.0;  // GeV
const int    kNW   =  500; 
const double kWmin =  0.0;  // GeV
const double kWmax =  5.0;  // GeV

// Wcut parameter for the RES/DIS joining algorithm
const double kWcut  = 1.7;  // GeV

// Number of modes included in the fit
const int kNModes = 7;
const int kNExclusiveModes = 6;

// datasets included in the fit for each mode
const char * kDataSets[kNModes] = 
{
// mode 0: nu_mu CC inclusive 
"ANL_12FT,2;ANL_12FT,4;BEBC,0;BEBC,2;BEBC,5;BEBC,8;BNL_7FT,0;BNL_7FT,4;CCFR,2;CCFRR,0;CHARM,0;CHARM,4;FNAL_15FT,1;FNAL_15FT,2;Gargamelle,0;Gargamelle,10;Gargamelle,12;IHEP_ITEP,0;IHEP_ITEP,2;IHEP_JINR,0;SKAT,0;MINOS,0",
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

// (Roo)TeX label for each mode
const char * kLabel[kNModes] = 
{
  "#nu_{#mu} N CC inclusive",                           // mode 0
  "#nu_{#mu} p #rightarrow #mu^{-} p #pi^{+}" ,         // mode 1
  "#nu_{#mu} n #rightarrow #mu^{-} n #pi^{+}",          // mode 2
  "#nu_{#mu} n #rightarrow #mu^{-} p #pi^{0}",          // mode 3
  "#nu_{#mu} p #rightarrow #mu^{-} n #pi^{+} #pi^{+}",  // mode 4
  "#nu_{#mu} p #rightarrow #mu^{-} p #pi^{+} #pi^{0}",  // mode 5
  "#nu_{#mu} n #rightarrow #mu^{-} p #pi^{+} #pi^{-}"   // mode 6
};

//..............................................................................................
//
// GENIE cross-section functor used in fit.
//
class XSecFunc
{
public:
   XSecFunc() {}
  ~XSecFunc() {}

  void Init(GSimFiles * genie_inputs)
  {
     fInLogE = true;

     const int n = 100;
     double xmin = (fInLogE) ? TMath::Log10(kEmin) : kEmin;
     double xmax = (fInLogE) ? TMath::Log10(kEmax) : kEmax;

     //  Get input sample and cross-section calcs
     assert(genie_inputs);
     int imodel = 0;
     fXSecFile = genie_inputs->XSecFile(imodel);
     assert(fXSecFile);
     fEventSamples = genie_inputs->EvtChain(imodel);
     assert(fEventSamples);

     // Set event tree branch address
     fEventSamples->SetBranchStatus("gmcrec", 1);
     fNev = fEventSamples->GetEntries();
     if (fNev<0) {
        LOG("gvldtest", pERROR) << "Number of events = 0";
        exit(1);
     }
     LOG("gvldtest", pNOTICE)
       << "Found " << fNev << " entries in the event tree";
     fMCRecord = 0;
     fEventSamples->SetBranchAddress("gmcrec", &fMCRecord);
     if (!fMCRecord) {
        LOG("gvldtest", pERROR) << "Null MC record";
        exit(1);
     }

     //
     TDirectory * xsec_dir_vp =
         (TDirectory *) fXSecFile->Get("nu_mu_H1");
     if(!xsec_dir_vp) {
        LOG("gvldtest", pERROR)
           << "Can't find cross-section directory: `nu_mu_H1'";
        exit(1);
     }
     TGraph * cc_incl_xsec_vp = (TGraph*) xsec_dir_vp->Get("tot_cc");
     if(!cc_incl_xsec_vp) {
        LOG("gvldtest", pERROR)
           << "Can't find inclusive CC cross-section calculation";
        exit(1);
     }
     TDirectory * xsec_dir_vn =
         (TDirectory *) fXSecFile->Get("nu_mu_n");
     if(!xsec_dir_vn) {
        LOG("gvldtest", pERROR)
           << "Can't find cross-section directory: `nu_mu_n'";
        exit(1);
     }
     TGraph * cc_incl_xsec_vn = (TGraph*) xsec_dir_vn->Get("tot_cc");
     if(!cc_incl_xsec_vn) {
        LOG("gvldtest", pERROR)
           << "Can't find inclusive CC cross-section calculation";
        exit(1);
     }

     fEvDistr[0] = 0; // not used
     for(int imode = 1; imode < kNModes; imode++) {
       fEvDistr[imode] = new TH1D( Form("hEv_%d",imode), "", n, xmin, xmax);
     }
     fEvDistrCCvp = new TH1D( "hEv_CCvp", "", n, xmin, xmax);
     fEvDistrCCvn = new TH1D( "hEv_CCvn", "", n, xmin, xmax);

     // Add weight calculation engines
     fRew.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
     fRew.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );

     // Fine-tune weight calculation engines
     GReWeightNuXSecCCRES * rwccres =
         dynamic_cast<GReWeightNuXSecCCRES *> (fRew.WghtCalc("xsec_ccres"));
     rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);

     GSystSet & systlist = fRew.Systematics();
   //systlist.Init(kXSecTwkDial_MaCCRES);
     systlist.Init(kXSecTwkDial_RvpCC1pi);
     systlist.Init(kXSecTwkDial_RvpCC2pi);
     systlist.Init(kXSecTwkDial_RvnCC1pi);
     systlist.Init(kXSecTwkDial_RvnCC2pi);
  }//init()

  void Update(double params[kNFitParams])
  {
     GSystSet & systlist = fRew.Systematics();
     systlist.Set(kXSecTwkDial_RvpCC1pi, params[kRvpCC1pi]); 
     systlist.Set(kXSecTwkDial_RvpCC2pi, params[kRvpCC2pi]);
     systlist.Set(kXSecTwkDial_RvnCC1pi, params[kRvnCC1pi]);
     systlist.Set(kXSecTwkDial_RvnCC2pi, params[kRvnCC2pi]);
     fRew.Reconfigure();

     for(Long64_t iev = 0; iev < fNev; iev++) {
       fEventSamples->GetEntry(iev);
       EventRecord & event = *(fMCRecord->event);
       Interaction * in = event.Summary();
       if(!in->ProcInfo().IsWeakCC()) continue;
       GHepParticle * neutrino = event.Probe();
       if(neutrino->Pdg() != kPdgNuMu) continue;
       GHepParticle * target = event.Particle(1);
       int tgtpdg = target->Pdg();
       double E = neutrino->P4()->E();
       double x = (fInLogE) ? TMath::Log10(E) : E;

       int npip = 0;
       int npi0 = 0;
       int npim = 0;
       int np   = 0;
       int nn   = 0;
       TObjArrayIter piter(&event);
       GHepParticle * p = 0;
       while ((p = (GHepParticle *) piter.Next())) {
          if(p->Status() != kIStStableFinalState) continue;
          if(p->Pdg() == kPdgPiP    ) npip++;
          if(p->Pdg() == kPdgPi0    ) npi0++;
          if(p->Pdg() == kPdgPiM    ) npim++;
          if(p->Pdg() == kPdgProton ) np++;
          if(p->Pdg() == kPdgNeutron) nn++;
       }//p

       int imode = -1;
       if      (tgtpdg == kPdgProton  && npip == 1 && npi0 == 0 && npim == 0 && np == 1 && nn == 0) imode = 1;
       else if (tgtpdg == kPdgNeutron && npip == 1 && npi0 == 0 && npim == 0 && np == 0 && nn == 1) imode = 2;
       else if (tgtpdg == kPdgProton  && npip == 0 && npi0 == 1 && npim == 0 && np == 1 && nn == 0) imode = 3;
       else if (tgtpdg == kPdgNeutron && npip == 2 && npi0 == 0 && npim == 0 && np == 0 && nn == 1) imode = 4;
       else if (tgtpdg == kPdgProton  && npip == 1 && npi0 == 1 && npim == 0 && np == 1 && nn == 0) imode = 5;
       else if (tgtpdg == kPdgProton  && npip == 1 && npi0 == 0 && npim == 1 && np == 1 && nn == 0) imode = 6;
       
       if(imode != -1) {
	 fEvDistr[imode]->Fill(x);
       }

       if      (tgtpdg == kPdgProton ) fEvDistrCCvp -> Fill(x);
       else if (tgtpdg == kPdgNeutron) fEvDistrCCvn -> Fill(x);

     }//iev

     fEvDistr[1]->Divide(fEvDistrCCvp);
     fEvDistr[2]->Divide(fEvDistrCCvn);
     fEvDistr[3]->Divide(fEvDistrCCvp);
     fEvDistr[4]->Divide(fEvDistrCCvn);
     fEvDistr[5]->Divide(fEvDistrCCvp);
     fEvDistr[6]->Divide(fEvDistrCCvp);
  }

  double operator() (int imode, double E)
  {
     return fCurrentXSec[imode]->Eval(E);
  }

private:

  TGraph *           fNominalXSec [kNModes];        ///< xsec = f(E;mode)
  TGraph *           fCurrentXSec [kNModes];        ///< xsec = f(E;mode)
  TFile *            fXSecFile;                     ///<
  TChain *           fEventSamples;                 ///<
  NtpMCEventRecord * fMCRecord;                     ///<
  Long64_t           fNev;                          ///<
  GReWeight          fRew;                          ///<
  TH1D *             fEvDistr [kNModes];            ///<
  TH1D *             fEvDistrCCvp;                  ///<
  TH1D *             fEvDistrCCvn;                  ///<
  bool               fInLogE;
};
//..............................................................................................

//
// Globals
//

// command-line arguments
string gOptDataFilename  = ""; // -d
string gOptGenieFileList = ""; // -g

// nominal value of fit parameters, as used in GENIE simulation, and best-fit value
double gRNominal[kNFitParams];
double gRBestFit[kNFitParams];

// data and xsec function used in fit
vector<TGraphAsymmErrors *> gXSecData[kNModes];
XSecFunc                    gXSecFunc;

//
// Function prototypes
//

void   GetCommandLineArgs (int argc, char ** argv);
void   Init               (void);
void   DoTheFit           (void);
void   FitFunc            (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
double Chisq              (double R[kNFitParams]);
void   Save               (string filename);

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  LOG("gtune", pNOTICE) << "Running Rijk tuning program...";

  GetCommandLineArgs (argc,argv);
  Init();
  DoTheFit();
  Save("rijk_tune");

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
    gXSecData[imode] = data_reader.Retrieve(datasets,kEmin,kEmax);
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
  LOG("gtune", pNOTICE)
    << "\nNominal R(vp;CC;1pi) value used in simulation: " << gRNominal[kRvpCC1pi]
    << "\nNominal R(vn;CC;1pi) value used in simulation: " << gRNominal[kRvnCC1pi]
    << "\nNominal R(vp;CC;2pi) value used in simulation: " << gRNominal[kRvpCC2pi]
    << "\nNominal R(vn;CC;2pi) value used in simulation: " << gRNominal[kRvnCC2pi];

  // Configure cross-section functor
  gXSecFunc.Init(&genie_inputs);
}
//____________________________________________________________________________
void DoTheFit(void)
{
  float        init [kNFitParams] = { 0.10,       0.30,       1.00,       1.00       };
  float        min  [kNFitParams] = { 0.00,       0.00,       0.00,       0.00       };
  float        max  [kNFitParams] = { 2.00,       2.00,       2.00,       2.00       };
  float        step [kNFitParams] = { 0.05,       0.05,       0.05,       0.05       };
  const char * name [kNFitParams] = { "RvpCC1pi", "RvnCC1pi", "RvpCC2pi", "RvnCC2pi" };

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
// MINUIT fit function with signature expected by TVirtualFitter::SetFCN()

  f = Chisq(par);
}
//____________________________________________________________________________
double Chisq(double * par)
{
  LOG("gtune", pNOTICE)
    << "Setting: "
    << "\n R(vp;CC;1pi) = " << par[kRvpCC1pi]
    << "\n R(vn;CC;1pi) = " << par[kRvnCC1pi]
    << "\n R(vp;CC;2pi) = " << par[kRvpCC2pi]
    << "\n R(vn;CC;2pi) = " << par[kRvnCC2pi];

  // update cross-section func with new params
  gXSecFunc.Update(par);

  double chisq = 0;

  // loop over all modes included in the fit
  for(int imode = 0; imode < kNModes; imode++) {

    // loop over graphs in current data-set 
    // (one graph per experiment/publication in this data set)
    unsigned int ndatasets = gXSecData[imode].size();
    for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {

       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       assert(data);

       LOG("gtune", pDEBUG) 
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
            double xsec_model    = gXSecFunc(E,imode);
            double delta = (xsec_data>0) ? (xsec_data - xsec_model) / xsec_data_err : 0.;
            chisq += delta*delta;

            LOG("gtune", pDEBUG)
               << " > pnt " << ip+1 << "/" << np 
               << " @ E = " << E << " GeV : "
               << "Data = "  << xsec_data << " +/- " << xsec_data_err << " x1E-38 cm^2/GeV/nucleon, "
               << "Model = " << xsec_model << " x1E-38 cm^2/GeV/nucleon "
               << " ==> Running chisq = " << chisq;

	  } else {

            LOG("gtune", pDEBUG)
               << " > pnt " << ip+1 << "/" << np 
               << " @ E = " << E << " GeV : ** not in fit range **";
          }

       } // graph points
    } // graph
  } // data set

  LOG("gtune", pNOTICE) << "Chisq = " << chisq;

  return chisq;
}
//____________________________________________________________________________
void Save(string filename)
{
// Write-out fit results

  // Set GENIE style
  utils::style::SetDefaultStyle();

  // Generate xsec prediction for nominal and best-fit param values
  TGraph * gr_xsec_bestfit[kNModes];
  TGraph * gr_xsec_nominal[kNModes];
  const int n = 100;
  const double dE = (kEmax - kEmin)/(n-1);
  double E[n];
  for(int i=0; i<n; i++) { E[i] = kEmin + i * dE; }
  double xsec_nominal[kNModes][n];
  gXSecFunc.Update(gRNominal);
  for(int imode=0; imode<kNModes; imode++) {
     for(int i=0; i<n; i++) { 
       xsec_nominal[imode][i] = gXSecFunc(imode,E[i]);
     }
  }
  double xsec_bestfit[kNModes][n];
  gXSecFunc.Update(gRBestFit);
  for(int imode=0; imode<kNModes; imode++) {
     for(int i=0; i<n; i++) { 
       xsec_bestfit[imode][i] = gXSecFunc(imode,E[i]);
     }
  }

  for(int imode=0; imode<kNModes; imode++) {
    gr_xsec_bestfit[imode] = new TGraph(n,E,xsec_bestfit[imode]);
    gr_xsec_nominal[imode] = new TGraph(n,E,xsec_nominal[imode]);
    gr_xsec_bestfit[imode]->SetLineStyle(kSolid);
    gr_xsec_bestfit[imode]->SetLineWidth(2);
    gr_xsec_nominal[imode]->SetLineStyle(kDashed);
    gr_xsec_nominal[imode]->SetLineWidth(2);
    gr_xsec_nominal[imode]->SetLineColor(kRed);
  }

  // Save data, nominal and best-fit MC and chisq plots in a ps file
  TCanvas * c = new TCanvas("c","",20,20,500,650);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  TPostScript * ps = new TPostScript(Form("%s.ps",filename.c_str()), 111);
  ps->NewPage();
  c->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE R_{ijk} tune");
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

     //
     // Draw frame
     //
     TH1F * hframe = 0;
     double xmin =  9999999;
     double xmax = -9999999;
     double ymin =  9999999;
     double ymax = -9999999;
     // Get x,y range from data
     unsigned int ndatasets = gXSecData[imode].size();
     for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       if(!data) continue;
       xmin  = TMath::Min(xmin, (data->GetX())[TMath::LocMin(data->GetN(),data->GetX())]);
       xmax  = TMath::Max(xmax, (data->GetX())[TMath::LocMax(data->GetN(),data->GetX())]);
       ymin  = TMath::Min(ymin, (data->GetY())[TMath::LocMin(data->GetN(),data->GetY())]);
       ymax  = TMath::Max(ymax, (data->GetY())[TMath::LocMax(data->GetN(),data->GetY())]);
     }
     // Also take into account the y range from the model predictions in the x range of the data
     TGraph * bestfit = gr_xsec_bestfit[imode];
     TGraph * nominal = gr_xsec_nominal[imode];
     if(!bestfit) continue;     
     if(!nominal) continue;     
     for(int k = 0; k < bestfit->GetN(); k++) {
         double x = (bestfit->GetX())[k];
         if(x < xmin || x > xmax) continue;
         ymin = TMath::Min(ymin, (bestfit->GetY())[k]);
         ymax = TMath::Max(ymax, (bestfit->GetY())[k]);
     }
     for(int k = 0; k < nominal->GetN(); k++) {
         double x = (nominal->GetX())[k];
         if(x < xmin || x > xmax) continue;
         ymin = TMath::Min(ymin, (nominal->GetY())[k]);
         ymax = TMath::Max(ymax, (nominal->GetY())[k]);
     }
     hframe = (TH1F*) c->GetPad(1)->DrawFrame(0.8*xmin, 0.4*ymin, 1.2*xmax, 1.2*ymax);
     hframe->GetXaxis()->SetTitle("E_{#nu} (GeV)");
     hframe->GetYaxis()->SetTitle("#sigma_{#nu} (1E-38 cm^{2}/GeV/nucleon)");
     hframe->Draw();
     //
     // Draw data and GENIE predictions and add a legend
     //
     for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {
       TGraphAsymmErrors * data = gXSecData[imode][idataset];
       if(!data) continue;
       data->Draw("P");
     }
     bestfit->Draw("l");
     nominal->Draw("l");
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
  }
  ps->Close();

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
