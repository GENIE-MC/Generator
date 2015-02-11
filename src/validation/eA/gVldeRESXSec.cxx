//____________________________________________________________________________
/*!

\program gvld_e_res_xsec

\brief   Compares GENIE with resonance electro-production cross-section data.

         The data come from the JLab Hall C archive maintained at:
         https://hallcweb.jlab.org/resdata/database/
         A local copy may be found in:
         $GENIE/data/validation/eA/xsec/differential/res/
         The archive contains ~12k data points.

         Syntax:
           gvld_e_res_xsec 
               [-d                             data_archive]
               [-s                             datasets_to_plot]
                --resonance-xsec-model         genie_model_name
                --non-resonance-bkg-xsec-model genie_model_name
               [--fit-params                   list_of_fit_params]
               [--fit-W2min                    W2min]
               [--fit-W2max                    W2max]

         Options:

           [] Denotes an optional argument.

           -d Full path to the electron cross-section data archive.
              By default, will pick the one at:
              $GENIE/data/validation/eA/xsec/differential/res/eRES.root

           -s Info on which datasets to plot. 
              By default, will pick the one at:
              $GENIE/data/validation/eA/xsec/differential/res/datasets.txt

           --resonance-xsec-model 
             Specify GENIE resonance electro-production cross-section model.

           --non-resonance-bkg-xsec-model 
             Specify GENIE electron DIS cross-section model.

           --fit-params
             Specify resonance cross-section model fit params.
             This is an optional argument. Don't specify a fit parameter list 
             if you only want to see a comparison of data with the nominal model.
             If a fit is performed, both the nominal and best-fit predictions are shown.
             An arbitrary number of fit parameters can be specified.
             The parameters must be present in the cross-section algorithm configuration.
             For each parameter, specify the name, the nominal/initial value, the minimum
             value, the maximum value and the step size in a comma separated list.
             Repeat for other fit parameters using a `:' to separate between them.
             For example:
		--fit-params norm,1.0,0.8,1.2,0.01:Mv,0.840,0.6,1.0,0.01

           --fit-W2min and 
           --fit-W2max
             Specify the W^2 range included in the fit.
             The default W^2 range is 1 GeV^2 to 2 GeV^2 (range including
             the Delta(1232) resonance).
 
         Example:

            % gvld_e_res_xsec 
                  --resonance-xsec-model genie::ReinSehgalRESPXSec/Default
                  -c genie::QPMDISPXSec/Default
         This is the default, so the plain command
            % gvld_e_res_xsec
         will do same thing

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created Oct 16, 2009 

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TLegend.h>
#include <TBox.h>
#include <TVirtualFitter.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Registry/Registry.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/KineUtils.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"

using std::ostringstream;
using std::ifstream;
using std::string;
using std::vector;

using namespace genie;
using namespace genie::constants;

//.................................................................................
// Utility class to hold info on plotted datasets
//
class eResDataSetDescription
{
public:
  eResDataSetDescription(
     int tgtpdg, string expt, double E, double theta) :
    fTgtPdg (tgtpdg), 
    fExpt   (expt), 
    fE      (E), 
    fTheta  (theta)
  {   
  }
  eResDataSetDescription() 
  {   
  }
  int    TgtPdg   (void) const { return fTgtPdg; }
  int    TgtZ     (void) const { return pdg::IonPdgCodeToZ(fTgtPdg); }
  int    TgtA     (void) const { return pdg::IonPdgCodeToA(fTgtPdg); }
  string TgtName  (void) const { 
    // all data are either in Hydrogen or Deuterium
    if(fTgtPdg == 1000010010) return "Hydrogen";
    if(fTgtPdg == 1000010020) return "Deuterium";
    return "Other";
  }
  string Expt     (void) const { return fExpt;  }
  double E        (void) const { return fE;     }
  double Theta    (void) const { return fTheta; }
  string LabelTeX (void) const { 
    ostringstream label;
    label << fExpt << " (" << this->TgtName() << "), ";
    label << "E = " << fE << " GeV, ";
    label << "#theta = " << fTheta << "^{o}";
    return label.str();
  }
private:
  int    fTgtPdg;  //
  string fExpt;    //
  double fE;       //
  double fTheta;   //
};
//.................................................................................

//
// constants
//

// defaults
const char * kDefDataArchiveFilename = "data/validation/eA/xsec/differential/res/eRES.root";  
const char * kDefDataSetsFilename    = "data/validation/eA/xsec/differential/res/datasets.txt";  

// plot config
const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

// number of resonances and GENIE resonance IDs
const int kNRes=18;  
Resonance_t kResId[kNRes] = {
   kP33_1232, kS11_1535, kD13_1520, kS11_1650,
   kD13_1700, kD15_1675, kS31_1620, kD33_1700,
   kP11_1440, kP33_1600, kP13_1720, kF15_1680,
   kP31_1910, kP33_1920, kF35_1905, kF37_1950,
   kP11_1710, kF17_1970 
};

// current program draws predictions only for the explicit resonance-production
// model at W<Wcut
const bool kDrawHatchcedScalingRegion = true; 

const double kWcut = 1.7; // Wcut from UserPhysicsOptions.xml

//
// globals
//

string gOptDataArchiveFilename = "";  // -d argument
string gOptDataSetsFilename    = "";  // -s argument
string gOptRESModelName        = "genie::ReinSehgalRESPXSec/Default";  // --resonance-xsec-model argument
string gOptDISModelName        = "genie::QPMDISPXSec/Default";  // --non-resonance-bkg-xsec-model argument
string gOptFitParams           = "";  // --fit-params argument
double gOptW2min               = 1.0; // --fit-W2min argument
double gOptW2max               = 2.0; // --fit-W2max argument

TFile *        gResDataFile  = 0;
TTree *        gResDataTree  = 0;
TPostScript *  gPS           = 0;
TCanvas *      gC            = 0;

XSecAlgorithmI * gRESXSecModel = 0; // resonance cross-section model
XSecAlgorithmI * gDISXSecModel = 0; // non-resonance bkg cross-section model

vector<eResDataSetDescription *> gDataSets; // list of plotted datasets

bool gOptFitEnabled = false;

Registry * fNominalResonanceParams = 0; // nominal configuration of resonance xsec model
Registry * fBestFitResonanceParams = 0; // best-fit configuration of resonance xsec model

vector<string> gFitParamName; // specified names of fit params
vector<double> gFitParamNom;  // specified nominal/initial values of fit params
vector<double> gFitParamMin;  // specified min values of fit params
vector<double> gFitParamMax;  // specified max values of fit params
vector<double> gFitParamStep; // specified step size of fit params

//
// function prototypes
//

void             Init               (void);
void             End                (void);
TGraphErrors *   Data               (unsigned int iset);
vector<TGraph *> Model              (unsigned int iset, unsigned int imodel);
void             Draw               (unsigned int iset);
void             DoTheFit           (void);
void             FitFunc            (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);
double           Chisquare          (double *par);
void             GetCommandLineArgs (int argc, char ** argv);
void             PrintSyntax        (void);

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);

  Init();

  // Fit resonance data, if fit params were specified in the command line
  if(gOptFitEnabled) {
     DoTheFit();
  }

  // Loop over data sets and plot data and corresponding GENIE predictions.
  // If a fit was enabled, plot both the nominal and best-fit predictions.
  for(unsigned int iset = 0; iset < gDataSets.size(); iset++) 
  {
    LOG("gvldtest", pNOTICE) 
      << "Producing plots for: " << gDataSets[iset]->LabelTeX();
    Draw(iset);
  }

  End();

  LOG("gvldtest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("gvldtest", pNOTICE) << "Initializing...";

  // Set GENIE style
  utils::style::SetDefaultStyle();

  //
  // Get TTree with electron-scattering data
  //
  if( ! utils::system::FileExists(gOptDataArchiveFilename) ) {
      LOG("gvldtest", pFATAL) 
         << "Can not find file: " << gOptDataArchiveFilename;
      gAbortingInErr = true;
      exit(1);
  }
  gResDataFile = new TFile(gOptDataArchiveFilename.c_str(),"read");  
  gResDataTree = (TTree *) gResDataFile->Get("resnt");
  if(!gResDataTree) {
      LOG("gvldtest", pFATAL) 
         << "Can not find TTree `resnt' in file: " << gOptDataArchiveFilename;
      gAbortingInErr = true;
      exit(1);
  }

  //
  // Read information on which data-sets to plot
  //
  LOG("gvldtest", pDEBUG) 
    << "Reading dataset summary info from: " << gOptDataSetsFilename;
  ifstream summary_file(gOptDataSetsFilename.c_str());
  if (!summary_file.good() ) {
      LOG("gvldtest", pFATAL) 
         << "Can't open data summary file: " << gOptDataSetsFilename;
      gAbortingInErr = true;
      exit(1);
  }
  while(1) {
      // skip header lines staring with #
      if(summary_file.peek() == '#') {
         summary_file.ignore(1000, '\n');
       } else {
         int    target = 0;
         string expt   = "";
	 double E      = 0;
         double theta  = 0;
         summary_file >> target >> expt >> E >> theta;
         summary_file.ignore(1000, '\n');
         if(summary_file.eof()) break;            

         LOG("gvldtest", pDEBUG) 
            << "target: " << target << ", experiment: " << expt 
            << ", E = " << E << " GeV, theta = " << theta << " deg";

         eResDataSetDescription * dataset = 
                 new eResDataSetDescription(target,expt,E,theta);
         gDataSets.push_back(dataset);
       }
  }
  summary_file.close();
  LOG("gvldtest", pNOTICE)
     << "Read "  << gDataSets.size() << " datasets";

  //
  // Get specified cross-section models
  //
  AlgFactory * algf = AlgFactory::Instance();
  gRESXSecModel = 0;
  if(gOptRESModelName != "none"){
     vector<string> modelv = utils::str::Split(gOptRESModelName,"/");
     assert(modelv.size()==2);
     string model_name = modelv[0];
     string model_conf = modelv[1];
     gRESXSecModel =
        dynamic_cast<XSecAlgorithmI *> (
            algf->AdoptAlgorithm(model_name, model_conf));
  }
  gDISXSecModel = 0;
  if(gOptDISModelName != "none"){
     vector<string> modelv = utils::str::Split(gOptDISModelName,"/");
     assert(modelv.size()==2);
     string model_name = modelv[0];
     string model_conf = modelv[1];
     gDISXSecModel =
        dynamic_cast<XSecAlgorithmI *> (
            algf->AdoptAlgorithm(model_name, model_conf));
  }

  //
  if(gOptFitEnabled) {
     assert(gRESXSecModel);
     gRESXSecModel->AdoptSubstructure(); 
     // print configuration of resonce model
     const Registry & r = gRESXSecModel->GetConfig();
     LOG("gvldtest", pNOTICE) << "Resonance model configuration: \n" << r;
     // find which parameters to fit
     vector<string> fitparams = utils::str::Split(gOptFitParams,":");
     assert(fitparams.size() >= 1);
     vector<string>::iterator iter = fitparams.begin();
     for( ; iter != fitparams.end(); ++iter) {
        string fp = *iter;
        vector<string> fptokens = utils::str::Split(fp,",");
        assert(fptokens.size() == 5);
        gFitParamName.push_back(fptokens[0]);
        gFitParamNom. push_back(atof(fptokens[1].c_str()));
        gFitParamMin. push_back(atof(fptokens[2].c_str()));
        gFitParamMax. push_back(atof(fptokens[3].c_str()));
        gFitParamStep.push_back(atof(fptokens[4].c_str()));
        LOG("gvldtest", pNOTICE) << "Fit parameter " << fptokens[0] << " was added";
     }
  }

  // Create plot canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // Get local time to tag outputs
  string lt_for_filename   = utils::system::LocalTimeAsString("%02d.%02d.%02d_%02d.%02d.%02d");
  string lt_for_cover_page = utils::system::LocalTimeAsString("%02d/%02d/%02d %02d:%02d:%02d");

  // Create output postscript file
  string filename  = Form("genie-e_res_data_comp-%s.ps",lt_for_filename.c_str());
  gPS = new TPostScript(filename.c_str(), 111);

  // Add cover page
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText("GENIE comparison with (e,e') data on 1H and 2D in the resonance region");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(lt_for_cover_page.c_str());
  hdr.AddText(" ");
  hdr.Draw();
  gC->Update();
}
//_________________________________________________________________________________
void End(void)
{
  LOG("gvldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;

  gResDataFile->Close();
}
//_________________________________________________________________________________
vector<TGraph *> Model(unsigned int iset, unsigned int imodel)
{
// Get GENIE predictions for the `iset' dataset 

  vector<TGraph *> model;

  LOG("gvldtest", pNOTICE) 
    << "Getting GENIE prediction (model ID = " 
    << imodel << ", data set ID = " << iset << ")";

  bool calc_res = (gRESXSecModel != 0);
  bool calc_dis = (gDISXSecModel != 0);
  bool calc_tot = calc_res && calc_dis;

  double M  = kNucleonMass;
  double M2 = M*M;

  double E     = gDataSets[iset]->E();
  double theta = gDataSets[iset]->Theta();
  double costh = TMath::Cos(2*kPi*theta/360.);

  LOG("gvldtest", pINFO) 
     << " E = " << E 
     << ", theta = " << theta << " (cos(theta) = " << costh << ")";

  int Z = gDataSets[iset]->TgtZ();
  int A = gDataSets[iset]->TgtA();
  int N = A-Z;
  bool   tgt_has_p = (Z>0);
  bool   tgt_has_n = (N>0);
  double frac_p = (double) Z / (double) A;
  double frac_n = (double) N / (double) A;

  const int n = 500;

  double W2_array[n];
  double d2sigRES_dEpdOmega_array[n]; 
  double d2sigDIS_dEpdOmega_array[n];
  double d2sigTOT_dEpdOmega_array[n];

  double Epmin = 0.01;
  double Epmax = E;
  double dEp = (Epmax-Epmin)/(n-1);

  for(int i=0; i<n; i++) {
     double Ep = Epmin + i*dEp;
     double Q2 = 2*E*Ep*(1-costh);
     double W2 = M2 + 2*M*(E-Ep)-Q2;
     double W  = TMath::Sqrt( TMath::Max(0.,W2) );

     double x  = 0;
     double y  = 0;
     utils::kinematics::WQ2toXY(Ep,M,W,Q2,x,y);          

     d2sigRES_dEpdOmega_array[i] = 0; 
     d2sigDIS_dEpdOmega_array[i] = 0;
     d2sigTOT_dEpdOmega_array[i] = 0;

     if(W2 <= 0) {
        LOG("gvldtest", pDEBUG) 
          << "Ep = " << Ep << ", Q2 = " << Q2 << ", W2 = " << W2
          << "... Skipping point";
        W2_array[i] = 0.;
        continue;
     }
     W2_array[i] = W2;

     // Will calculate  d^2 sigma / dW dQ^2 and then convert to d^2sigma / dE' dOmega
//     double jacobian = (E*Ep)*(M+2*E*(1-costh))/(kPi*W);
     double jacobian = (E*Ep*M)/(kPi*W);

     //
     // Calculate resonance cross-section
     //

     if(calc_res) {

        double d2sigRES_dWdQ2   = 0.;
        double d2sigRES_dWdQ2_p = 0.;
        double d2sigRES_dWdQ2_n = 0.;

        //
        // e+p -> e+Resonance
        //
        if(tgt_has_p) {
           Interaction * ep_res = 
                Interaction::RESEM(1000010010, kPdgProton,  kPdgElectron, E);
           ep_res->KinePtr()->SetW (W);
           ep_res->KinePtr()->SetQ2(Q2);
           // loop over resonances
           for(int ires=0; ires<kNRes; ires++) {
              ep_res->ExclTagPtr()->SetResonance(kResId[ires]);
              double xsec = gRESXSecModel->XSec(ep_res,kPSWQ2fE) / units::nb;
              xsec = TMath::Max(0., xsec);
              d2sigRES_dWdQ2_p += xsec;
              LOG("gvldtest", pINFO) 
                  << "d2xsec_dWdQ2(ep; " << utils::res::AsString(kResId[ires])
                  << "; E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
                  << "; Ep = " << Ep << " GeV, theta = " << theta << " deg) = " 
                  << xsec << " nbarn/GeV^3";
           }//res
           delete ep_res;
        }//has_p

        //
        // e+n -> e+Resonance
        //
        if(tgt_has_n) {
           Interaction * en_res = 
                Interaction::RESEM(1000000010, kPdgNeutron, kPdgElectron, E);
           en_res->KinePtr()->SetW (W);
           en_res->KinePtr()->SetQ2(Q2);
           // loop over resonances
           for(int ires=0; ires<kNRes; ires++) {
              en_res->ExclTagPtr()->SetResonance(kResId[ires]);
              double xsec = gRESXSecModel->XSec(en_res,kPSWQ2fE) / units::nb;
              xsec = TMath::Max(0., xsec);
              d2sigRES_dWdQ2_n += xsec;
              LOG("gvldtest", pINFO) 
                  << "d2xsec_dWdQ2(en; " << utils::res::AsString(kResId[ires])
                  << "; E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
                  << "; Ep = " << Ep << " GeV, theta = " << theta << " deg) = " 
                  << xsec << " nbarn/GeV^3";
           }//res
           delete en_res;
        }//has_n

        // sum over resonances and averaged over p,n
        d2sigRES_dWdQ2 = (frac_p * d2sigRES_dWdQ2_p + frac_n * d2sigRES_dWdQ2_n);

        // convert to d^2sigma / dE' dOmega and store
        double d2sigRES_dEpdOmega = jacobian * d2sigRES_dWdQ2;
        if(TMath::IsNaN(d2sigRES_dEpdOmega)) {
           LOG("gvldtest", pWARN) << "Got a NaN!";
           d2sigRES_dEpdOmega = 0;
        }
        d2sigRES_dEpdOmega_array[i] = TMath::Max(0., d2sigRES_dEpdOmega);

     }//calc_res

     //
     // Calculate DIS cross-section
     //

     if(calc_dis) {
        double d2sigDIS_dWdQ2   = 0.;
        double d2sigDIS_dWdQ2_p = 0.;
        double d2sigDIS_dWdQ2_n = 0.;
        if(tgt_has_p) {
           // Note: Not setting quark ID
           // If the quark ID is set, the code returns (eg for neutrino CC) the vq->lq' cross-section.
           // But if a quark ID is not specified then the code loops over all relevant
           // valence and sea quark species and returns (eg for neutrino CC) the vN->lX cross-section.
           Interaction * ep_dis = 
               Interaction::DISEM(1000010010, kPdgProton,  kPdgElectron, E); 
           ep_dis->KinePtr()->SetW (W);
           ep_dis->KinePtr()->SetQ2(Q2);
           ep_dis->KinePtr()->Setx (x);
           ep_dis->KinePtr()->Sety (y);
           // Note:
           d2sigDIS_dWdQ2_p = gDISXSecModel->XSec(ep_dis,kPSWQ2fE) / units::nb;
           d2sigDIS_dWdQ2_p = TMath::Max(0., d2sigDIS_dWdQ2_p);
           delete ep_dis;
        }
        if(tgt_has_n) {
           Interaction * en_dis = 
               Interaction::DISEM(1000000010, kPdgNeutron, kPdgElectron, E);
           en_dis->KinePtr()->SetW (W);
           en_dis->KinePtr()->SetQ2(Q2);
           en_dis->KinePtr()->Setx (x);
           en_dis->KinePtr()->Sety (y);
           d2sigDIS_dWdQ2_n = gDISXSecModel->XSec(en_dis,kPSWQ2fE) / units::nb;
           d2sigDIS_dWdQ2_n = TMath::Max(0., d2sigDIS_dWdQ2_n);
           delete en_dis;
        }
        d2sigDIS_dWdQ2 = (frac_p * d2sigDIS_dWdQ2_p + frac_n * d2sigDIS_dWdQ2_n);

        // convert to d^2sigma / dE' dOmega and store
        double d2sigDIS_dEpdOmega = jacobian * d2sigDIS_dWdQ2;
        if(TMath::IsNaN(d2sigDIS_dEpdOmega)) {
           LOG("gvldtest", pWARN) << "Got a NaN!";
           d2sigDIS_dEpdOmega = 0;
        }
        d2sigDIS_dEpdOmega_array[i] = TMath::Max(0., d2sigDIS_dEpdOmega);

     }//calc_dis

     //
     // Calculate total cross-section
     //
     if(calc_tot) {
        d2sigTOT_dEpdOmega_array[i] = 
           d2sigRES_dEpdOmega_array[i] +
           d2sigDIS_dEpdOmega_array[i]; 
     }
   
  }//i

  // Create graphs & store them in array
  if(calc_tot) {
     TGraph * gr = new TGraph(n,W2_array,d2sigTOT_dEpdOmega_array);
     const char * title = 0;
     int lstyle = kSolid;
     if(gOptFitEnabled) {
       if      (imodel==0) { title = "Total (Best-fit)"; }
       else if (imodel==1) { title = "Total (Nominal)"; lstyle = kDashed; }
     } else {
       title = "Total";
     }
     gr->SetLineColor(kBlack);
     gr->SetLineStyle(lstyle);
     gr->SetTitle(title);
     model.push_back(gr);
  }
  if(calc_res) {
     TGraph * gr = new TGraph(n,W2_array,d2sigRES_dEpdOmega_array);
     const char * title = 0;
     int lstyle = kSolid;
     if(gOptFitEnabled) {
       if      (imodel==0) { title = "Resonance (Best-fit)"; }
       else if (imodel==1) { title = "Resonance (Nominal)"; lstyle = kDashed; }
     } else {
       title = "Resonance";
     }
     gr->SetLineColor(kRed);
     gr->SetLineStyle(lstyle);
     gr->SetTitle(title);
     model.push_back(gr);
  }
  if(calc_dis) {
     TGraph * gr = new TGraph(n,W2_array,d2sigDIS_dEpdOmega_array);
     const char * title = 0;
     int lstyle = kSolid;
     if(gOptFitEnabled) {
       if      (imodel==0) { title = "Non-resonance bkg / DIS (Best-fit)"; }
       else if (imodel==1) { title = "Non-resonance bkg / DIS (Nominal)"; lstyle = kDashed; }
     } else {
       title = "Non-resonance bkg / DIS";
     }
     gr->SetLineColor(kBlue);
     gr->SetLineStyle(lstyle);
     gr->SetTitle(title);
     model.push_back(gr);
  }
  
  return model;
}
//_________________________________________________________________________________
TGraphErrors * Data(unsigned int iset)
{
  const double dE      = 1.0E-3;
  const double dtheta  = 2.5E-2;

  double E     = gDataSets[iset]->E();
  double theta = gDataSets[iset]->Theta();

  int Z = gDataSets[iset]->TgtZ();
  int A = gDataSets[iset]->TgtA();

  const char * selection = 
    Form("E > %f && E < %f && theta > %f && theta < %f && Z == %d && A == %d",
         E     - dE,
         E     + dE,
         theta - dtheta,
         theta + dtheta,
         Z,A);

  gResDataTree->Draw("W2:xsec:xsec_err", selection, "goff");
  
  int n = gResDataTree->GetSelectedRows();

  LOG("gvldtest", pNOTICE) 
    << "Found " << n << " data points in the xsec archive";

  if(n == 0) return 0; // return null graph

  // Data returned by TTree::Draw() are not necessarily ordered in W
  // Do the ordering here before building the graph
  int    *  idx = new int   [n];
  double *  xv  = new double[n];
  double *  yv  = new double[n];
  double *  dyv = new double[n];

  TMath::Sort(n,gResDataTree->GetV1(),idx,false);

  for(int i=0; i<n; i++) {
     int ii = idx[i];
     xv [i] = (gResDataTree->GetV1())[ii];
     yv [i] = (gResDataTree->GetV2())[ii];
     dyv[i] = (gResDataTree->GetV3())[ii];
  }

  TGraphErrors * gr = new TGraphErrors(n,xv,yv,0,dyv);
  utils::style::Format(gr, 1,1,1,1,8,0.4);

  delete [] idx;
  delete [] xv;
  delete [] yv;
  delete [] dyv;

  return gr;
}
//_________________________________________________________________________________
void Draw(unsigned int iset)
{
  // get all measurements for the current channel from the NuValidator MySQL dbase
  TGraphErrors * data = Data(iset);
  if(!data) return;

  // get the corresponding GENIE model prediction
  vector< vector<TGraph *> > model(2);
  if(gOptFitEnabled) {
    gRESXSecModel->Configure(*fBestFitResonanceParams); 
    model[0] = Model(iset,0);
    gRESXSecModel->Configure(*fNominalResonanceParams); 
    model[1] = Model(iset,1);
  } else {
    model[0] = Model(iset,0);
  }

  int plots_per_page = kNCx * kNCy;
  int iplot = 1 + iset % plots_per_page;

  if(iplot == 1) {
     gPS->NewPage();
     gC -> Clear();
     gC -> Divide(kNCx,kNCy);
  }

  gC -> GetPad(iplot) -> Range(0,0,100,100);
  gC -> GetPad(iplot) -> SetFillColor(0);
  gC -> GetPad(iplot) -> SetBorderMode(0);
  gC -> GetPad(iplot) -> cd();

  TH1F * hframe = 0;
  double xmin =  9999999999, scale_xmin = 0.50;
  double xmax = -9999999999, scale_xmax = 1.05;
  double ymin =  9999999999, scale_ymin = 0.40;
  double ymax = -9999999999, scale_ymax = 1.30;
  // get W2 and d2sigma/dEdOmega range in the the data
  xmin  = ( data->GetX() )[TMath::LocMin(data->GetN(),data->GetX())];
  xmax  = ( data->GetX() )[TMath::LocMax(data->GetN(),data->GetX())];
  ymin  = ( data->GetY() )[TMath::LocMin(data->GetN(),data->GetY())];
  ymax  = ( data->GetY() )[TMath::LocMax(data->GetN(),data->GetY())];
  xmin  = TMath::Max(xmin, 0.5); // some data go very low 
  // take also into account the d2sigma/dEdOmega range in the the models
  // (for the W2 range of the dataset) just in case data and MC are not that similar...
  for(unsigned int imodel=0; imodel < model.size(); imodel++) {
    for(unsigned int imode=0; imode < model[imodel].size(); imode++) {
      TGraph * mm = model[imodel][imode];
      if(mm) {
        for(int k=0; k<mm->GetN(); k++) {
          double x = (mm->GetX())[k];
          if(x < xmin || x > xmax) continue;
          ymin = TMath::Min(ymin, (mm->GetY())[k]);
          ymax = TMath::Max(ymax, (mm->GetY())[k]);
        }//k
      }//mm
    }//imode
  }//imodel
  LOG("gvldtest", pDEBUG) 
    << "Plot range:" 
    << "W^{2} = [" << xmin << ", " << xmax << "] GeV^{2}, "
    << "d^{2}#sigma / d#Omega dE = [" << ymin << ", " << ymax << "]  nb/sr/GeV";

  hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(
        scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
  hframe->Draw();
  hframe->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
  hframe->GetYaxis()->SetTitle("d^{2}#sigma / d#Omega dE (nb/sr/GeV)");

  // draw data and GENIE models
  data->Draw("P");
  for(unsigned int imodel=0; imodel < model.size(); imodel++) {
    for(unsigned int imode=0; imode < model[imodel].size(); imode++) {
       TGraph * mm = model[imodel][imode];
       if(!mm) continue;
       mm->Draw("L");
    }
  }

  // add legend
  double lymin = (gOptFitEnabled) ? 0.65 : 0.75;
  TLegend * legend = new TLegend(0.20, lymin, 0.50, 0.85);
  legend->SetLineStyle(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.025);
  legend->SetHeader("GENIE");
  for(unsigned int imodel=0; imodel < model.size(); imodel++) {
    for(unsigned int imode=0; imode < model[imodel].size(); imode++) {
       TGraph * mm = model[imodel][imode];
       if(!mm) continue;
       legend->AddEntry(mm,mm->GetTitle(),"L");
    }
  }
  legend->Draw();

  // scaling region
  TBox * scaling_region = 0;
  if(kDrawHatchcedScalingRegion) {
    double W2c = kWcut*kWcut;
    if(W2c > scale_xmin*xmin && W2c < scale_xmax*xmax) {
       scaling_region = new TBox(
           W2c, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
//       scaling_region->SetFillColor(kRed);
//       scaling_region->SetFillStyle(3905);
       scaling_region->SetFillStyle(0);
       scaling_region->SetLineColor(kRed);
       scaling_region->SetLineStyle(kDashed);
       scaling_region->Draw();
    }
  }

  // some data show the elastic peak - mark the are to avoid confusion
  if(xmin < 1) {
    double Wm2 = 1.21; // between the QE and Delta peaks
    TBox * qe_peak = new TBox(
       scale_xmin*xmin, scale_ymin*ymin, Wm2, scale_ymax*ymax);
     qe_peak->SetFillColor(kBlue);
     qe_peak->SetFillStyle(3005);
     qe_peak->Draw();
  }

  // title
  TLatex * title = new TLatex(
     scale_xmin*xmin + 0.2*(scale_xmax*xmax-scale_xmin*xmin),
    1.01*scale_ymax*ymax,gDataSets[iset]->LabelTeX().c_str());
  title->SetTextSize(0.027);
  title->Draw();

  gC->GetPad(iplot)->Update();
  gC->Update();
}
//_________________________________________________________________________________
void DoTheFit(void)
{
  if(!gOptFitEnabled) return;

  // Store nominal params
  fNominalResonanceParams = new Registry(gRESXSecModel->GetConfig());

  // Find number of fit params and get a fitter
  const unsigned int nfitparams = gFitParamName.size();

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * fitter = TVirtualFitter::Fitter(0,nfitparams);
      
  double arglist[100];
  arglist[0] = -1;
  fitter->ExecuteCommand("SET PRINT",arglist,1);

  // Set parameters
  for(unsigned int i = 0; i < nfitparams; i++) 
  {
    string name    = gFitParamName[i];
    double nominal = gFitParamNom [i];
    double min     = gFitParamMin [i];
    double max     = gFitParamMax [i];
    double step    = gFitParamStep[i];

    LOG("gtune", pNOTICE)
        << "** Setting fit param " << i
          << " (" << name << ") initial value = " << nominal
             << ", range = [" << min << ", " << max <<"]";

    fitter->SetParameter(i, name.c_str(), nominal, step, min, max);
  }
  fitter->SetFCN(FitFunc);

  // MINUIT minimization step

  LOG("gtune", pNOTICE)
        << "** Running fitter considering data-points in the W^2 range = [" 
        << gOptW2min << ", " << gOptW2min << "] GeV^2";

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

  // Get and store best-fit values
  fBestFitResonanceParams = new Registry(gRESXSecModel->GetConfig());

  // Print results
  double amin = 0;
  fitter->PrintResults(3,amin);

}
//_________________________________________________________________________________
void FitFunc (
        Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
// MINUIT fit function with signature expected by TVirtualFitter::SetFCN()

  f = Chisquare(par);
}
//_________________________________________________________________________________
double Chisquare(double *par)
{
  //
  // Reconfigure resonance cross-section algorithm
  //
  Registry r(gRESXSecModel->GetConfig());
  r.UnLock();
  const unsigned int nfitparams = gFitParamName.size();
  for(unsigned int i = 0; i < nfitparams; i++) {
    string pname  = gFitParamName[i];
    double pvalue = par[i];
    r.Set(pname,pvalue);
  }
  gRESXSecModel->Configure(r);
  LOG("gvldtest", pNOTICE) 
      << "Computing \\chi^{2} using RES xsec config: " << r;

  //
  // Compute chisq using updated model
  //
  double chisq = 0;

  // loop over data sets 
  for(unsigned int iset = 0; iset < gDataSets.size(); iset++) 
  {
    LOG("gvldtest", pINFO) 
      << "Including dataset: " << gDataSets[iset]->LabelTeX();

    // get the data
    TGraphErrors * data = Data(iset);
    if(!data) continue;

    // get the corresponding GENIE model prediction
    vector<TGraph *> model = Model(iset,0);

    // loop over data points and calculate contribution to chisq
    for(int ipoint = 0; ipoint < data->GetN(); ipoint++) {
        double xd = data->GetX() [ipoint]; // W^2
        // consider only the data points in the specified W^2 range
	if(xd < gOptW2min) continue;        
	if(xd > gOptW2max) continue;
        double yd     = data->GetY() [ipoint];
        double yd_err = data->GetEY()[ipoint];
        double ym     = model[0]->Eval(xd);
        chisq += TMath::Power((yd-ym)/yd_err, 2.);
    }
    delete data;
    for(unsigned int i=0; i<model.size(); i++) {
	delete model[i];
        model[i]=0;
    }
    model.clear();

  }//iset

  LOG("gvldtest", pNOTICE) << "Chisq = " << chisq;

  return chisq;
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  if(parser.OptionExists('d')){
     string filename = parser.ArgAsString('d');
     gOptDataArchiveFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataArchiveFilename;
        gOptDataArchiveFilename = filename;
     } else { 
        LOG("gvldtest", pFATAL) 
          << "\n Please make sure that $GENIE is defined, or use the -d option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  if(parser.OptionExists('s')){
     string filename = parser.ArgAsString('s');
     gOptDataSetsFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataSetsFilename;
        gOptDataSetsFilename = filename;
     } else { 
        LOG("gvldtest", pFATAL) 
          << "\n Please make sure that $GENIE is defined, or use the -s option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  // Get GENIE model names to be used
  if(parser.OptionExists("resonance-xsec-model")){
     gOptRESModelName = parser.Arg("resonance-xsec-model");
  }   
  if(parser.OptionExists("non-resonance-bkg-xsec-model")){
     gOptDISModelName = parser.Arg("non-resonance-bkg-xsec-model");
  }

  if(parser.OptionExists("fit-params")) {
     gOptFitEnabled = true;
     gOptFitParams = parser.Arg("fit-params");
  }

  if(parser.OptionExists("fit-W2min")) {
     string sW2 = parser.Arg("fit-W2min");
     gOptW2min = atof(sW2.c_str());
  }

  if(parser.OptionExists("fit-W2max")) {
     string sW2 = parser.Arg("fit-W2max");
     gOptW2max = atof(sW2.c_str());
  }

}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" 
    << "\n"
    << " gvld_e_res_xsec \n"
    << "       --resonance-xsec-model         genie_model \n"
    << "       --non-resonance-bkg-xsec-model genie_model \n"
    << "      [--fit-params parameters_to_fit_for_if_any] \n"
    << "      [-d data_archive] \n"
    << "      [-s datasets_to_plot] \n";
}
//_________________________________________________________________________________
