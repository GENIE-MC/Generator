//____________________________________________________________________________
/*!

\program gvld_sf

\brief   Compares GENIE F2 and xF3 structure-function models with data.

         Syntax:
           gvld_sf 
                [-d data_archive] 
                [-s list_of_data_mc_comparisons]
                [--resonance-xsec-model genie::model_name/model_config]
                [--dis-xsec-model       genie::model_name/model_config]
                [--dis-charm-xsec-model genie::model_name/model_config]

         Options:
           [] Denotes an optional argument.

           -d Location of the neutrino cross-section data archive.
              By default, the program will look-up the one located in:
              $GENIE/data/validation/vA/sf/structFunc.root

           -s List of data/MC comparisons to perform.
              By default, the program will look-up the one located in:
              $GENIE/data/validation/vA/sf/comparisons_xbj.txt

           --resonance-xsec-model
              Specify GENIE resonance cross-section model.

           --dis-xsec-model
              Specify GENIE DIS cross-section model.

           --dis-charm-xsec-model
              Specify GENIE DIS charm production cross-section model.

         Example:

            % gvld_sf \
                --resonance-xsec-model genie::ReinSehgalRESPXSec/Default \
                --dis-xsec-model genie::QPMDISPXSec/Default \
		--dis-charm-xsec-model genie::AivazisCharmPXSecLO/CC-Default

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Jelena Ilic <jelena.ilic \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 06, 2008 

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

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"
#include "validation/StructFunc/StructFunc.h"

using std::ostringstream;
using std::ifstream;
using std::string;
using std::vector;

using namespace genie;
using namespace genie::constants;
using namespace genie::mc_vs_data;

//.................................................................................
// Utility class to hold info on data/MC comparisons
//
class SFCompInfo
{
public:
  SFCompInfo(
    string fi, string datasets, double x, double Q2min, double Q2max, 
    int lepton_pdg, int target_pdg, double scale) :
   fFi(fi),
   fDataSetKeys(datasets),
   fxBj(x),
   fQ2min(Q2min),
   fQ2max(Q2max),
   fLeptonPDG(lepton_pdg),
   fTargetPDG(target_pdg),
   fScale(scale)
   { 
   }
  ~SFCompInfo() 
   { 
   }

   string Fi          (void) const { return fFi;          }
   string DataSetKeys (void) const { return fDataSetKeys; }
   double xBjorken    (void) const { return fxBj;         }
   double Q2Min       (void) const { return fQ2min;       }
   double Q2Max       (void) const { return fQ2max;       }
   double Scale       (void) const { return fScale;       }
   int    LeptonPDG   (void) const { return fLeptonPDG;   }
   int    TargetPDG   (void) const { return fTargetPDG;   }

   string Label (void) const 
   { 
      return Form("x_{Bj} = %.4f", fxBj);
   }

private:
  string fFi;
  string fDataSetKeys;
  double fxBj;
  double fQ2min;
  double fQ2max;
  int    fLeptonPDG;
  int    fTargetPDG;
  double fScale;

};
//.................................................................................

//
// constants
//

// defaults
const char * kDefDataArchiveFilename = "data/validation/vA/sf/structFunc.root";  
const char * kDefCompInfoFilename    = "data/validation/vA/sf/comparisons_xbj.txt";  

const double kQ2min =   0.1; // GeV^2
const double kQ2max = 150.0; // GeV^2

// plot config
const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

const int kNMaxDataSets = 20; // max number of datasets in single plot

const int kDataPointStyle[kNMaxDataSets] = 
{ 
  20,    20,       20,    20,
  21,    21,       21,    21,
  24,    24,       24,    24,
  25,    25,       25,    25,
  29,    29,       29,    29
};
const int kDataPointColor[kNMaxDataSets] = 
{
  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kRed,  kGreen+1, kBlue, kMagenta+1, 
  kRed,  kGreen+1, kBlue, kMagenta+1
};

const int kNumOfSummaryPlotColors = 5;

const int SummaryPlotColor[kNumOfSummaryPlotColors] =
{
  kBlack, kRed, kGreen+1, kBlue, kMagenta+1
};

//
// globals
//

string gOptDataArchiveFilename = ""; // -d command-line argument
string gOptCompInfoFilename    = ""; // -s command-line argument
string gOptRESModelName        = ""; // --resonance-xsec-model command-line argument
string gOptDISModelName        = ""; // --dis-xsec-model command-line argument
string gOptDISCharmModelName   = ""; // --dis-charm-xsec-model command-line argument

TFile *        gSFDataFile  = 0;
TTree *        gSFDataTree  = 0;
TPostScript *  gPS          = 0;
TCanvas *      gC           = 0;

const XSecAlgorithmI * gRESXSecModel      = 0; // resonance cross-section model
const XSecAlgorithmI * gDISXSecModel      = 0; // DIS cross-section model
const XSecAlgorithmI * gDISCharmXSecModel = 0; // DIS cross-section model (charm production)
vector<SFCompInfo *>   gComparisons;      // info on which comparisons to perform
StructFunc             gStructFuncCalc;   // utility class extracting structure functions from the cross-section model

// function prototypes
void     Init               (void);
void     End                (void);
void     Draw               (void);
TH1F *   DrawFrame          (TGraph * gr0, TGraph * gr1, TCanvas * c);
void     GetCommandLineArgs (int argc, char ** argv);
void     PrintSyntax        (void);

TGraph *               Model(unsigned int icomp, unsigned int imodel, double scale = 1.0);
vector<TGraphErrors *> Data (unsigned int icomp, double scale = 1.0);

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);
  Init();
  Draw();
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
  // Get TTree with structure-function data
  //
  if( ! utils::system::FileExists(gOptDataArchiveFilename) ) {
     LOG("gvldtest", pFATAL) 
       << "Can not find file: " << gOptDataArchiveFilename;
     gAbortingInErr = true;
     exit(1);
  }
  gSFDataFile = new TFile(gOptDataArchiveFilename.c_str(),"read");  
  gSFDataTree = (TTree *) gSFDataFile->Get("sfnt");
  if(!gSFDataTree) {
      LOG("gvldtest", pFATAL) 
        << "Can not find TTree `sfnt' in file: " << gOptDataArchiveFilename;
      gAbortingInErr = true;
      exit(1);
  }

  //
  // Get cross-section models
  //
  AlgFactory * algf = AlgFactory::Instance();
  gRESXSecModel = 0;
  if(gOptRESModelName != "none"){
      vector<string> modelv = utils::str::Split(gOptRESModelName,"/");
      assert(modelv.size()==2);
      string model_name = modelv[0];
      string model_conf = modelv[1];
      gRESXSecModel =dynamic_cast<const XSecAlgorithmI *> (
         algf->GetAlgorithm(model_name, model_conf));
  }
  gDISXSecModel = 0;
  if(gOptDISModelName != "none"){
       vector<string> modelv = utils::str::Split(gOptDISModelName,"/");
       assert(modelv.size()==2);
       string model_name = modelv[0];
       string model_conf = modelv[1];
       gDISXSecModel = dynamic_cast<const XSecAlgorithmI *> (
          algf->GetAlgorithm(model_name, model_conf));
  }
  gDISCharmXSecModel = 0;
  if(gOptDISCharmModelName != "none"){
       vector<string> modelv = utils::str::Split(gOptDISCharmModelName,"/");
       assert(modelv.size()==2);
       string model_name = modelv[0];
       string model_conf = modelv[1];
       gDISCharmXSecModel = dynamic_cast<const XSecAlgorithmI *> (
          algf->GetAlgorithm(model_name, model_conf));
  }

  gStructFuncCalc.SetResonanceXSecModel (gRESXSecModel);
  gStructFuncCalc.SetDISXSecModel       (gDISXSecModel);
  gStructFuncCalc.SetDISCharmXSecModel  (gDISCharmXSecModel);

  //
  // Read info on comparisons to perform
  //
  LOG("gvldtest", pNOTICE) 
        << "Reading dataset summary info from: " << gOptCompInfoFilename;
  ifstream summary_file(gOptCompInfoFilename.c_str());
  if (!summary_file.good() ) {
       LOG("gvldtest", pFATAL) 
          << "Can't open data summary file: " << gOptCompInfoFilename;
       gAbortingInErr = true;
       exit(1);
  }
  while(1) {
      // skip header lines staring with #
      if(summary_file.peek() == '#') {
          summary_file.ignore(1000, '\n');
      } else {
          string Fi="", datasets="";
          double xbj = 0.;
          int lpdg=0, tpdg=0;
          summary_file >> Fi >> xbj >> lpdg >> tpdg >> datasets;
          summary_file.ignore(1000, '\n');
          if(summary_file.eof()) break;            
          SFCompInfo * comparison = 
              new SFCompInfo(Fi, datasets, xbj, kQ2min, kQ2max, lpdg, tpdg, 1.);
          gComparisons.push_back(comparison);
      }//new non# line
  }//end of summary file
  summary_file.close();

  LOG("gvldtest", pNOTICE) 
      << "Read "  << gComparisons.size() << " datasets";

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
  string filename  = Form("genie-structfunc_data_comp-%s.ps",lt_for_filename.c_str());
  gPS = new TPostScript(filename.c_str(), 111);

  // Add cover page
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE comparison with F2, xF3 structure function data");
  hdr.AddText(" ");
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

  gSFDataFile->Close();
}
//_________________________________________________________________________________
TGraph * Model(unsigned int icomp, unsigned int imodel, double scale)
{
  LOG("gvldtest", pNOTICE) 
      << "Getting GENIE prediction (comparison ID = " 
      << icomp << ", model ID = " << imodel << ")";

  int    lepton_pdg = gComparisons[icomp]->LeptonPDG();
  int    target_pdg = gComparisons[icomp]->TargetPDG();
  double xBj        = gComparisons[icomp]->xBjorken();
  double Q2min      = gComparisons[icomp]->Q2Min();
  double Q2max      = gComparisons[icomp]->Q2Max();
  string Fi         = gComparisons[icomp]->Fi();

  LOG("gvldtest", pNOTICE) 
      << "lepton = " << lepton_pdg << " + target = " << target_pdg << ", " 
      << Fi << " for x_{Bj} = " << xBj 
      << " and Q2 in " << Q2min << " to " << Q2max << " GeV^2 range";

  const int knQ2 = 250;
  const double logQ2min = TMath::Log10(Q2min);
  const double logQ2max = TMath::Log10(Q2max);
  const double dlogQ2   = (logQ2max-logQ2min)/(knQ2-1);

  double * Q2_array = new double[knQ2]; 
  double * Fi_array = new double[knQ2]; 

  for (int i = 0; i < knQ2; i++) 
  {
    double Q2  = TMath::Power(10., logQ2min + i * dlogQ2);
    double F1  = 0.;
    double F2  = 0.;
    double xF3 = 0.;
    gStructFuncCalc.ExtractF1F2xF3(xBj, Q2, lepton_pdg, target_pdg, F1, F2, xF3);
    F1  = TMath::Max(0., F1 );
    F2  = TMath::Max(0., F2 );
    xF3 = TMath::Max(0., xF3);
    LOG("gvldtest", pNOTICE) 
       << "x_Bj = " << xBj << ", Q2 = " << Q2 << " GeV^2 : "
       << "F1 = " << F1 << ", F2 = " << F2 << ", xF3 = " << xF3;
    Q2_array[i] = Q2;
    Fi_array[i] = 0.;
    if(Fi == "F2" ) Fi_array[i] = scale * F2;
    if(Fi == "xF3") Fi_array[i] = scale * xF3;
  }

  TGraph * grFi = new TGraph(knQ2,Q2_array,Fi_array);
  grFi->SetLineColor(kBlack);
  grFi->SetLineStyle(kSolid);
  grFi->SetLineWidth(1);

  delete [] Q2_array;
  delete [] Fi_array;

  return grFi;
}
//_________________________________________________________________________________
vector<TGraphErrors *> Data(unsigned int icomp, double scale)
{
  LOG("gvldtest", pNOTICE) 
      << "Retrieving expt data (comparison ID = " << icomp << ")";

  double epsilon = 1E-5;

//int    lpdg  = gComparisons[icomp]->LeptonPDG();
//int    tpdg  = gComparisons[icomp]->TargetPDG();
  double xBj   = gComparisons[icomp]->xBjorken();
//double Q2min = gComparisons[icomp]->Q2Min();
//double Q2max = gComparisons[icomp]->Q2Max();
  string Fi    = gComparisons[icomp]->Fi();
  string keys  = gComparisons[icomp]->DataSetKeys();

  vector<string> keyv = utils::str::Split(keys,";");
  unsigned int ndatasets = keyv.size();

  vector<TGraphErrors *> data(ndatasets);

  for(unsigned int idataset = 0; idataset < ndatasets; idataset++) {

      const char * selection = Form("(x>%f) && (x<%f) && (Ftype==\"%s\") && (dataset==\"%s\")",
           xBj-epsilon, xBj+epsilon, Fi.c_str(), keyv[idataset].c_str());

      gSFDataTree->Draw("Q2:F:dFp", selection, "goff");
      int np = gSFDataTree->GetSelectedRows();
      LOG("gvldtest", pNOTICE) 
          << "Found " << np << " points in the data archive";

      // Data returned by TTree::Draw() are not necessarily ordered in W
      // Do the ordering here before building the graph
      int    *  idx = new int   [np];
      double *  xv  = new double[np];
      double *  yv  = new double[np];
      double *  dyv = new double[np];
      TMath::Sort(np,gSFDataTree->GetV1(),idx,false);
      for(int i = 0; i < np; i++) {
         int ii = idx[i];
         xv [i] = (gSFDataTree->GetV1())[ii];
         yv [i] = (gSFDataTree->GetV2())[ii] * scale;
         dyv[i] = (gSFDataTree->GetV3())[ii] * scale;
      }

      TGraphErrors * gr = new TGraphErrors(np,xv,yv,0,dyv);
      int sty = kDataPointStyle[idataset];
      int col = kDataPointColor[idataset];
      utils::style::Format(gr,col,kSolid,1,col,sty,0.7);
      string title = keyv[idataset].substr(0, keyv[idataset].find_first_of(","));
      gr->SetTitle(title.c_str());
      data.push_back(gr);

      delete [] idx;
      delete [] xv;
      delete [] yv;
      delete [] dyv;
  }
  return data;
}
//_________________________________________________________________________________
void Draw(void)
{
  //
  // Draw summary plots first, with all F2 (xF3) datasets and GENIE predictions 
  // shown on the same plot.
  //

  // > F2
  {
    gPS-> NewPage();
    gC -> SetFillColor(0);
    gC -> SetBorderMode(0);
    gC -> SetLogx(1);
    gC -> SetLogy(1);

    TH1F * hframe = 0;
    double xmin =  kQ2min;
    double xmax =  kQ2max*12;
    double ymin =  1E-2;
    double ymax =  2E+3;
    hframe = (TH1F*) gC->DrawFrame(xmin, ymin, xmax, ymax);
    hframe->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    hframe->GetYaxis()->SetTitle("F2");

    // loop over data sets and plot data and corresponding GENIE predictions
    for(unsigned int icomp = 0; icomp < gComparisons.size(); icomp++) 
    {
      LOG("gvldtest", pNOTICE) 
         << "Producing plots for: " << gComparisons[icomp]->Label();
      int color = SummaryPlotColor[(icomp%kNumOfSummaryPlotColors)];
      double xBj = gComparisons[icomp]->xBjorken();
      string Fi  = gComparisons[icomp]->Fi();
      if(Fi != "F2") continue;
      double scale = 1.0;
      if(xBj > 0.000 && xBj <= 0.008) scale = 500.0;
      if(xBj > 0.008 && xBj <= 0.013) scale = 350.0;
      if(xBj > 0.013 && xBj <= 0.016) scale = 250.0;
      if(xBj > 0.016 && xBj <= 0.020) scale = 150.0;
      if(xBj > 0.020 && xBj <= 0.030) scale = 100.0;
      if(xBj > 0.030 && xBj <= 0.040) scale =  60.0;
      if(xBj > 0.040 && xBj <= 0.049) scale =  40.0;
      if(xBj > 0.049 && xBj <= 0.060) scale =  30.0;
      if(xBj > 0.060 && xBj <= 0.075) scale =  20.0;
      if(xBj > 0.075 && xBj <= 0.085) scale =  15.0;
      if(xBj > 0.085 && xBj <= 0.095) scale =  10.0;
      if(xBj > 0.095 && xBj <= 0.115) scale =   6.0;
      if(xBj > 0.115 && xBj <= 0.130) scale =   4.0;
      if(xBj > 0.130 && xBj <= 0.150) scale =   3.0;
      if(xBj > 0.150 && xBj <= 0.177) scale =   2.0;
      if(xBj > 0.177 && xBj <= 0.182) scale =   1.6;
      if(xBj > 0.182 && xBj <= 0.190) scale =   1.3;
      if(xBj > 0.190 && xBj <= 0.230) scale =   1.1;
      // Get expt data from GENIE archive
      vector<TGraphErrors *> data = Data(icomp,scale);
      if(data.size() == 0) continue;
      // Get the corresponding GENIE model prediction
      TGraph* model = Model(icomp,0,scale);
      for(unsigned int i = 0; i< data.size(); i++) {
         if(!data[i]) continue;
         data[i]->SetMarkerColor(color);
         data[i]->SetLineColor(color);
         data[i]->Draw("P");
       //legend->AddEntry(data[i], data[i]->GetTitle(), "P");
      }  
      if(model) {
         model->SetLineColor(color);
         model->Draw("L");
       //legend->AddEntry(model, "GENIE", "L"); 
         int n = model->GetN();
         double xl = (model->GetX())[n-1];
         double yl = (model->GetY())[n-1];
         TLatex * t = new TLatex();
         t->SetTextColor(color);
         t->SetTextFont(42);
         t->SetTextSize(0.02);
         t->DrawLatex(xl+25,yl,Form("x=%.4f, F2 #times %.1f",xBj,scale));
      }
    }
    gC->Update();
  }

  // > xF3
  {
    gPS-> NewPage();
    gC -> SetFillColor(0);
    gC -> SetBorderMode(0);
    gC -> SetLogx(1);
    gC -> SetLogy(1);
    TH1F * hframe = 0;
    double xmin =  kQ2min;
    double xmax =  kQ2max*12;
    double ymin =  1E-2;
    double ymax =  2E+4;
    hframe = (TH1F*) gC->DrawFrame(xmin, ymin, xmax, ymax);
    hframe->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    hframe->GetYaxis()->SetTitle("xF3");
    // loop over data sets and plot data and corresponding GENIE predictions
    for(unsigned int icomp = 0; icomp < gComparisons.size(); icomp++) 
    {
      LOG("gvldtest", pNOTICE) 
         << "Producing plots for: " << gComparisons[icomp]->Label();
      int color = SummaryPlotColor[(icomp%kNumOfSummaryPlotColors)];
      double xBj = gComparisons[icomp]->xBjorken();
      string Fi  = gComparisons[icomp]->Fi();
      if(Fi != "xF3") continue;
      double scale = 1.0;
      if(xBj > 0.000 && xBj <= 0.008) scale = 20000.0;
      if(xBj > 0.008 && xBj <= 0.013) scale =  6000.0;
      if(xBj > 0.013 && xBj <= 0.016) scale =  3500.0;
      if(xBj > 0.016 && xBj <= 0.020) scale =  1500.0;
      if(xBj > 0.020 && xBj <= 0.030) scale =  1000.0;
      if(xBj > 0.030 && xBj <= 0.040) scale =   600.0;
      if(xBj > 0.040 && xBj <= 0.049) scale =   400.0;
      if(xBj > 0.049 && xBj <= 0.060) scale =   300.0;
      if(xBj > 0.060 && xBj <= 0.075) scale =   200.0;
      if(xBj > 0.075 && xBj <= 0.085) scale =   150.0;
      if(xBj > 0.085 && xBj <= 0.095) scale =   100.0;
      if(xBj > 0.095 && xBj <= 0.115) scale =    50.0;
      if(xBj > 0.115 && xBj <= 0.130) scale =    22.0;
      if(xBj > 0.130 && xBj <= 0.150) scale =    10.0;
      if(xBj > 0.150 && xBj <= 0.177) scale =     5.0;
      if(xBj > 0.177 && xBj <= 0.182) scale =     3.0;
      if(xBj > 0.182 && xBj <= 0.230) scale =     1.8;
      if(xBj > 0.230 && xBj <= 0.280) scale =     1.2;
      // Get expt data from GENIE archive
      vector<TGraphErrors *> data = Data(icomp,scale);
      if(data.size() == 0) continue;
      // Get the corresponding GENIE model prediction
      TGraph* model = Model(icomp,0,scale);
      for(unsigned int i = 0; i< data.size(); i++) {
         if(!data[i]) continue;
         data[i]->SetMarkerColor(color);
         data[i]->SetLineColor(color);
         data[i]->Draw("P");
       //legend->AddEntry(data[i], data[i]->GetTitle(), "P");
      }  
      if(model) {
         model->SetLineColor(color);
         model->Draw("L");
       //legend->AddEntry(model, "GENIE", "L"); 
         int n = model->GetN();
         double xl = (model->GetX())[n-1];
         double yl = (model->GetY())[n-1];
         TLatex * t = new TLatex();
         t->SetTextColor(color);
         t->SetTextFont(42);
         t->SetTextSize(0.02);
         t->DrawLatex(xl+25,yl,Form("x=%.4f, xF3 #times %.1f",xBj,scale));
      }
    }
    gC->Update();
  }

  //
  // Now also show each comparison on a separate plot
  //
  for(unsigned int icomp = 0; icomp < gComparisons.size(); icomp++) 
  {
    LOG("gvldtest", pNOTICE) 
        << "Producing plots for: " << gComparisons[icomp]->Label();
    double xBj = gComparisons[icomp]->xBjorken();
    string Fi  = gComparisons[icomp]->Fi();
    // Get expt data from GENIE archive
    vector<TGraphErrors *> data = Data(icomp);
    if(data.size() == 0) return;
    // Get the corresponding GENIE model prediction
    TGraph* model = Model(icomp,0);
    // Create new ps page if needed and get pad
    int plots_per_page = kNCx * kNCy;
    int iplot = 1 + icomp % plots_per_page;
    if(iplot == 1) {
       gPS-> NewPage();
       gC -> Clear();
       gC -> Divide(kNCx,kNCy);
    }
    gC -> GetPad(iplot) -> Range(0,0,100,100);
    gC -> GetPad(iplot) -> SetFillColor(0);
    gC -> GetPad(iplot) -> SetBorderMode(0);
    gC -> GetPad(iplot) -> SetLogx(1);
    gC -> GetPad(iplot) -> SetLogy(1);
    gC -> GetPad(iplot) -> cd();
    // Draw frame
    TH1F * hframe = 0;
    double xmin =  9999999;
    double xmax = -9999999;
    double ymin =  9999999;
    double ymax = -9999999;
    for(unsigned int i = 0; i< data.size(); i++) {
       if(!data[i]) continue;
       xmin  = TMath::Min(xmin, (data[i]->GetX())[TMath::LocMin(data[i]->GetN(),data[i]->GetX())]);
       xmax  = TMath::Max(xmax, (data[i]->GetX())[TMath::LocMax(data[i]->GetN(),data[i]->GetX())]);
       ymin  = TMath::Min(ymin, (data[i]->GetY())[TMath::LocMin(data[i]->GetN(),data[i]->GetY())]);
       ymax  = TMath::Max(ymax, (data[i]->GetY())[TMath::LocMax(data[i]->GetN(),data[i]->GetY())] );
    }
    hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(0.5*xmin, 0.3*ymin, 1.5*xmax, 5.0*ymax);
    hframe->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    hframe->GetYaxis()->SetTitle(Fi.c_str());
    hframe->Draw();
    // Draw data set & GENIE model and add legend
    TLegend * legend = new TLegend(0.20, 0.65, 0.55, 0.90);
    legend->SetLineStyle(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetHeader(Form("x_{Bj}=%.4f",xBj));
    for(unsigned int i = 0; i< data.size(); i++) {
       if(!data[i]) continue;
       data[i]->Draw("P");
       legend->AddEntry(data[i], data[i]->GetTitle(), "P");
    }  
    if(model) {
       model->Draw("L");
       legend->AddEntry(model, "GENIE", "L"); 
    }
    legend->Draw();
    gC->GetPad(iplot)->Update();
    gC->Update();
  }

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
     gOptCompInfoFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefCompInfoFilename;
        gOptCompInfoFilename = filename;
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
  if(parser.OptionExists("dis-xsec-model")){
     gOptDISModelName = parser.Arg("dis-xsec-model");
  }
  if(parser.OptionExists("dis-charm-xsec-model")){
     gOptDISCharmModelName = parser.Arg("dis-charm-xsec-model");
  }

}
//_________________________________________________________________________________
void PrintSyntax(void)
{

}
//_________________________________________________________________________________
