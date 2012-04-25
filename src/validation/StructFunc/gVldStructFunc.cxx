//____________________________________________________________________________
/*!

\program gvld_sf

\brief   Compares GENIE F2 and xF3 structure-function models with data.

         Syntax:
           gvld_sf 
                [-r model_for_resonance_xsec]
                [-c model_for_dis_continuum_xsec]
                [-d data_archive] 

         Options:
           [] Denotes an optional argument.

          -r Specify GENIE resonance cross-section model.
          -c Specify GENIE DIS cross-section model.

          -d Location of the neutrino cross-section data archive.
             By default, the program will look-up the one located in:
             $GENIE/data/validation/sf/??

         Example:

            % gvld_sf
                  -r genie::ReinSeghalRESPXSec/Default
                  -c genie::QPMDISPXSec/Default


		      
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

         Jelena Ilic <jelena.ilic \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 06, 2008 

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
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
  SFCompInfo(string label, string datasets) :
   fLabel       (label),
   fDataSetKeys (datasets)
  { 
  }
 ~SFCompInfo() { }

  string Label       (void) const { return fLabel;         }
  string DataSetKeys (void) const { return fDataSetKeys;   }

private:

  string fLabel;
  string fDataSetKeys;
};
//.................................................................................

//
// constants
//

// defaults
const char * kDefDataArchiveFilename = "data/validation/vA/sf/SF.root";  
const char * kDefDataSetsFilename    = "data/validation/vA/sf/datasets.txt";  

//
// globals
//

string gOptDataArchiveFilename = ""; // -d command-line argument
string gOptDataSetsFilename    = ""; // -s command-line argument
string gOptRESModelName        = ""; // -r command-line argument
string gOptDISModelName        = ""; // -c command-line argument

TFile *        gSFDataFile  = 0;
TTree *        gSFDataTree  = 0;
TPostScript *  gPS          = 0;
TCanvas *      gC           = 0;

const XSecAlgorithmI * gRESXSecModel = 0; // resonance cross-section model
const XSecAlgorithmI * gDISXSecModel = 0; // DIS cross-section model

vector<SFCompInfo *> gComparisons; // info on which comparisons to perform

StructFunc gStructFuncCalc; // utility class extracting structure functions from the cross-section model

// function prototypes
void     Init               (void);
void     End                (void);
void     Draw               (unsigned int icomp);
TH1F *   DrawFrame          (TGraph * gr0, TGraph * gr1, TCanvas * c);
void     GetCommandLineArgs (int argc, char ** argv);
void     PrintSyntax        (void);

vector<TGraph *>       Model(unsigned int icomp, unsigned int imodel);
vector<TGraphErrors *> Data (unsigned int icomp);

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
 GetCommandLineArgs (argc,argv);

  Init();

  // loop over data sets and plot data and corresponding GENIE predictions
  for(unsigned int icomp = 0; icomp < gComparisons.size(); icomp++) 
  {
    LOG("gvldtest", pNOTICE) 
      << "Producing plots for: " << gComparisons[icomp]->Label();
    Draw(icomp);
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
     gRESXSecModel =
        dynamic_cast<const XSecAlgorithmI *> (
            algf->GetAlgorithm(model_name, model_conf));
  }
  gDISXSecModel = 0;
  if(gOptDISModelName != "none"){
     vector<string> modelv = utils::str::Split(gOptDISModelName,"/");
     assert(modelv.size()==2);
     string model_name = modelv[0];
     string model_conf = modelv[1];
     gDISXSecModel =
        dynamic_cast<const XSecAlgorithmI *> (
            algf->GetAlgorithm(model_name, model_conf));
  }

  gStructFuncCalc.SetResonanceXSecModel (gRESXSecModel);
  gStructFuncCalc.SetDISXSecModel       (gDISXSecModel);
  gStructFuncCalc.SetDISCharmXSecModel  (0);

  // Create plot canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

 // Create output postscript file
  string localtime = utils::system::LocalTimeAsString("%d.%d.%d_%d.%d.%d"); 
  string filename  = Form("genie-sf_data_comp-%s.ps",localtime.c_str());
  gPS = new TPostScript(filename.c_str(), 111);

  // Add cover page
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE comparison with F2, xF3 structure function data");
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
vector<TGraph *> Model(unsigned int icomp, unsigned int imodel)
{
  vector<TGraph *> models;

  return models;
}
//_________________________________________________________________________________
vector<TGraphErrors *> Data(unsigned int icomp)
{
  vector<TGraphErrors *> data;

  return data;
}
//_________________________________________________________________________________
void Draw(unsigned int icomp)
{
  bool inlogy = false; // log y scale? Eventually get it from SFCompInfo 

  // get all measurements for the current channel from the NuValidator MySQL dbase
  vector<TGraphErrors *> data = Data(icomp);
  if(data.size() == 0) return;

  // get the corresponding GENIE model prediction
  vector<TGraph*> models = Model(icomp,0);

  gPS->NewPage();

  gC->Clear();
  gC->Divide(2,1);
  gC->GetPad(1)->SetPad("mplots_pad","",0.01,0.25,0.99,0.99);
  gC->GetPad(2)->SetPad("legend_pad","",0.01,0.01,0.99,0.24);
  gC->GetPad(1)->SetFillColor(0);
  gC->GetPad(1)->SetBorderMode(0);
  gC->GetPad(2)->SetFillColor(0);
  gC->GetPad(2)->SetBorderMode(0);
  gC->GetPad(1)->cd();
  gC->GetPad(1)->SetBorderMode(0);
  gC->GetPad(1)->SetLogx(1);
  gC->GetPad(1)->SetLogy(1);

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
  double ymax_scale = (inlogy) ? 2. : 1.4;
  hframe = (TH1F*) gC->GetPad(1)->DrawFrame(0.5*xmin, 0.4*ymin, 1.2*xmax, ymax_scale*ymax);
  hframe->GetXaxis()->SetTitle("???");
  hframe->GetYaxis()->SetTitle("???");
  hframe->Draw();

  // add legend
  TLegend * legend = new TLegend(0.01, 0.01, 0.99, 0.99);
  legend->SetLineStyle(0);
  legend->SetFillColor(0);
  legend->SetTextSize(0.06);

  // draw current data set
  for(unsigned int i = 0; i< data.size(); i++) {
    if(!data[i]) continue;
    data[i]->Draw("P");
    legend->AddEntry(data[i], data[i]->GetTitle(), "LP");
  }  

  // have model prediction to plot?
  // just one for the time-being - superimposing different models will be added asap
  if(models.size() == 1) {
       models[0]->Draw("L");
       legend->AddEntry(models[0], "GENIE", "L"); 
  }

  gC->GetPad(2)->cd();
  legend->Draw();

  gC->GetPad(2)->Update();
  gC->Update();
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
  if(parser.OptionExists('r')){
     gOptRESModelName = parser.ArgAsString('r');
  }   
  if(parser.OptionExists('c')){
     gOptDISModelName = parser.ArgAsString('c');
  }

}
//_________________________________________________________________________________
void PrintSyntax(void)
{

}
//_________________________________________________________________________________
