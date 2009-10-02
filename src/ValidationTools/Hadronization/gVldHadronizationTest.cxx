//____________________________________________________________________________
/*!

\program gvld_hadronization_test

\brief   Hadronization validation program
         (comparing GENIE models with neutrino bubble chamber data).

\syntax  gvld_hadronization_test -d mc_data_files

         where 'mc_data_files' is a string specifying a list of models and
         corresponding GHEP event files.
         The general syntax is:
         model_name:file_name.root,another_file_name.root,...+another_model_name:file_name.root,...

         Example:
         Assume that you have 2 hadronization models: `AGKY_v1' and `AGKY_v2'.
         You used the first one to generate the GHEP event files v1_1.root,
         v1_2.root and v1_3.root. Then you used the second one to generate
         the GHEP event files v2_1.root and v2_2.root. To process these files 
         and compare both models against bubble chamber data type:
         shell% gvld_hadronization_test -d \
            AGKY_v1:v1_1.root,v1_2.root,v1_3.root+AGKY_v2:v2_1.root,v2_2.root
         
\author  Tingjun Yang

\created March 1, 2009

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <map>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>

#include "Messenger/Messenger.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
#include "ValidationTools/Hadronization/HadPlots.h"
#include "ValidationTools/Hadronization/HadPlotter.h"

using namespace std;
using namespace genie;
using namespace genie::utils;
using namespace genie::vld_hadronization;

void GetCommandLineArgs (int argc, char** argv);
void PrintSyntax        (void);
void SetStyle           (bool bw=false);

map<string, vector<string> > gOptInputFiles; // model name -> list of files

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // get the command line arguments
  GetCommandLineArgs(argc, argv);

  // set style
  SetStyle();

  vector<HadPlots *> hpv;
  vector<HadPlots *>::iterator hpvit;

  // loop over models considered at the current validatio test
  for(map<string, vector<string> >::const_iterator mapit = 
          gOptInputFiles.begin(); mapit!=gOptInputFiles.end(); ++mapit) {

      // book plots for current model
      string model_name = mapit->first;
      LOG("VldHadro", pNOTICE)
         << "Booking plots for model: " << model_name;

      HadPlots * hadplot = new HadPlots(model_name);

      // include all data files for current model
      vector<string> filevec = mapit->second;
      for(vector<string>::const_iterator vecit=filevec.begin(); vecit!=filevec.end(); ++vecit) {
         string file_name = *vecit;
         LOG("VldHadro", pNOTICE)
            << " Loading data from file:.....: " << file_name;
         hadplot->LoadData(file_name);
      }

      // store
      hpv.push_back(hadplot);
  }

  // loop over models and analyze the input data files
  for(hpvit = hpv.begin(); hpvit != hpv.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    curr_model->Analyze();
  }

  // plot all bubble chamber data and start superimposing model predictions
  HadPlotter plotter;
  for(hpvit = hpv.begin(); hpvit != hpv.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    plotter.AddPlots(*curr_model);
  }
  plotter.ShowPlots();

  LOG("VldHadro", pNOTICE) << "Done!";

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char** argv)
{
  // help?
  bool help = genie::utils::clap::CmdLineArgAsBool(argc,argv,'h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  try {
    string input_data = genie::utils::clap::CmdLineArgAsString(argc,argv,'d');
    LOG("VldHadro", pNOTICE) << input_data;

    // split input string per model
    vector<string> mvec = utils::str::Split(input_data, "+");
    vector<string>::iterator mvecit = mvec.begin();
    for( ; mvecit != mvec.end(); ++mvecit) {
       string mcurr = *mvecit;

       // split into model name and list of files
       vector<string> mcurvec = utils::str::Split(mcurr, ":");

       assert(mcurvec.size()==2);
       string model_name = mcurvec[0];
       string file_list  = mcurvec[1];

       // split file list
       vector<string> filevec = utils::str::Split(file_list, ",");

       gOptInputFiles.insert(
         map<string, vector<string> >::value_type(model_name, filevec)
       );

    } // input models

  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
       PrintSyntax();
       exit(1);
    }
  }

  //
  // Report
  //
  map<string, vector<string> >::const_iterator mapit;
  vector<string>::const_iterator vecit;

  for(mapit=gOptInputFiles.begin(); mapit!=gOptInputFiles.end(); ++mapit) {
      string model_name = mapit->first;
      vector<string> filevec = mapit->second;
      LOG("VldHadro", pDEBUG)
        << " - Model name: " << model_name;
      for(vecit=filevec.begin(); vecit!=filevec.end(); ++vecit) {
         LOG("VldHadro", pDEBUG)
           << " |---> Including MC data file ....... " << *vecit;
      }
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("VldHadro", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gvld_hadronization_test -d model_name/file1.root,file2.root,...;another_model_name/file1.root,...\n";
}
//____________________________________________________________________________
void SetStyle(bool bw)
{
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Turn off all borders
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetDrawBorder(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetTitleBorderSize(0);

  // Set the size of the default canvas
  gStyle->SetCanvasDefH(600);
  gStyle->SetCanvasDefW(730);
  gStyle->SetCanvasDefX(10);
  gStyle->SetCanvasDefY(10);

  //set marker style
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1);

  // Set Line Widths
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFuncWidth(2);
  gStyle->SetHistLineWidth(3);
  gStyle->SetFuncColor(2);
  gStyle->SetFuncWidth(3);

  if(bw){
    gStyle->SetFuncWidth(1);
    gStyle->SetHistLineWidth(1);
    gStyle->SetFuncColor(1);
    gStyle->SetFuncWidth(1);
  }   

  // Set margins -- I like to shift the plot a little up and to the
  // right to make more room for axis labels
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadBottomMargin(0.20);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.03);

  // Set tick marks and turn off grids
  gStyle->SetNdivisions(505,"xyz");

  // Set Data/Stat/... and other options
  gStyle->SetOptDate(0);
//gStyle->SetDateX(0.1);
//gStyle->SetDateY(0.1);
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatFormat("6.2f");
  gStyle->SetFitFormat("8.4f");
  gStyle->SetOptFit(1);
  gStyle->SetStatH(0.20);
  gStyle->SetStatStyle(0);
  gStyle->SetStatW(0.30);
  gStyle->SetStatX(0.845);
  gStyle->SetStatY(0.845);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleX(0.15);
  gStyle->SetTitleW(0.75);

  // Adjust size and placement of axis labels
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.005,"x");
  gStyle->SetLabelOffset(0.005,"y");
  gStyle->SetLabelOffset(0.005,"z");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.2,"xz");
  gStyle->SetTitleOffset(1,"y");
  // Set paper size for life in the US
  gStyle->SetPaperSize(TStyle::kUSLetter);

  gStyle->SetTitleY(.90);
  gStyle->SetPalette(1);

  if(bw){
    const int ncol=7;
    double red[ncol];
    double green[ncol];
    double blue[ncol];
    double stops[ncol];
    double dcol = -1/double(ncol);
    double gray = 1;
    for (int j = 0; j < ncol; j++) {
      //   ...... Define color with RGB equal to : gray, gray, gray
      stops[j]=double(j)/double(ncol-1);
      red[j]=gray;
      blue[j]=gray;
      green[j]=gray;
      
      gray += dcol;
    }
    UInt_t totcol=50;
    TColor::CreateGradientColorTable(ncol,stops,red,green,blue,totcol); 
  }
  gStyle->SetLegendBorderSize(0);
  gROOT->ForceStyle();
}
//____________________________________________________________________________
