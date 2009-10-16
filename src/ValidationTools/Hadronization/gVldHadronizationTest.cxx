//____________________________________________________________________________
/*!

\program gvld_hadronization_test

\brief   Hadronization validation program
         (comparing GENIE models with neutrino bubble chamber data).

\syntax  gvld_hadronization_test -g genie_inputs.xml [-f fmt]

         Options:

          [] Denotes an optional argument.
          
          -f Output plot format (0: eps, 1: gif) [default: gif]          
          -g An input XML file for specifying a GHEP event list for each
             model to be considered in the hadronization benchmark test.

             The input file should look like:

             <?xml version="1.0" encoding="ISO-8859-1"?>
             <vld_inputs>
               <model name="a_model_name">
                 <evt_file format="ghep"> /model_1/evtfile0.root </evt_file>
                 <evt_file format="ghep"> /model_1/evtfile1.root </evt_file>
                 <evt_file format="ghep"> /model_1/evtfile2.root </evt_file>
                 ...
               </model>

               <model name="another_model_name">
                 <evt_file format="ghep"> /model_2/evtfile0.root </evt_file>
                 <evt_file format="ghep"> /model_2/evtfile1.root </evt_file>
                 <evt_file format="ghep"> /model_2/evtfile2.root </evt_file>
                 ...
               </model>
               ...
             </vld_inputs>


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
#include "Utils/VldTestInputs.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
#include "ValidationTools/Hadronization/HadPlots.h"
#include "ValidationTools/Hadronization/HadPlotter.h"

using namespace std;
using namespace genie;
using namespace genie::utils;
using namespace genie::utils::vld;
using namespace genie::vld_hadronization;

// prototypes
void SetStyle              (bool bw=false);
void LoadFilesAndBookPlots (void);
void Analyze               (void);
void Plot                  (void);
void End                   (void);
void GetCommandLineArgs    (int argc, char** argv);
void PrintSyntax           (void);

// globals & user inputs
vector<HadPlots *> gHP;
VldTestInputs      gOptGenieInputs(false,10);;
int                gFmt = 1;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc, argv);
  SetStyle();
  LoadFilesAndBookPlots();
  Analyze();
  Plot();
  End();

  LOG("vldtest", pNOTICE) << "Done!";

  return 0;
}
//____________________________________________________________________________
void SetStyle(bool bw)
{
  gROOT->SetStyle("Plain");

  gStyle -> SetPadTickX (1);
  gStyle -> SetPadTickY (1);

  // Turn off all borders
  //
  gStyle -> SetCanvasBorderMode (0);
  gStyle -> SetFrameBorderMode  (0);
  gStyle -> SetPadBorderMode    (0);
  gStyle -> SetDrawBorder       (0);
  gStyle -> SetCanvasBorderSize (0);
  gStyle -> SetFrameBorderSize  (0);
  gStyle -> SetPadBorderSize    (0);
  gStyle -> SetTitleBorderSize  (0);

  // Set the size of the default canvas
  //
  gStyle -> SetCanvasDefH (600);
  gStyle -> SetCanvasDefW (730);
  gStyle -> SetCanvasDefX  (10);
  gStyle -> SetCanvasDefY  (10);

  // Set marker style
  //
  gStyle -> SetMarkerStyle (20);
  gStyle -> SetMarkerSize   (1);

  // Set line widths
  //
  gStyle -> SetFrameLineWidth (1);
  gStyle -> SetFuncWidth      (2);
  gStyle -> SetHistLineWidth  (3);
  gStyle -> SetFuncColor      (2);
  gStyle -> SetFuncWidth      (3);

  // Set margins 
  //
  gStyle -> SetPadTopMargin    (0.10);
  gStyle -> SetPadBottomMargin (0.20);
  gStyle -> SetPadLeftMargin   (0.15);
  gStyle -> SetPadRightMargin  (0.03);

  // Set tick marks and turn off grids
  //
  gStyle -> SetNdivisions (505,"xyz");

  // Adjust size and placement of axis labels
  //
  gStyle -> SetLabelSize   (0.050,  "xyz");
  gStyle -> SetLabelOffset (0.005,  "x"  );
  gStyle -> SetLabelOffset (0.005,  "y"  );
  gStyle -> SetLabelOffset (0.005,  "z"  );
  gStyle -> SetTitleSize   (0.060,  "xyz");
  gStyle -> SetTitleOffset (1.200,  "xz" );
  gStyle -> SetTitleOffset (1.000,  "y"  );

  // Set Data/Stat/... and other options
  //
  gStyle -> SetOptDate          (0);
  gStyle -> SetOptFile          (0);
  gStyle -> SetOptStat          (0);
  gStyle -> SetStatFormat       ("6.2f");
  gStyle -> SetFitFormat        ("8.4f");
  gStyle -> SetOptFit           (1);
  gStyle -> SetStatH            (0.20);
  gStyle -> SetStatStyle        (0);
  gStyle -> SetStatW            (0.30);
  gStyle -> SetStatX            (0.845);
  gStyle -> SetStatY            (0.845);
  gStyle -> SetOptTitle         (0);
  gStyle -> SetTitleX           (0.15);
  gStyle -> SetTitleW           (0.75);
  gStyle -> SetTitleY           (0.90);
  gStyle -> SetPalette          (1);
  gStyle -> SetLegendBorderSize (0);


  // Set paper size for life in the US or EU
  //
  gStyle -> SetPaperSize (TStyle::kA4);       //<-- tartes aux fraises
//gStyle -> SetPaperSize (TStyle::kUSLetter); //<-- donuts

  // In B&W (papers)
  //
  if(bw){
    const int ncol = 7;

    double red   [ncol];
    double green [ncol];
    double blue  [ncol];
    double stops [ncol];

    double dcol = -1/double(ncol);
    double gray = 1;
    for (int j = 0; j < ncol; j++) {
      // Define color with RGB equal to : gray, gray, gray
      stops[j] = double(j)/double(ncol-1);
      red  [j] = gray;
      blue [j] = gray;
      green[j] = gray;
      
      gray += dcol;
    }
    UInt_t totcol=50;
    TColor::CreateGradientColorTable(ncol,stops,red,green,blue,totcol); 

    gStyle -> SetFuncWidth     (1);
    gStyle -> SetHistLineWidth (1);
    gStyle -> SetFuncColor     (1);
    gStyle -> SetFuncWidth     (1);
  }//bw

  gROOT->ForceStyle();
}
//____________________________________________________________________________
void LoadFilesAndBookPlots(void)
{
  vector<HadPlots *>::iterator hpvit;

  // loop over models considered at the current validation test
  for(int imodel=0; imodel < gOptGenieInputs.NModels(); imodel++) {

    string model_name = gOptGenieInputs.ModelTag(imodel);
    vector<string> & event_filenames = gOptGenieInputs.EvtFileNames(imodel);

      LOG("vldtest", pNOTICE)
            << "Booking plots for model: " << model_name;

      HadPlots * hadplot = new HadPlots(model_name);

      // include all data files for current model
      vector<string>::const_iterator file_iter = event_filenames.begin();
      for( ; file_iter != event_filenames.end(); ++file_iter) {
         string filename = *file_iter;
         LOG("vldtest", pNOTICE)
            << " Loading data from file:.....: " << filename;
         hadplot->LoadData(filename);
      }// file_iter

      // store
      gHP.push_back(hadplot);
  }
}
//____________________________________________________________________________
void Analyze(void)
{
  vector<HadPlots *>::iterator hpvit = gHP.begin();
  for( ; hpvit != gHP.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    curr_model->Analyze();
  }
}
//____________________________________________________________________________
void Plot(void)
{
  // plot all bubble chamber data and start superimposing model predictions
  bool in_eps = (gFmt==0);
  HadPlotter plotter(in_eps);
  vector<HadPlots *>::iterator hpvit = gHP.begin();
  for( ; hpvit != gHP.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    plotter.AddPlots(*curr_model);
  }
  plotter.ShowPlots();
}
//____________________________________________________________________________
void End(void)
{
  vector<HadPlots *>::iterator hpvit = gHP.begin();
  for( ; hpvit != gHP.end(); ++hpvit) {
    HadPlots * curr_model = *hpvit;
    if(curr_model) {
       delete curr_model;
       curr_model = 0;
    }
  }
  gHP.clear();
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

  // get GENIE inputs
  try {
     string inputs = utils::clap::CmdLineArgAsString(argc,argv,'g');
     bool ok = gOptGenieInputs.LoadFromFile(inputs);
     if(!ok) { 
        LOG("gvldtest", pFATAL) << "Could not read: " << inputs;
        exit(1);
     }
  } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
     }
  }

  // output plot format
  try {
     int format = utils::clap::CmdLineArgAsInt(argc,argv,'f');
     if(format==0 || format==1) {
       gFmt=format;
     }
  } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
     }
  }

  if(gOptGenieInputs.NModels()==0) {
    LOG("gvldtest", pFATAL) << "Not input model data to analyze";
    exit(1);
  }

  LOG("gvldtest", pFATAL) << "Input data: ";  
  LOG("gvldtest", pFATAL) << gOptGenieInputs;

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("vldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gvld_hadronization_test -g genie_inputs.xml\n [-f format]";
}
//____________________________________________________________________________
