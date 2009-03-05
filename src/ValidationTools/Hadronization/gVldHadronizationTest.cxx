//____________________________________________________________________________
/*!

\program gvld_hadronization_test

\brief   Hadronization test (validating GENIE hadronization models with
         neutrino bubble chamber data)

\syntax  ./gvld_hadronization_test

\author  Tingjun Yang

\created March 1, 2009

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

#include "ValidationTools/Hadronization/HadPlots.h"
#include "ValidationTools/Hadronization/HadPlotter.h"

using namespace std;
using namespace genie;
using namespace genie::vld;

void SetStyle(bool bw=false);

int main()
{
  SetStyle();

  HadPlots agky("AGKY");
  agky.LoadData("../../../macros/gntp.0.ghep.root");
  agky.LoadData("../../../macros/gntp.1.ghep.root");
  agky.LoadData("../../../macros/gntp.2.ghep.root");
  agky.LoadData("../../../macros/gntp.3.ghep.root");
  agky.Analyze();

  HadPlotter plotter;;
  plotter.AddPlots(agky);
  plotter.ShowPlots();

  return 0;

}

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
   TColor::CreateGradientColorTable(ncol,stops,
				    red,green,blue,totcol);
   
}
gStyle->SetLegendBorderSize(0);
gROOT->ForceStyle();
}
