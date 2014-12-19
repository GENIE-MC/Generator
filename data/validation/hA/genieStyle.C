#include <TStyle.h>
#include <TLatex.h>


void set_root_env(){

  //TStyle* genieStyle = new TStyle("genieStyle", "GENIE Style");

  //set the background color to white
 gStyle->SetFillColor(10);
 gStyle->SetFrameFillColor(10);
 gStyle->SetCanvasColor(10);
 gStyle->SetPadColor(10);
 gStyle->SetTitleFillColor(0);
 gStyle->SetStatColor(10);

//dont put a colored frame around the plots
 gStyle->SetFrameBorderMode(0);
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetPadBorderMode(0);
 gStyle->SetLegendBorderSize(3);

//use the primary color palette
 gStyle->SetPalette(1,0);

//set the default line color for a histogram to be black
 gStyle->SetHistLineColor(kBlack);

//set the default line color for a fit function to be red
 gStyle->SetFuncColor(kRed);

//make the axis labels black
 gStyle->SetLabelColor(kBlack,"xyz");

//set the default title color to be black
 gStyle->SetTitleColor(kBlack);
 
//set the margins
 gStyle->SetPadBottomMargin(0.17);
 gStyle->SetPadTopMargin(0.11);
 gStyle->SetPadRightMargin(0.08);
 gStyle->SetPadLeftMargin(0.17);

//set axis label and title text sizes
 gStyle->SetLabelFont(42,"xyz");
 gStyle->SetLabelSize(0.06,"xyz");
 gStyle->SetLabelOffset(0.015,"xyz");
 gStyle->SetTitleFont(42,"xyz");
 gStyle->SetTitleSize(0.05,"xyz");
 gStyle->SetTitleOffset(1.4,"y");
 gStyle->SetTitleOffset(1.3,"x");
 gStyle->SetTitleX(0.27);
 gStyle->SetStatFont(42);
 gStyle->SetStatFontSize(0.07);
 gStyle->SetTitleBorderSize(1);
 gStyle->SetStatBorderSize(0);
 gStyle->SetTextFont(42);
 gStyle->SetTitleW(0.5);
 gStyle->SetTitleH(0.1);

//set line widths
 gStyle->SetFrameLineWidth(2);
 gStyle->SetFuncWidth(2);
 gStyle->SetHistLineWidth(2);

//set the number of divisions to show
 gStyle->SetNdivisions(506, "xy");
 //gStyle->SetPadTickX(-50202);

//turn off xy grids
 gStyle->SetPadGridX(0);
 gStyle->SetPadGridY(0);

//set the tick mark style
 gStyle->SetPadTickX(1);
 gStyle->SetPadTickY(1);

//turn off stats
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(0);

//marker/line settings
//gStyle->SetMarkerStyle(20);
 gStyle->SetMarkerSize(.95);//0.7
 gStyle->SetLineWidth(2); 
 gStyle->SetErrorX(0);
 gStyle->SetHistLineStyle(3);

//done
 gStyle->cd();
 gROOT->ForceStyle();
 gStyle->ls();

}

void add_plot_label( char* label, double x, double y, double size = 0.05, int color = 1, int font = 62, int align = 22 ){

  TLatex *latex = new TLatex( x, y, label );
  latex->SetNDC();
  latex->SetTextSize(size);
  latex->SetTextColor(color);
  latex->SetTextFont(font);
  latex->SetTextAlign(align);
  latex->Draw();

}


