//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - Aug 26, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Dec 08, 2008 - CA
   NuValidator package refurbishment. Removed the neugen3 dependency. 
   Moved all sources to $GENIE/src/ValidationTools/NuVld.
   Some clases have been renamed.
*/
//____________________________________________________________________________ 

#include <TRootEmbeddedCanvas.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGProgressBar.h>
#include <TMath.h>
#include <TH1F.h>

#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/DBNuXSecTableRow.h"
#include "ValidationTools/NuVld/DBElDiffXSecTableRow.h"
#include "ValidationTools/NuVld/DBSFTableRow.h"
#include "ValidationTools/NuVld/GuiTableRenderer.h"
#include "ValidationTools/NuVld/GuiSysLogSingleton.h"
#include "ValidationTools/NuVld/GuiBrowserSingleton.h"
#include "ValidationTools/NuVld/NuVldUserData.h"

using namespace genie;
using namespace genie::nuvld;

ClassImp(GuiTableRenderer)

//______________________________________________________________________________
GuiTableRenderer::GuiTableRenderer() :
fEmbCanvas(0)
{
  fXmin = -1;
  fXmax = -1;

  fPlotVar  = "";
  fDrawOpt  = "";
  fErrorOpt = "";
  
  fLegend  = 0;
  
  fScaleWithE       = false;
  fIsMultigraph     = false;
  fHaveCustomXRange = false;
}
//______________________________________________________________________________
GuiTableRenderer::GuiTableRenderer(TRootEmbeddedCanvas * emb_canvas) :
fEmbCanvas(emb_canvas)
{
  fXmin = -1;
  fXmax = -1;

  fPlotVar  = "";
  fDrawOpt  = "";
  fErrorOpt = "";

  fLegend  = 0;

  fScaleWithE       = false;
  fIsMultigraph     = false;
  fHaveCustomXRange = false;
}
//______________________________________________________________________________
GuiTableRenderer::~GuiTableRenderer()
{

}
//______________________________________________________________________________
void GuiTableRenderer::SwitchEmbdeddedCanvas(TRootEmbeddedCanvas * emb_canvas)
{
  fEmbCanvas = emb_canvas;
}
//______________________________________________________________________________
void GuiTableRenderer::SetScaleWithEnergy(bool tf)
{
  fScaleWithE = tf;
}
//______________________________________________________________________________
void GuiTableRenderer::SetMultigraph(bool tf)
{
  fIsMultigraph = tf;
}
//______________________________________________________________________________
void GuiTableRenderer::SetCustomXRange(double xmin, double xmax)
{
  fXmin = xmin;
  fXmax = xmax;

  fHaveCustomXRange = true;
}
//______________________________________________________________________________
void GuiTableRenderer::ResetCustomXRange(void)
{
  fXmin = -1;
  fXmax = -1;

  fHaveCustomXRange = false;
}
//______________________________________________________________________________
void GuiTableRenderer::SetPlotVariable(string plot_var)
{
  fPlotVar = plot_var;
}
//______________________________________________________________________________
void GuiTableRenderer::SetErrorOption(string err_opt)
{
  fErrorOpt = err_opt;
}
//______________________________________________________________________________
void GuiTableRenderer::SetDrawOption(string draw_opt)
{
  fDrawOpt = draw_opt;
}
//______________________________________________________________________________
void GuiTableRenderer::SetExternalLegend(TLegend * legend)
{
  fLegend = legend;
}
//______________________________________________________________________________
void GuiTableRenderer::PrintDrawingOptions(void)
{
  LOG("NuVld", pDEBUG)
     << "GuiTableRenderer options: \n"
     << "Embedded canvas ....... " << ( (fEmbCanvas) ? fEmbCanvas->GetName() : "" ) << "\n"
     << "Scale with energy ..... " << fScaleWithE << "\n"
     << "Set multigraph ........ " << fIsMultigraph << "\n"
     << "Custom x range ........ " << fXmin << ", " << fXmax << "\n"
     << "Plot variable ......... " << fPlotVar << "\n"
     << "Error options ......... " << fErrorOpt << "\n";     
}
//______________________________________________________________________________
void GuiTableRenderer::DrawXSecTable(DBTable<DBNuXSecTableRow> * table)
{
  GuiSysLogSingleton *  syslog  = GuiSysLogSingleton::Instance();
  
  fEmbCanvas->GetCanvas()->cd();

  // if >0 entries plot them, else display warning
  if( table->NRows() > 0) {

    // Check if the plot should be color-coded...

    if(fIsMultigraph) {
       this->DrawMultiGraphInCanvas(table);

    } else {      
       this->DrawGraphInCanvas(table);
    } 

    fEmbCanvas->GetCanvas()->Update();

    syslog->ProgressBar()->SetPosition(100);

    syslog->StatusBar()->SetText( "Plot is shown in 'Plotter' tab", 0 );
    syslog->StatusBar()->SetText( " ",    1 );

    syslog->ProgressBar()->SetPosition(0);

  } else {
    
    syslog->StatusBar() -> SetText( "The table you want to draw is empty", 1 );
    syslog->Log()       -> AddLine( "The table you want to draw is empty"    );

    syslog->ProgressBar()->SetPosition(0);

  } // lines > 0
}
//______________________________________________________________________________
void GuiTableRenderer::DrawXSecTable(DBTable<DBElDiffXSecTableRow> * table)
{
  GuiSysLogSingleton *  syslog  = GuiSysLogSingleton::Instance();

  fEmbCanvas->GetCanvas()->cd();

  //x-variable

  // if >0 entries plot them, else display warning
  if( table->NRows() > 0) {

    // Check if the plot should be color-coded...

    if(fIsMultigraph) {
       this->DrawMultiGraphInCanvas(table);
    } else {
       this->DrawGraphInCanvas(table);
    } 

    fEmbCanvas->GetCanvas()->Update();

    syslog->ProgressBar()->SetPosition(100);

    syslog->StatusBar()->SetText( "Plot is shown in 'Plotter' tab", 0 );

    syslog->ProgressBar()->SetPosition(0);

  } else {

    syslog->StatusBar() -> SetText( "The table you want to draw is empty", 1 );
    syslog->Log()       -> AddLine( "The table you want to draw is empty"    );

    syslog->ProgressBar()->SetPosition(0);
    
  } // lines > 0
}
//______________________________________________________________________________
void GuiTableRenderer::DrawXSecTable(DBTable<DBSFTableRow> * table)
{
  GuiSysLogSingleton *  syslog  = GuiSysLogSingleton::Instance();

  fEmbCanvas->GetCanvas()->cd();

  // if >0 entries plot them, else display warning
  if( table->NRows() > 0) {

    // Check if the plot should be color-coded...

    if(fIsMultigraph) {
       this->DrawMultiGraphInCanvas(table);

    } else {
       this->DrawGraphInCanvas(table);
    }

    fEmbCanvas->GetCanvas()->Update();

    syslog->ProgressBar()->SetPosition(100);

    syslog->StatusBar()->SetText( "Plot is shown in 'Plotter' tab", 0 );
    syslog->StatusBar()->SetText( " ",    1 );

    syslog->ProgressBar()->SetPosition(0);

  } else {

    syslog->StatusBar() -> SetText( "The table you want to draw is empty", 1 );
    syslog->Log()       -> AddLine( "The table you want to draw is empty"    );

    syslog->ProgressBar()->SetPosition(0);

  } // lines > 0
}
//______________________________________________________________________________
void GuiTableRenderer::DrawGraphInCanvas(DBTable<DBNuXSecTableRow> * table)
{
  string draw_opt = fErrorOpt + ( (fScaleWithE) ? "-scale-with-E" : "") ;
  
  TGraphAsymmErrors * graph = table->GetGraph( draw_opt.c_str() );

  TCanvas * c = fEmbCanvas->GetCanvas();
  
  if(graph != 0) {
    double xmin = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    double xmax = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    double ymin = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    double ymax = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];

    // check if there is a user-defined energy Range
    if(fHaveCustomXRange) {
         xmin = fXmin;
         xmax = fXmax;
    }
    
    c->Clear();
    c->Divide(2,1);

    c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.99);
    c->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);
      
    c->GetPad(1)->SetFillColor(0);
    c->GetPad(1)->SetBorderMode(0);
    c->GetPad(2)->SetFillColor(0);
    c->GetPad(2)->SetBorderMode(0);

    c->GetPad(1)->cd();

    c->GetPad(1)->SetBorderMode(0);

    TH1F * hframe = (TH1F*) c->GetPad(1)->DrawFrame(.98*xmin, .98*ymin, 1.02*xmax, 1.02*ymax);
    hframe->Draw();

    graph->SetMarkerStyle(3);
    graph->SetMarkerSize(1);

    graph->Draw("P");

    c->GetPad(1)->cd();

    hframe->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    if(fScaleWithE)
          hframe->GetYaxis()->SetTitle("#sigma/E_{#nu} (10^{-38} cm^{2}/GeV/nucleon)");
    else  hframe->GetYaxis()->SetTitle("#sigma (10^{-38} cm^{2}/nucleon)");

    hframe->GetXaxis()->SetTitleFont(32);
    hframe->GetXaxis()->SetTitleColor(1);
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleFont(32);
    hframe->GetYaxis()->SetTitleColor(1);
    hframe->GetYaxis()->SetTitleSize(0.04);
    hframe->GetYaxis()->SetTitleOffset(1);
    hframe->GetYaxis()->SetLabelSize(0.03);

    // check if it is better to plot it in log scale
    if(xmin > 0 && xmax > xmin) {
       double Range = log10(xmax/xmin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogx();
       else                     c->GetPad(1)->SetLogx(0);
    }
    if(ymin > 0 && ymax > ymin) {
       double Range = log10(ymax/ymin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogy();
       else                     c->GetPad(1)->SetLogy(0);
    }
  }
  c->GetPad(1)->Update();
  c->Update();
}
//____________________________________________________________________________
void GuiTableRenderer::DrawMultiGraphInCanvas(DBTable<DBNuXSecTableRow> * table)
{
  string draw_opt = fErrorOpt + ( (fScaleWithE) ? "-scale-with-E" : "") ;

  MultiGraph * mgraph = table->GetMultiGraph( draw_opt.c_str() );

  TH1F * hframe = 0;

  double xmin = 0, xmax = 0, ymin = 0, ymax = 0;

  TGraphAsymmErrors * graph = table->GetGraph( draw_opt.c_str() );

  if(graph != 0) {
    xmin  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    xmax  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    ymin  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    ymax  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];

    delete graph;
  }
  
  TCanvas * c = fEmbCanvas->GetCanvas();

  if(mgraph != 0) {

    // check if there is a user-defined energy Range
    if(fHaveCustomXRange) {
         xmin = fXmin;
         xmax = fXmax;
    }

    SLOG("NuVld", pDEBUG)
                  << " x = [" << xmin << ", " << xmax << "],"
                  << " y = [" << ymin << ", " << ymax << "] ";

    if(fLegend) {

       // fill in the legend but do not plot it. Return the TLegend to be plotted in
       // an external 'key-canvas'. Since there is no legend, use the full canvas.

       c->Clear();
       c->Divide(2,1);

       c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.99);
       c->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);

       c->GetPad(1)->SetFillColor(0);
       c->GetPad(1)->SetBorderMode(0);
       c->GetPad(2)->SetFillColor(0);
       c->GetPad(2)->SetBorderMode(0);

       c->GetPad(1)->cd();

       c->GetPad(1)->SetBorderMode(0);
       
       hframe = (TH1F*) c->GetPad(1)->DrawFrame(.98*xmin, .98*ymin, 1.02*xmax, 1.02*ymax);
       hframe->Draw();

       for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {

          LOG("NuVld", pDEBUG)
                   << "** Drawing graph........." << igraph
                               << " with N-Points = " << mgraph->GetGraph(igraph)->GetN();
                  
          mgraph->GetGraph(igraph)->Draw("P");
       }
       
       c->GetPad(1)->Update();

       delete fLegend;
       fLegend = new TLegend(0.01, 0.01, 0.99, 0.99);
       fLegend->SetFillColor(0);
       mgraph->FillLegend("LP", fLegend);

       TCanvas * cext = new TCanvas("cext","Multi-Graph Legend",0,0,400,700);
       cext->cd();
       fLegend->SetTextSize(0.03);
       fLegend->Draw();
       cext->Update();
       
    }  else {
      
       // no input legend - make one and superimpose it on the MultiGraph plotted on
       // the TRootEmbeddedCanvas. Split the canvas to fit both the plot & the legend.

       c->Clear();
       c->Divide(2,1);

       c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.60,0.99);
       c->GetPad(2)->SetPad("legend_pad","",0.61,0.01,0.99,0.99);

       c->GetPad(1)->cd();

       c->GetPad(1)->SetFillColor(0);
       c->GetPad(1)->SetBorderMode(0);
       c->GetPad(2)->SetFillColor(0);
       c->GetPad(2)->SetBorderMode(0);

       hframe = (TH1F*) c->GetPad(1)->DrawFrame(.98*xmin, .98*ymin, 1.02*xmax, 1.02*ymax);
       hframe->Draw();

       for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {

            LOG("NuVld", pDEBUG)
                   << "** Drawing graph........." << igraph
                               << " with N-Points = " << mgraph->GetGraph(igraph)->GetN();

            mgraph->GetGraph(igraph)->Draw("P");
       }
       
       c->GetPad(1)->Update();

       SLOG("NuVld", pDEBUG) << "Drawing legend";

       c->GetPad(2)->cd();

       if(fLegend) delete fLegend;
       fLegend = new TLegend(0.01, 0.01, 0.99, 0.99);
       fLegend->SetFillColor(0);
       mgraph->FillLegend("LP", fLegend);
       fLegend->SetTextSize(0.03);
       fLegend->Draw();

       c->GetPad(2)->Update();
    }

    c->GetPad(1)->cd();

    hframe->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    if(fScaleWithE)
          hframe->GetYaxis()->SetTitle("#sigma/E_{#nu} (10^{-38} cm^{2}/GeV/nucleon)");
    else  hframe->GetYaxis()->SetTitle("#sigma (10^{-38}  cm^{2}/nucleon)");

    hframe->GetXaxis()->SetTitleFont(32);
    hframe->GetXaxis()->SetTitleColor(1);
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleFont(32);
    hframe->GetYaxis()->SetTitleColor(1);
    hframe->GetYaxis()->SetTitleSize(0.04);
    hframe->GetYaxis()->SetTitleOffset(1);
    hframe->GetYaxis()->SetLabelSize(0.03);

    // check if it is better to plot it in log scale
    if(xmin > 0 && xmax > xmin) {
       double Range = log10(xmax/xmin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogx();
       else                     c->GetPad(1)->SetLogx(0);
    }
    if(ymin > 0 && ymax > ymin) {
       double Range = log10(ymax/ymin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogy();
       else                     c->GetPad(1)->SetLogy(0);
    }

    c->GetPad(1)->Update();
    c->Update();
  }
}
//____________________________________________________________________________
void GuiTableRenderer::DrawGraphInCanvas(
        genie::nuvld::DBTable<DBElDiffXSecTableRow> * table)
{
  TGraphAsymmErrors * graph = table->GetGraph(fDrawOpt.c_str(), fPlotVar.c_str());

  TCanvas * c = fEmbCanvas->GetCanvas();

  if(graph != 0) {
    double xmin = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    double xmax = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    double ymin = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    double ymax = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];

    // check if there is a user-defined energy Range
    if(fHaveCustomXRange) {
        xmin = fXmin;
        xmax = fXmax;
    }

    c->Clear();
    c->Divide(2,1);

    c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.95);
    c->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);

    c->GetPad(1)->SetFillColor(0);
    c->GetPad(1)->SetBorderMode(0);
    c->GetPad(2)->SetFillColor(0);
    c->GetPad(2)->SetBorderMode(0);

    c->GetPad(1)->cd();
    c->GetPad(1)->SetBorderMode(0);
    
    TH1F * hframe = (TH1F*) c->GetPad(1)->DrawFrame(.8*xmin, .8*ymin, 1.2*xmax, 1.2*ymax);
    hframe->Draw();

    graph->SetMarkerStyle(7);
    graph->SetMarkerSize(1);

    graph->Draw("P");

    c->GetPad(1)->cd();

    hframe->GetXaxis()->SetTitle(fPlotVar.c_str());
    hframe->GetYaxis()->SetTitle("d#sigma/dEd#Omega (nb/GeV*sr)");

    hframe->GetXaxis()->SetTitleFont(32);
    hframe->GetXaxis()->SetTitleColor(1);
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleFont(32);
    hframe->GetYaxis()->SetTitleColor(1);
    hframe->GetYaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleOffset(1.1);

    // check if it is better to plot it in log scale
    if(xmin > 0 && xmax > xmin) {
        double Range = log10(xmax/xmin);
        if(Range > kRangeInLog)  c->GetPad(1)->SetLogx();
        else                     c->GetPad(1)->SetLogx(0);
    }
    if(ymin > 0 && ymax > ymin) {
        double Range = log10(ymax/ymin);
        if(Range > kRangeInLog)  c->GetPad(1)->SetLogy();
        else                     c->GetPad(1)->SetLogy(0);
    }
  }

  c->GetPad(1)->Update();
  c->Update();
}
//____________________________________________________________________________
void GuiTableRenderer::DrawMultiGraphInCanvas (
         genie::nuvld::DBTable<DBElDiffXSecTableRow> * table)
{
  MultiGraph * mgraph = table->GetMultiGraph(fDrawOpt.c_str(), fPlotVar.c_str());

  TCanvas * c = fEmbCanvas->GetCanvas();
  
  TH1F * hframe = 0;

  double xmin = 0, xmax = 0, ymin = 0, ymax = 0;

  TGraphAsymmErrors * graph = table->GetGraph(fDrawOpt.c_str(), fPlotVar.c_str());

  if(graph != 0) {
    xmin  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    xmax  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    ymin  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    ymax  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];

    delete graph;
  }
  
  if(mgraph != 0) {

   // check if there is a user-defined energy Range
    if(fHaveCustomXRange) {
        xmin = fXmin;
        xmax = fXmax;
    }

    SLOG("NuVld", pDEBUG)
                  << " x = [" << xmin << ", " << xmax << "],"
                  << " y = [" << ymin << ", " << ymax << "] ";

    if(fLegend) {

       // fill in the legend but do not plot it. Return the TLegend to be plotted in
       // an external 'key-canvas'. Since there is no legend, use the full canvas.

       c->Clear();
       c->Divide(2,1);

       c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.99);
       c->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);

       c->GetPad(1)->SetFillColor(0);
       c->GetPad(1)->SetBorderMode(0);
       c->GetPad(2)->SetFillColor(0);
       c->GetPad(2)->SetBorderMode(0);

       c->GetPad(1)->cd();
       c->GetPad(1)->SetBorderMode(0);

       hframe = (TH1F*) c->GetPad(1)->DrawFrame(.98*xmin, .98*ymin, 1.02*xmax, 1.02*ymax);
       hframe->Draw();

       for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {
            mgraph->GetGraph(igraph)->SetMarkerStyle(7);
            mgraph->GetGraph(igraph)->SetMarkerSize(1);
            mgraph->GetGraph(igraph)->Draw("P");
       }

       c->GetPad(1)->Update();

       delete fLegend;
       fLegend = new TLegend(0.01, 0.01, 0.99, 0.99);
       fLegend->SetFillColor(0);
       mgraph->FillLegend("LP", fLegend);

       TCanvas * cext = new TCanvas("cext","Multi-Graph Legend",0,0,400,700);
       cext->cd();
       fLegend->SetTextSize(0.03);
       fLegend->Draw();
       cext->Update();
       
    }  else {

       // no input legend - make one and superimpose it on the MultiGraph plotted on
       // the TRootEmbeddedCanvas. Split the canvas to fit both the plot & the legend.

       c->Clear();
       c->Divide(2,1);

       c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.60,0.99);
       c->GetPad(2)->SetPad("legend_pad","",0.61,0.01,0.99,0.99);

       c->GetPad(1)->cd();

       c->GetPad(1)->SetFillColor(0);
       c->GetPad(1)->SetBorderMode(0);
       c->GetPad(2)->SetFillColor(0);
       c->GetPad(2)->SetBorderMode(0);

       hframe = (TH1F*) c->GetPad(1)->DrawFrame(.98*xmin, .98*ymin, 1.02*xmax, 1.02*ymax);
       hframe->Draw();

       for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {

            mgraph->GetGraph(igraph)->SetMarkerStyle(7);
            mgraph->GetGraph(igraph)->SetMarkerSize(1);
            mgraph->GetGraph(igraph)->Draw("P");
       }

       c->GetPad(1)->Update();

       SLOG("NuVld", pDEBUG) << "Drawing legend";

       c->GetPad(2)->cd();

       if(fLegend) delete fLegend;
       fLegend = new TLegend(0.01, 0.01, 0.99, 0.99);
       fLegend->SetFillColor(0);
       mgraph->FillLegend("LP", fLegend);
       fLegend->SetTextSize(0.03);
       fLegend->Draw();

       c->GetPad(2)->Update();
    }

    c->GetPad(1)->cd();

    hframe->GetXaxis()->SetTitle( fPlotVar.c_str() );
    hframe->GetYaxis()->SetTitle("d#sigma/dEd#Omega (nb/GeV*sr)");
    hframe->GetXaxis()->SetTitleFont(32);
    hframe->GetXaxis()->SetTitleColor(1);
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleFont(32);
    hframe->GetYaxis()->SetTitleColor(1);
    hframe->GetYaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleOffset(1.1);


    // check if it is better to plot it in log scale
    if(xmin > 0 && xmax > xmin) {
       double Range = log10(xmax/xmin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogx();
       else                     c->GetPad(1)->SetLogx(0);
    }

    if(ymin > 0 && ymax > ymin) {
       double Range = log10(ymax/ymin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogy();
       else                     c->GetPad(1)->SetLogy(0);
    }

    c->GetPad(1)->Update();
    c->Update();
  }
}
//____________________________________________________________________________
void GuiTableRenderer::DrawGraphInCanvas(
        genie::nuvld::DBTable<DBSFTableRow> * table)
{
  TGraphAsymmErrors * graph = table->GetGraph(fDrawOpt.c_str(), fPlotVar.c_str());

  TCanvas * c = fEmbCanvas->GetCanvas();

  if(graph != 0) {
    double xmin = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    double xmax = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    double ymin = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    double ymax = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];

    // check if there is a user-defined energy Range
    if(fHaveCustomXRange) {
        xmin = fXmin;
        xmax = fXmax;
    }

    c->Clear();
    c->Divide(2,1);

    c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.95);
    c->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);

    c->GetPad(1)->SetFillColor(0);
    c->GetPad(1)->SetBorderMode(0);
    c->GetPad(2)->SetFillColor(0);
    c->GetPad(2)->SetBorderMode(0);

    c->GetPad(1)->cd();
    c->GetPad(1)->SetBorderMode(0);

    TH1F * hframe = (TH1F*) c->GetPad(1)->DrawFrame(.8*xmin, .8*ymin, 1.2*xmax, 1.2*ymax);
    hframe->Draw();

    graph->SetMarkerStyle(7);
    graph->SetMarkerSize(1);

    graph->Draw("P");

    c->GetPad(1)->cd();

    hframe->GetXaxis()->SetTitle(fPlotVar.c_str());
    hframe->GetYaxis()->SetTitle("Structure Function");

    hframe->GetXaxis()->SetTitleFont(32);
    hframe->GetXaxis()->SetTitleColor(1);
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleFont(32);
    hframe->GetYaxis()->SetTitleColor(1);
    hframe->GetYaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleOffset(1.1);

    // check if it is better to plot it in log scale
    if(xmin > 0 && xmax > xmin) {
        double Range = log10(xmax/xmin);
        if(Range > kRangeInLog)  c->GetPad(1)->SetLogx();
        else                     c->GetPad(1)->SetLogx(0);
    }
    if(ymin > 0 && ymax > ymin) {
        double Range = log10(ymax/ymin);
        if(Range > kRangeInLog)  c->GetPad(1)->SetLogy();
        else                     c->GetPad(1)->SetLogy(0);
    }
  }

  c->GetPad(1)->Update();
  c->Update();
}
//____________________________________________________________________________
void GuiTableRenderer::DrawMultiGraphInCanvas (
       genie::nuvld::DBTable<DBSFTableRow> * table)
{
  MultiGraph * mgraph = table->GetMultiGraph(fDrawOpt.c_str(), fPlotVar.c_str());

  TCanvas * c = fEmbCanvas->GetCanvas();

  TH1F * hframe = 0;

  double xmin = 0, xmax = 0, ymin = 0, ymax = 0;

  TGraphAsymmErrors * graph = table->GetGraph(fDrawOpt.c_str(), fPlotVar.c_str());

  if(graph != 0) {
    xmin  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    xmax  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    ymin  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    ymax  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];

    delete graph;
  }
  
  if(mgraph != 0) {
    
   // check if there is a user-defined energy Range
    if(fHaveCustomXRange) {
        xmin = fXmin;
        xmax = fXmax;
    }

    SLOG("NuVld", pDEBUG)
                  << " x = [" << xmin << ", " << xmax << "],"
                  << " y = [" << ymin << ", " << ymax << "] ";

    if(fLegend) {

       // fill in the legend but do not plot it. Return the TLegend to be plotted in
       // an external 'key-canvas'. Since there is no legend, use the full canvas.

       c->Clear();
       c->Divide(2,1);

       c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.98,0.99);
       c->GetPad(2)->SetPad("legend_pad","",0.98,0.01,0.99,0.99);

       c->GetPad(1)->SetFillColor(0);
       c->GetPad(1)->SetBorderMode(0);
       c->GetPad(2)->SetFillColor(0);
       c->GetPad(2)->SetBorderMode(0);

       c->GetPad(1)->cd();
       c->GetPad(1)->SetBorderMode(0);

       hframe = (TH1F*) c->GetPad(1)->DrawFrame(.98*xmin, .98*ymin, 1.02*xmax, 1.02*ymax);
       hframe->Draw();

       for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {
            mgraph->GetGraph(igraph)->SetMarkerStyle(7);
            mgraph->GetGraph(igraph)->SetMarkerSize(1);
            mgraph->GetGraph(igraph)->Draw("P");
       }

       c->GetPad(1)->Update();

       delete fLegend;
       fLegend = new TLegend(0.01, 0.01, 0.99, 0.99);
       fLegend->SetFillColor(0);
       mgraph->FillLegend("LP", fLegend);

       TCanvas * cext = new TCanvas("cext","Multi-Graph Legend",0,0,400,700);
       cext->cd();
       fLegend->SetTextSize(0.03);
       fLegend->Draw();
       cext->Update();

    }  else {

       // no input legend - make one and superimpose it on the MultiGraph plotted on
       // the TRootEmbeddedCanvas. Split the canvas to fit both the plot & the legend.

       c->Clear();
       c->Divide(2,1);

       c->GetPad(1)->SetPad("mplots_pad","",0.01,0.01,0.60,0.99);
       c->GetPad(2)->SetPad("legend_pad","",0.61,0.01,0.99,0.99);

       c->GetPad(1)->cd();

       c->GetPad(1)->SetFillColor(0);
       c->GetPad(1)->SetBorderMode(0);
       c->GetPad(2)->SetFillColor(0);
       c->GetPad(2)->SetBorderMode(0);

       hframe = (TH1F*) c->GetPad(1)->DrawFrame(.98*xmin, .98*ymin, 1.02*xmax, 1.02*ymax);
       hframe->Draw();

       for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {

            mgraph->GetGraph(igraph)->SetMarkerStyle(7);
            mgraph->GetGraph(igraph)->SetMarkerSize(1);
            mgraph->GetGraph(igraph)->Draw("P");
       }

       c->GetPad(1)->Update();

       SLOG("NuVld", pDEBUG) << "Drawing legend";

       c->GetPad(2)->cd();

       if(fLegend) delete fLegend;
       fLegend = new TLegend(0.01, 0.01, 0.99, 0.99);
       fLegend->SetFillColor(0);
       mgraph->FillLegend("LP", fLegend);
       fLegend->SetTextSize(0.03);
       fLegend->Draw();

       c->GetPad(2)->Update();
    }

    c->GetPad(1)->cd();

    hframe->GetXaxis()->SetTitle( fPlotVar.c_str() );
    hframe->GetYaxis()->SetTitle("Structure Function");
    hframe->GetXaxis()->SetTitleFont(32);
    hframe->GetXaxis()->SetTitleColor(1);
    hframe->GetXaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleFont(32);
    hframe->GetYaxis()->SetTitleColor(1);
    hframe->GetYaxis()->SetTitleSize(0.05);
    hframe->GetYaxis()->SetTitleOffset(1.1);


    // check if it is better to plot it in log scale
    if(xmin > 0 && xmax > xmin) {
       double Range = log10(xmax/xmin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogx();
       else                     c->GetPad(1)->SetLogx(0);
    }

    if(ymin > 0 && ymax > ymin) {
       double Range = log10(ymax/ymin);
       if(Range > kRangeInLog)  c->GetPad(1)->SetLogy();
       else                     c->GetPad(1)->SetLogy(0);
    }

    c->GetPad(1)->Update();
    c->Update();
  }
}
//____________________________________________________________________________
