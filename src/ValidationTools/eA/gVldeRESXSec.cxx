//____________________________________________________________________________
/*!

\program gvld_e_res_xsec

\brief   Compares GENIE with electron scattering data in the resonance region

         Syntax:
           gvld_e_res_xsec 
                [-h host] [-u user] [-p passwd] [-m genie_model]

         Options:

           [] Denotes an optional argument.

           -h NuVld MySQL URL (eg mysql://localhost/NuScat).
           -u NuVld MySQL username.
           -p NuVld MySQL password.
           -m Specifies a GENIE model 

         Example:

           .. .. .. ..

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created Oct 16, 2009 

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TChain.h>

#include "Conventions/GBuild.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
#include "Utils/StringUtils.h"
#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBStatus.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::nuvld;

/* 
..............................................................................
ELECTRON SCATTERING CROSS SECTION DATA IN THE RESONANCE REGION
..............................................................................
*/
const int    kElXSecDataSets = 4;
const char * kElXSecDataSetLabel[kElXSecDataSets] = 
{
/* 0 */ "JLAB, E = 2.445 GeV, #theta = 20.0 deg",
/* 1 */ "JLAB, E = 2.445 GeV, #theta = 30.0 deg",
/* 2 */ "JLAB, E = 2.445 GeV, #theta = 35.5 deg",
/* 3 */ "JLAB, E = 2.445 GeV, #theta = 70.0 deg"
//
// ... add more ...
//
};
const char * kElXSecKeyList[kElXSecDataSets] = {
/* 0 */ "JLAB,0",
/* 1 */ "JLAB,0",
/* 2 */ "JLAB,0",
/* 3 */ "JLAB,0"
//
// ... add more ...
//
};
float kElXSecEnergy[kElXSecDataSets] = {
/* 0 */ 2.445,
/* 1 */ 2.445,
/* 2 */ 2.445,
/* 3 */ 2.445
//
// ... add more ...
//
};
float kElXSecTheta[kElXSecDataSets] = {
/* 0 */ 20.00,
/* 1 */ 30.00,
/* 2 */ 38.50,
/* 3 */ 70.01
//
// ... add more ...
//
};

typedef DBQueryString                 DBQ;
typedef DBTable<DBElDiffXSecTableRow> DBT;

// function prototypes
void     Init               (void);
void     Plot               (void);
void     End                (void);
void     SetStyle           (void);
void     AddCoverPage       (void);
bool     Connect            (void);
DBQ      FormQuery          (const char * key_list, float energy, float theta);
DBT *    Data               (int iset);
TGraph * Model              (int iset, int imodel);
void     Draw               (int iset);
TH1F *   DrawFrame          (TGraph * gr0, TGraph * gr1);
TH1F *   DrawFrame          (double xmin, double xmax, double ymin, double yman);
void     Format             (TGraph* gr, int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
void     GetCommandLineArgs (int argc, char ** argv);
void     PrintSyntax        (void);

// command-line arguments
string         gOptDbURL;
string         gOptDbUser;
string         gOptDbPasswd;

// dbase information
const char * kDefDbURL = "mysql://localhost/NuScat";  

// globals
bool            gCmpWithData  = true;
DBI *           gDBI          = 0;
TPostScript *   gPS           = 0;
TCanvas *       gC            = 0;

const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);

  Init();
  Plot();
  End();

  LOG("gvldtest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void Plot(void)
{
#ifdef __GENIE_MYSQL_ENABLED__

  // connect to the NuValidator MySQL dbase
  bool ok = Connect();
  if(!ok) {
    return;
  }
 
  // loop over data sets
  for(int iset = 0; iset < kElXSecDataSets; iset++) 
  {
    Draw(iset);
  }
#endif
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("vldtest", pNOTICE) << "Initializing...";

  // genie style
  SetStyle();

  // canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // output postscript file
  gPS = new TPostScript("genie_eres_vs_data.ps", 111);

  // cover page
  AddCoverPage();
}
//_________________________________________________________________________________
void AddCoverPage(void)
{
  // header
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE Resonance Electro-Production vs Data");
  hdr.AddText(" ");
  hdr.Draw();
  gC->Update();
}
//_________________________________________________________________________________
void End(void)
{
  LOG("vldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;
}
//_________________________________________________________________________________
// Corresponding GENIE prediction for the `iset' data set 
//.................................................................................
TGraph * Model(int iset, int imodel)
{
  LOG("vldtest", pNOTICE) 
    << "Getting GENIE prediction (model ID = " 
    << imodel << ", data set ID = " << iset << ")";

  switch(iset) {

    case (0) :
    case (1) :
    case (2) :
    case (3) :
    {
       return 0;
       break;
    }

    default:
    {
       return 0;
       break;
    }
  }
  return 0;
}
//_________________________________________________________________________________
// Download cross section data from NuVld MySQL dbase 
//.................................................................................
bool Connect(void)
{
  if(!gCmpWithData) return true;

  // Get a data-base interface
  TSQLServer * sql_server = TSQLServer::Connect(
      gOptDbURL.c_str(),gOptDbUser.c_str(),gOptDbPasswd.c_str());

  if(!sql_server) return false;
  if(!sql_server->IsConnected()) return false;

  gDBI = new DBI(sql_server);
  return true;
}
//_________________________________________________________________________________
DBQ FormQuery(const char * key_list, float energy, float theta)
{
// forms a DBQueryString for extracting neutrino cross section data from the input 
// key-list and for the input energy range
//  
  ostringstream query_string;
  
  query_string 
    << "KEY-LIST:" << key_list
    << "$CUTS:E_min=" << energy-0.001 << ";E_max=" << energy+0.001 
    << ";Theta_min=" << theta-0.001 << ";Theta_max=" << theta+0.001 
    << "$DRAW_OPT:none$DB-TYPE:eN-Diff-XSec";
  
  DBQ query(query_string.str());
  
  return query;
}
//_________________________________________________________________________________
DBT * Data(int iset)
{
  if(!gCmpWithData) return 0;

  DBT * dbtable = new DBT;

  const char * keylist = kElXSecKeyList[iset];
  float        energy  = kElXSecEnergy [iset];
  float        theta   = kElXSecTheta  [iset];

  DBQ query = FormQuery(keylist, energy, theta);
  assert( gDBI->FillTable(dbtable, query) == eDbu_OK );

  return dbtable;
}
//_________________________________________________________________________________
void Draw(int iset)
{
  // get all measurements for the current channel from the NuValidator MySQL dbase
  DBT * dbtable = Data(iset);

  // get the corresponding GENIE model prediction
  //
  // ...
  //

  if( /* not genie data && */ !dbtable) return;

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

//  TLatex * title = new TLatex(1,600,kElXSecDataSetLabel[iset]);
//  title->SetNDC(false);
//  title->Draw();

  TH1F * hframe = 0;
  bool have_frame = false;

  // have data points to plot?
  if(dbtable) {
    TGraphAsymmErrors * graph = dbtable->GetGraph("err","W2");

    // create frame from the data point range
    double xmin  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    double xmax  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    double ymin  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    double ymax  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];
    hframe = (TH1F*) gC->GetPad(1)->DrawFrame(0.5*xmin, 0.4*ymin, 1.2*xmax, 1.2*ymax);
    hframe->SetTitle(kElXSecDataSetLabel[iset]);
//    hframe->Draw();
    have_frame = true;

    //
    // draw current data set
    //
    MultiGraph * mgraph = dbtable->GetMultiGraph("err","W2");
    for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {
       mgraph->GetGraph(igraph)->Draw("P");
    }
  }//dbtable?

  // have model prediction to plot?
/*
  if( have_genie_prediction? ) {
     if(!have_frame) {
        // the data points have not been plotted
        // create a frame from this graph range
        double xmin  = ...
        double xmax  = ...
        double ymin  = ... // get from output graph
        double ymax  = ...
        hframe = (TH1F*) gC->GetPad(1)->DrawFrame(0.5*xmin, 0.4*ymin, 1.2*xmax, 1.2*ymax);
        hframe->Draw();
     }
     for(loop over models, if multiple one are used) {
       TGraph * plot = get plot for current model;
       if(plot) {
         int lsty = kLStyle[imodel];     
         Format(plot,1,lsty,2,1,1,1);
         plot->Draw("L");
       }
     }
  }//model?
*/

  hframe->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
  hframe->GetYaxis()->SetTitle("d^{2}#sigma / d#Omega dE (nb/sr/GeV)");

  gC->GetPad(iplot)->Update();
  gC->Update();

  if(dbtable) {
    delete dbtable;
  }
}
//_________________________________________________________________________________
// Formatting
//.................................................................................
void SetStyle(void)
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
}
//_________________________________________________________________________________
void Format(
    TGraph* gr, int lcol, int lsty, int lwid, int mcol, int msty, double msiz)
{
  if(!gr) return;

  if (lcol >= 0) gr -> SetLineColor   (lcol);
  if (lsty >= 0) gr -> SetLineStyle   (lsty);
  if (lwid >= 0) gr -> SetLineWidth   (lwid);

  if (mcol >= 0) gr -> SetMarkerColor (mcol);
  if (msty >= 0) gr -> SetMarkerStyle (msty);
  if (msiz >= 0) gr -> SetMarkerSize  (msiz);
}
//_________________________________________________________________________________
// Parsing command-line arguments, check/form filenames, etc
//.................................................................................
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  gCmpWithData = true;

  // get DB URL
  try {
     gOptDbURL = utils::clap::CmdLineArgAsString(argc,argv,'h');
  } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
       gOptDbURL = kDefDbURL;
     }
  }

  // get DB username
  try {
     gOptDbUser = utils::clap::CmdLineArgAsString(argc,argv,'u');
  } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
       gCmpWithData = false;
     }
  }

  // get DB passwd
  try {
     gOptDbPasswd = utils::clap::CmdLineArgAsString(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
       gCmpWithData = false;
     }
  }

}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gvld_nuxsec_vs_world_data [-h host] [-u user] [-p passwd] -f files\n";
}
//_________________________________________________________________________________

