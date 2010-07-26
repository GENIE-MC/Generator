//____________________________________________________________________________
/*!

\program gvld_sf_vs_world_data

\brief   A GENIE utility that generates structure function comparison plots

         Syntax:
           gvld_sf_vs_world_data [-h host] [-u user] [-p passwd] [-g inputs]

         Options:
           [] Denotes an optional argument.
           -h NuVld MySQL URL (eg mysql://localhost/NuScat)
           -u NuVld MySQL username 
           -p NuVld MySQL password
           -g GENIE inputs
		      
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 06, 2008 

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TPostScript.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TStyle.h>
#include <TLegend.h>

#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "ValidationTools/NuVld/DBI.h"
#include "ValidationTools/NuVld/DBStatus.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::nuvld;

/*
..............................................................................
STRUCTURE FUNCTION DATA
..............................................................................
ID   DESCRIPTION
0    neutrino F2 ....... [all, x = 0.007]
1    neutrino F2 ....... [all, x = 0.013]
2    neutrino F2 ....... [all, x = 0.015]
3    neutrino F2 ....... [all, x = 0.018]
4    neutrino F2 ....... [all, x = 0.025]
5    neutrino F2 ....... [all, x = 0.035]
6    neutrino F2 ....... [all, x = 0.045]
7    neutrino F2 ....... [all, x = 0.050]
8    neutrino F2 ....... [all, x = 0.070]
9    neutrino F2 ....... [all, x = 0.080]
10   neutrino F2 ....... [all, x = 0.090]
11   neutrino F2 ....... [all, x = 0.110]
12   neutrino F2 ....... [all, x = 0.125]
13   neutrino F2 ....... [all, x = 0.140]
14   neutrino F2 ....... [all, x = 0.175]
15   neutrino F2 ....... [all, x = 0.180]
16   neutrino F2 ....... [all, x = 0.225]
17   neutrino F2 ....... [all, x = 0.275]
18   neutrino F2 ....... [all, x = 0.350]
19   neutrino F2 ....... [all, x = 0.450]
20   neutrino F2 ....... [all, x = 0.550]
21   neutrino F2 ....... [all, x = 0.650]
22   neutrino F2 ....... [all, x = 0.750]
23   charged lepton F2.. [all, x = 0.049-0.051]
24   xF3................ [all, x = 0.049-0.051]
..............................................................................
*/
const int   kSFDataSets = 25;
const float kSFQ2min  = 1;
const float kSFQ2max  = 400;
const int   kSFRawDis = 2;

const char * kSFDataSetLabel[kSFDataSets] = {
/* 0  */ "neutrino F2  (all/Fe56, x = 0.007)",
/* 1  */ "neutrino F2  (all/Fe56, x = 0.013)",
/* 2  */ "neutrino F2  (all/Fe56, x = 0.015)",
/* 3  */ "neutrino F2  (all/Fe56, x = 0.018)",
/* 4  */ "neutrino F2  (all/Fe56, x = 0.025)",
/* 5  */ "neutrino F2  (all/Fe56, x = 0.035)",
/* 6  */ "neutrino F2  (all/Fe56, x = 0.045)",
/* 7  */ "neutrino F2  (all/Fe56, x = 0.050)",
/* 8  */ "neutrino F2  (all/Fe56, x = 0.070)",
/* 9  */ "neutrino F2  (all/Fe56, x = 0.080)",
/* 10 */ "neutrino F2  (all/Fe56, x = 0.090)",
/* 11 */ "neutrino F2  (all/Fe56, x = 0.110)",
/* 12 */ "neutrino F2  (all/Fe56, x = 0.125)",
/* 13 */ "neutrino F2  (all/Fe56, x = 0.140)",
/* 14 */ "neutrino F2  (all/Fe56, x = 0.175)",
/* 15 */ "neutrino F2  (all/Fe56, x = 0.180)",
/* 16 */ "neutrino F2  (all/Fe56, x = 0.225)",
/* 17 */ "neutrino F2  (all/Fe56, x = 0.275)",
/* 18 */ "neutrino F2  (all/Fe56, x = 0.350)",
/* 19 */ "neutrino F2  (all/Fe56, x = 0.450)",
/* 20 */ "neutrino F2  (all/Fe56, x = 0.550)",
/* 21 */ "neutrino F2  (all/Fe56, x = 0.650)",
/* 22 */ "neutrino F2  (all/Fe56, x = 0.750)",
/* 23 */ "charged lepton F2 (all, x = 0.049-0.051) ",
/* 24 */ "xF3 (all, x = 0.049-0.051) "
};

const char * kSFKeyList[kSFDataSets] = {
/* 0  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 1  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 2  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 3  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 4  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 5  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 6  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 7  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 8  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 9  */ "CCFR,0;CDHS,0;NUTEV,0",
/* 10 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 11 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 12 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 13 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 14 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 15 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 16 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 17 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 18 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 19 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 20 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 21 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 22 */ "CCFR,0;CDHS,0;NUTEV,0",
/* 23 */ "BCDMS,0;EMC,0;NMC,0;SLAC,0",
/* 24 */ "CCFR,1;CDHS,1;NUTEV,1;WA59,1"
};
float kSFx[kSFDataSets] = {
/* 0  */ 0.007,
/* 1  */ 0.013,
/* 2  */ 0.015,
/* 3  */ 0.018,
/* 4  */ 0.025,
/* 5  */ 0.035,
/* 6  */ 0.045,
/* 7  */ 0.050,
/* 8  */ 0.070,
/* 9  */ 0.080,
/* 10 */ 0.090,
/* 11 */ 0.110,
/* 12 */ 0.125,
/* 13 */ 0.140,
/* 14 */ 0.175,
/* 15 */ 0.180,
/* 16 */ 0.225,
/* 17 */ 0.275,
/* 18 */ 0.350,
/* 19 */ 0.450,
/* 20 */ 0.550,
/* 21 */ 0.650,
/* 22 */ 0.750,
/* 23 */ 0.050,
/* 24 */ 0.050
};
int kSFA[kSFDataSets] = {
/* 0  */ 56,
/* 1  */ 56,
/* 2  */ 56,
/* 3  */ 56,
/* 4  */ 56,
/* 5  */ 56,
/* 6  */ 56,
/* 7  */ 56,
/* 8  */ 56,
/* 9  */ 56,
/* 10 */ 56,
/* 11 */ 56,
/* 12 */ 56,
/* 13 */ 56,
/* 14 */ 56,
/* 15 */ 56,
/* 16 */ 56,
/* 17 */ 56,
/* 18 */ 56,
/* 19 */ 56,
/* 20 */ 56,
/* 21 */ 56,
/* 22 */ 56,
/* 23 */ 56,
/* 24 */ 56
};
float kSFScale[kSFDataSets] = {
/* 0  */ 1.0,
/* 1  */ 1.0,
/* 2  */ 1.0,
/* 3  */ 1.0,
/* 4  */ 1.0,
/* 5  */ 1.0,
/* 6  */ 1.0,
/* 7  */ 1.0,
/* 8  */ 1.0,
/* 9  */ 1.0,
/* 10 */ 1.0,
/* 11 */ 1.0,
/* 12 */ 1.0,
/* 13 */ 1.0,
/* 14 */ 1.0,
/* 15 */ 1.0,
/* 16 */ 1.0,
/* 17 */ 1.0,
/* 18 */ 1.0,
/* 19 */ 1.0,
/* 20 */ 1.0,
/* 21 */ 1.0,
/* 22 */ 1.0,
/* 23 */ 1.0,
/* 24 */ 1.0
};

typedef DBQueryString         DBQ;
typedef DBTable<DBSFTableRow> DBT;

// function prototypes
void     Init               (void);
void     Plot               (void);
void     End                (void);
void     AddCoverPage       (void);
bool     Connect            (void);
DBQ      FormQuery          (const char * key_list, float Q2min, float Q2max, float xmin, float xmax);
DBT *    Data               (int iset);
void     Draw               (int iset);
TH1F *   DrawFrame          (TGraph * gr0, TGraph * gr1, TCanvas * c);
void     Format             (TGraph* gr, int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
void     Draw               (TGraph* gr, const char * opt);
void     GetCommandLineArgs (int argc, char ** argv);
void     PrintSyntax        (void);

// command-line arguments
string gOptDbURL;
string gOptDbUser;
string gOptDbPasswd;

// default command line arguments
const char * kDefDbURL = "mysql://localhost/NuScat";  

// globals
bool            gCmpWithData  = true;
DBI *           gDBI          = 0;
TPostScript *   gPS           = 0;
TCanvas *       gC            = 0;
TLegend *       gLS           = 0;

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
  for(int iset = 0; iset < kSFDataSets; iset++) 
  {
    Draw(iset);
  }
#endif
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("vldtest", pNOTICE) << "Initializing...";;

  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  gLS = new TLegend(0.15,0.92,0.85,0.98);
  gLS -> SetFillColor(0);
  gLS -> SetBorderSize(1);

  // output file
  gPS = new TPostScript("genie_sf_vs_data.ps", 111);

  AddCoverPage();

  gC->SetLogx();
  gC->SetLogy();
}
//_________________________________________________________________________________
void AddCoverPage(void)
{
  // header
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE Structure Function Comparisons with World Data");
  hdr.AddText(" ");
/*
  hdr.AddText(" ");
  for(int imodel=0; imodel< gOptGenieInputs.NModels(); imodel++) {
    ostringstream stream;
    stream << "model tag: " << gOptGenieInputs.ModelTag(imodel)
           << "(" << kLStyleTxt[imodel] << ")";
    hdr.AddText(stream.str().c_str());
  }
  hdr.AddText(" ");
*/
  hdr.Draw();
  gC->Update();
}
//_________________________________________________________________________________
void End(void)
{
  LOG("vldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gLS;
  delete gPS;
}
//_________________________________________________________________________________
// Corresponding GENIE prediction for the `iset' data set 
//.................................................................................
TGraph * Model(int /*iset*/, int /*imodel*/)
{
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
DBQ FormQuery(
    const char * key_list, float Q2min, float Q2max, float xmin, float xmax)
{
// forms a DBQueryString for extracting structure function data from the input 
// key-list and for the Q2 range

  ostringstream query_string;

  query_string << "KEY-LIST:" << key_list;
  query_string << "$CUTS:Q2min=" << Q2min << ";Q2max=" << Q2max
               << ";xmin=" << xmin << ";xmax=" << xmax
               << "$DRAW_OPT:none$DB-TYPE:SF";

  DBQ query(query_string.str());

  return query;
}
//_________________________________________________________________________________
DBT * Data(int iset)
{
  DBT * dbtable = new DBT;

  const char * keylist = kSFKeyList[iset];
  float        xmin    = kSFx[iset] - 0.001;
  float        xmax    = kSFx[iset] + 0.001;

  DBQ query = FormQuery(keylist, kSFQ2min, kSFQ2max, xmin, xmax);
  
  assert( gDBI->FillTable(dbtable, query) == eDbu_OK );

  return dbtable;
}
//_________________________________________________________________________________
void Draw(int iset)
{
  // get all measurements for the current channel from the NuValidator MySQL dbase
  DBT * dbtable = Data(iset);

  // get the corresponding GENIE model prediction
  vector<TGraph*> models;
/*
  for(int imodel=0; imodel< gOptGenieInputs.NModels(); imodel++) {
       models.push_back(Model(iset,imodel));
  }
*/

  if(models.size()==0 && !dbtable) return;

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

  gLS->SetHeader(kSFDataSetLabel[iset]);

  TLegend * legend = new TLegend(0.01, 0.01, 0.99, 0.99);
  legend->SetFillColor(0);
  legend->SetTextSize(0.08);

  TH1F * hframe = 0;
  bool have_frame = false;

 // have data points to plot?
  if(dbtable) {
    TGraphAsymmErrors * graph = dbtable->GetGraph("err","Q2");

    // create frame from the data point range
    double xmin  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    double xmax  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    double ymin  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    double ymax  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];
    hframe = (TH1F*) gC->GetPad(1)->DrawFrame(0.5*xmin, 0.4*ymin, 1.2*xmax, 2.0*ymax);
    hframe->Draw();
    have_frame = true;

    //
    // draw current data set
    //
    MultiGraph * mgraph = dbtable->GetMultiGraph("err","Q2");
    for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {
       mgraph->GetGraph(igraph)->Draw("P");
    }
    mgraph->FillLegend("LP", legend);
  }//dbtable?

  // have model prediction to plot?
/*
  if(models.size()>0) {
     if(!have_frame) {
        // the data points have not been plotted
        // create a frame from this graph range
        double xmin  = ( models[0]->GetX() )[TMath::LocMin(models[0]->GetN(),models[0]->GetX())];
        double xmax  = ( models[0]->GetX() )[TMath::LocMax(models[0]->GetN(),models[0]->GetX())];
        double ymin  = ( models[0]->GetY() )[TMath::LocMin(models[0]->GetN(),models[0]->GetY())];
        double ymax  = ( models[0]->GetY() )[TMath::LocMax(models[0]->GetN(),models[0]->GetY())];
        hframe = (TH1F*) gC->GetPad(1)->DrawFrame(0.5*xmin, 0.4*ymin, 1.2*xmax, 2.0*ymax);
        hframe->Draw();
     }
     for(int imodel=0; imodel<gOptGenieInputs.NModels(); imodel++) {
       TGraph * plot = models[imodel];
       if(plot) {
         int lsty = kLStyle[imodel];     
         Format(plot,1,lsty,2,1,1,1);
         plot->Draw("L");
       }
     }
  }//model?
*/

/*
  hframe->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hframe->GetYaxis()->SetTitle("#sigma_{#nu} (1E-38 cm^{2})");
*/
  gLS->Draw();

  gC->GetPad(2)->cd();
  legend->Draw();

  gC->GetPad(2)->Update();
  gC->Update();

  if(dbtable) {
    delete dbtable;
  }
}
//_________________________________________________________________________________
// Formatting
//.................................................................................
TH1F* DrawFrame(TGraph * gr0, TGraph * gr1, TCanvas * c)
{
  TAxis * x0 = gr0 -> GetXaxis();
  TAxis * y0 = gr0 -> GetYaxis();
  double xmin = x0 -> GetXmin();
  double xmax = x0 -> GetXmax();
  double ymin = y0 -> GetXmin();
  double ymax = y0 -> GetXmax();
  if(gr1) {
     TAxis * x1 = gr1 -> GetXaxis();
     TAxis * y1 = gr1 -> GetYaxis();
     xmin = TMath::Min(xmin, x1 -> GetXmin());
     xmax = TMath::Max(xmax, x1 -> GetXmax());
     ymin = TMath::Min(ymin, y1 -> GetXmin());
     ymax = TMath::Max(ymax, y1 -> GetXmax());
  }
  xmin *= 0.5;
  xmax *= 1.5;
  ymin *= 0.5;
  ymax *= 1.5;
  xmin = TMath::Max(0.1, xmin);
  
  TH1F * hf = (TH1F*) c->DrawFrame(xmin, ymin, xmax, ymax);
  hf->GetXaxis()->SetTitle("E (GeV)");
  hf->GetYaxis()->SetTitle("#sigma (10^{-38} cm^{2})");
  hf->GetYaxis()->SetTitleSize(0.03);
  hf->GetYaxis()->SetTitleOffset(1.3);
  hf->GetXaxis()->SetLabelSize(0.03);
  hf->GetYaxis()->SetLabelSize(0.03);
  return hf;
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
void Draw(TGraph* gr, const char * opt)
{
  if(!gr) return;
  gr->Draw(opt);
}
//_________________________________________________________________________________
// Parsing command-line arguments, check/form filenames, etc
//.................................................................................
void GetCommandLineArgs(int argc, char ** argv)
{
  CmdLnArgParser parser(argc,argv);

  gCmpWithData = true;

  // get DB URL
  if( parser.OptionExists('h') ) {
     gOptDbURL = parser.ArgAsString('h');
  } else {
     gOptDbURL = kDefDbURL;
  }

  // get DB username
  if( parser.OptionExists('u') ) {
     gOptDbUser = parser.ArgAsString('u');
  } else {
     gCmpWithData = false;
  }

  // get DB passwd
  if( parser.OptionExists('p') ) {
     gOptDbPasswd = parser.ArgAsString('p');
  } else {
     gCmpWithData = false;
  }
}
//_________________________________________________________________________________
void PrintSyntax(void)
{

}
//_________________________________________________________________________________
