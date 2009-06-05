//____________________________________________________________________________
/*!

\program gvld_structfunc_comp

\brief   A GENIE utility that generates structure function comparison plots

         Syntax:
           gvld_structfunc_comp -f sf_file [-r reference_sf_file] 
                         [-d] [-h host] [-u user] [-p passwd]

         Options:
           [] Denotes an optional argument.
           -f Specifies a ROOT file with a GENIE SF ntuple.
	   -r Specifies a reference ROOT file with a reference GENIE SF ntuple.
           -d Enable comparisons with data.
           -h NuVld MySQL URL (eg mysql://localhost/NuScat)
           -u NuVld MySQL username 
           -p NuVld MySQL password
		      
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 06, 2008 

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
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
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
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
/* 24 */ "xF3               (all, x = 0.049-0.051) "
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

// function prototypes
void                    GetCommandLineArgs (int argc, char ** argv);
void                    PrintSyntax        (void);
bool                    CheckRootFilename  (string filename);
string                  OutputFileName     (string input_file_name);
void                    CompareWithData    (void);
TH1F *                  DrawFrame          (TGraph * gr0, TGraph * gr1, TCanvas * c);
void                    Format             (TGraph* gr, int lcol, int lsty, int lwid, int mcol, int msty, double msiz);
void                    Draw               (TGraph* gr, const char * opt);
DBQueryString           FormQuery          (const char * key_list, float Q2min, float Q2max, float xmin, float xmax);
DBTable<DBSFTableRow> * GetNuVldData       (int iset);


// command-line arguments
string gOptSfFilename_curr = "";  // (-f) input ROOT cross section file
string gOptSfFilename_ref0 = "";  // (-r) input ROOT cross section file (reference)
bool   gOptHaveRef;
bool   gOptCmpWithData;
string gOptDbURL;
string gOptDbUser;
string gOptDbPasswd;

// default command line arguments
const char * kDefDbURL = "mysql://localhost/NuScat";  

// globals
TFile *  gSfFile_curr = 0;
TFile *  gSfFile_ref0 = 0;
DBI *    dbi          = 0;

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);
  CompareWithData();

  LOG("gvldtest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void CompareWithData(void)
{
  // Get a data-base interface
  TSQLServer * sql_server = TSQLServer::Connect(
      gOptDbURL.c_str(),gOptDbUser.c_str(),gOptDbPasswd.c_str());
  assert(sql_server && sql_server->IsConnected());
  dbi = new DBI(sql_server);
  assert(dbi);
 
  TCanvas * c = new TCanvas("c","",20,20,500,650);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  c->SetGridx();
  c->SetGridy();
 
  TLegend * l  = new TLegend(0.80,0.20,0.99,0.99);
  TLegend * ls = new TLegend(0.15,0.85,0.55,0.95);
  l ->SetFillColor(0);
  ls->SetFillColor(0);
  ls->SetBorderSize(1);
  TPostScript * ps = new TPostScript("gsf.ps", 111);
   
  //
  // header
  //
  ps->NewPage();
  c->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText("GENIE structure function comparison plots");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText("Notes:");
  hdr.AddText("");
  hdr.AddText("");
  hdr.Draw();
  c->Update();
  TPavesText title(10,40,90,70,1,"ndc");
  
  double xmin = 9999;
  double xmax = -1;
  double ymin = 9999;
  double ymax = -1;

  //
  // summarize all comparisons in a single page
  //

  for(int iset = 0; iset < kSFDataSets; iset++) {
    DBTable<DBSFTableRow> * dbtable = GetNuVldData(iset);
    TGraphAsymmErrors * graph = dbtable->GetGraph("err","Q2");

    double xmin_c  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    double xmax_c  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    double ymin_c  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    double ymax_c  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];
  
    xmin = TMath::Min(xmin, xmin_c);
    xmax = TMath::Max(xmax, xmax_c);
    ymin = TMath::Min(ymin, ymin_c);
    ymax = TMath::Max(ymax, ymax_c);
  }

//  c->Clear();
//  c->Divide(2,1);
//  c->GetPad(1)->SetPad("mplots_pad","",0.01,0.25,0.99,0.99);
//  c->GetPad(2)->SetPad("legend_pad","",0.01,0.01,0.99,0.24);
//  c->GetPad(1)->SetFillColor(0);
//  c->GetPad(1)->SetBorderMode(0);
//  c->GetPad(2)->SetFillColor(0);
//  c->GetPad(2)->SetBorderMode(0);
//  c->GetPad(1)->cd();
//  c->GetPad(1)->SetBorderMode(0);
  
//  c->GetPad(1)->SetLogx();
//  c->GetPad(1)->SetLogy();
  c->SetLogx();


//  TH1F * hframe = (TH1F*) c->GetPad(1)->DrawFrame(.5*xmin, .4*ymin, 1.2*xmax, 1.42*ymax);
  TH1F * hframe = (TH1F*) c->DrawFrame(.5*xmin, .4*ymin, 1.2*xmax, 1.42*ymax);
  hframe->Draw();

  for(int iset = 0; iset < kSFDataSets; iset++) {
    DBTable<DBSFTableRow> * dbtable = GetNuVldData(iset);
    MultiGraph * mgraph = dbtable->GetMultiGraph("err","Q2");
    for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {
       mgraph->GetGraph(igraph)->Draw("P");
    }

    // Superimpose GENIE predictions
    //
    // ...

//    c->GetPad(2)->cd();
//    TLegend * legend = new TLegend(0.01, 0.01, 0.99, 0.99);
//    legend->SetFillColor(0);
//    mgraph->FillLegend("LP", legend);
//    legend->SetTextSize(0.08);
//    legend->Draw();    
  }

//  c->GetPad(2)->Update();
  c->Update();      

  //
  // show more detail: one page per SF set
  //

  for(int iset = 0; iset < kSFDataSets; iset++) {
    DBTable<DBSFTableRow> * dbtable = GetNuVldData(iset);
    TGraphAsymmErrors * graph = dbtable->GetGraph("err","Q2");

    xmin  = ( graph->GetX() )[TMath::LocMin(graph->GetN(),graph->GetX())];
    xmax  = ( graph->GetX() )[TMath::LocMax(graph->GetN(),graph->GetX())];
    ymin  = ( graph->GetY() )[TMath::LocMin(graph->GetN(),graph->GetY())];
    ymax  = ( graph->GetY() )[TMath::LocMax(graph->GetN(),graph->GetY())];

    c->Clear();
    c->Divide(2,1);
    c->GetPad(1)->SetPad("mplots_pad","",0.01,0.25,0.99,0.99);
    c->GetPad(2)->SetPad("legend_pad","",0.01,0.01,0.99,0.24);
    c->GetPad(1)->SetFillColor(0);
    c->GetPad(1)->SetBorderMode(0);
    c->GetPad(2)->SetFillColor(0);
    c->GetPad(2)->SetBorderMode(0);
    c->GetPad(1)->cd();
    c->GetPad(1)->SetBorderMode(0);
  
    c->GetPad(1)->SetLogx();
    //c->GetPad(1)->SetLogy();

    TH1F * hframe = (TH1F*) c->GetPad(1)->DrawFrame(.5*xmin, .4*ymin, 1.2*xmax, 1.42*ymax);
    hframe->Draw();

    MultiGraph * mgraph = dbtable->GetMultiGraph("err","Q2");
    for(unsigned int igraph = 0; igraph < mgraph->NGraphs(); igraph++) {
       mgraph->GetGraph(igraph)->Draw("P");
    }

    // Superimpose GENIE prediction
    //
    // ...

    c->GetPad(2)->cd();
    TLegend * legend = new TLegend(0.01, 0.01, 0.99, 0.99);
    legend->SetFillColor(0);
    mgraph->FillLegend("LP", legend);
    legend->SetTextSize(0.08);
    legend->Draw();
    
    c->GetPad(2)->Update();
    c->Update();        
  }


  ps->Close();
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
// Get the requested structure function data from the NuVld MySQL dbase
//.................................................................................
DBQueryString FormQuery(
           const char * key_list, float Q2min, float Q2max, float xmin, float xmax)
{
// forms a DBQueryString for extracting structure function data from the input 
// key-list and for the Q2 range

  ostringstream query_string;

  query_string << "KEY-LIST:" << key_list;
  query_string << "$CUTS:Q2min=" << Q2min << ";Q2max=" << Q2max
               << ";xmin=" << xmin << ";xmax=" << xmax
               << "$DRAW_OPT:none$DB-TYPE:SF";

  DBQueryString query(query_string.str());

  return query;
}
//_________________________________________________________________________________
DBTable<DBSFTableRow> * GetNuVldData(int iset)
{
  DBTable<DBSFTableRow> * dbtable = new DBTable<DBSFTableRow>;

  const char * keylist = kSFKeyList[iset];
  float        xmin    = kSFx[iset] - 0.001;
  float        xmax    = kSFx[iset] + 0.001;

  DBQueryString query = FormQuery(keylist, kSFQ2min, kSFQ2max, xmin, xmax);
  
  assert( dbi->FillTable(dbtable, query) == eDbu_OK );

  return dbtable;
}
//_________________________________________________________________________________
// Parsing command-line arguments, check/form filenames, etc
//.................................................................................
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing commad line arguments";

  // get input GENIE cross section file
  try {
    gOptSfFilename_curr = utils::clap::CmdLineArgAsString(argc,argv,'f');
    bool ok = CheckRootFilename(gOptSfFilename_curr.c_str());
    if(!ok) {
      PrintSyntax();
      exit(1);
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      PrintSyntax();
      exit(1);
    }
  }

  // get [reference] input GENIE cross section file
  try {
    gOptSfFilename_ref0 = utils::clap::CmdLineArgAsString(argc,argv,'r');
    bool ok = CheckRootFilename(gOptSfFilename_ref0.c_str());
    if(!ok) {
      PrintSyntax();
      exit(1);
    }
    gOptHaveRef = true;
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gvldtest", pNOTICE) << "No reference cross section file";
      gOptHaveRef = false;
    }
  }

  // check whether to compare with data
  gOptCmpWithData = genie::utils::clap::CmdLineArgAsBool(argc,argv,'d');

  if(gOptCmpWithData) {

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
       PrintSyntax();
       exit(1);
     }
   }

   // get DB passwd
   try {
     gOptDbPasswd = utils::clap::CmdLineArgAsString(argc,argv,'p');
   } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
       PrintSyntax();
       exit(1);
     }
   }
  } // -d enabled?
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gNuSampleTest  -f sample.root [-n nev] [-r reference_sample.root]\n";
}
//_________________________________________________________________________________
bool CheckRootFilename(string filename)
{
  if(filename.size() == 0) return false;

  bool is_accessible = ! (gSystem->AccessPathName(filename.c_str()));
  if (!is_accessible) {
   LOG("gvldtest", pERROR)
       << "The input ROOT file [" << filename << "] is not accessible";
   return false;
  }
  return true;
}
//_________________________________________________________________________________
string OutputFileName(string inpname)
{
// Builds the output filename based on the name of the input filename
// Perfors the following conversion: name.root -> name.nuxsec_test.ps

  unsigned int L = inpname.length();

  // if the last 4 characters are "root" (ROOT file extension) then
  // remove them
  if(inpname.substr(L-4, L).find("root") != string::npos) {
    inpname.erase(L-4, L);
  }

  ostringstream name;
  name << inpname << "nuxsec_test.ps";

  return gSystem->BaseName(name.str().c_str());
}
//_________________________________________________________________________________

