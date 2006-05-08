//____________________________________________________________________________
/*!

\program gsplt

\brief   GENIE utility program to produce cross section plots from the XML
         cross section spline data that GENIE MC drivers can load and / or
         create by themselves at job initialization.

         Syntax :
             gsplt -f xml_file -t target_pdg [-e emax] [-o output_root_file]

         Options :
           -f  the input XML file containing the cross section spline data
           -t  a target pdg code (format: 1aaazzz000)
           -e  the maximum energy (in plots)
           -o  if an output ROOT file is specified the splines data will be 
               saved there in a hierarchical directory structure.
               Note that the input ROOT file, if already exists, would not be 
               recreated and the output splines would be appended to the ones
               already stored.

         Example:

           gsplt -f ~/mydata/mysplines.xml -t 1056026000 -e 78.4

           will load the cross section splines from the XML file mysplines.xml,
           then will select the cross section splines that refer to an iron
           target (A=56,Z=26) and will plot xsec/E vs E up to 78.4 GeV (or the
           maximum energy in the loaded splines, if it is less that 78.4 GeV)
           The generated cross section plots will be saved in a postscript
           document named "xsec-splines-1056026000.ps"

         Notes:
         - To create the cross sections splines in XML format (for some target
           list or input geometry and for some input neutrino list) run the gmkspl
           GENIE application (see $GENIE/src/stdapp/gMakeSplines.cxx)
         - The event generation drivers can also generate cross section splines
           or complete a loaded cross section spline list with all the splines
           that would be needed but were not loaded. See at the event generation
           drivers for more information.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 15, 2005

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <string>
#include <sstream>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TPostScript.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TString.h>

#include "Conventions/XmlParserStatus.h"
#include "Conventions/Units.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGUtils.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;

//Prototypes:
void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//User-specified options:
string gOptXMLFilename;  // input XML filename
string gOptROOTFilename; // output ROOT filename
double gOptNuEnergy;     // Ev(max)
int    gOptTgtPdgCode;   // target PDG code

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- load the cross section spline list
  XSecSplineList * splist = XSecSplineList::Instance();
  XmlParserStatus_t ist = splist->LoadFromXml(gOptXMLFilename);
  assert(ist == kXmlOK);

  //-- get all available spline keys
  const vector<string> * keyv = splist->GetSplineKeys();
  vector<string>::const_iterator kiter;

  //-- figuring out the energy, xsec range
  double Emin  = 0.1;
  double Emax  = gOptNuEnergy;
  double XSmax = -9999;
  double XSmin =  9999;

  //-- string to match for selecting splines for the input target only
  ostringstream tgtmatch;
  tgtmatch << "tgt:" << gOptTgtPdgCode;

  //-- get all splines and draw them
  TCanvas * c = 0;
  TH1F *    h = 0;

  //-- create a postscript document
  ostringstream filename;
  filename << "xsec-splines-" << gOptTgtPdgCode << ".ps";
  TPostScript * ps = new TPostScript(filename.str().c_str(), 111);

  //-- create plot legend
  TLegend * legend = new TLegend(0.85,0.2,1.00,0.9);
  legend->SetFillColor(0);

  //-- book enough space for xsec plots
  TGraph * gr[keyv->size()];

  //-- define some marker styles / colors
  const int kNMarkers = 5;
  const int kNColors  = 6;
  int markers[kNMarkers] = {20, 28, 29, 27, 3};
  int colors [kNColors]  = {1, 2, 4, 6, 8, 28};

  //-- check whether to append the splines into a ROOT file
  TFile * froot = 0;
  bool save_in_root = gOptROOTFilename.size()>0;

  //-- check whether the file already exists in which case open
  //   it in APPEND mode, otherwise open it in RECREATE mode
  if(save_in_root) {
    bool exists = !(gSystem->AccessPathName(gOptROOTFilename.c_str()));

    if(exists) froot = new TFile(gOptROOTFilename.c_str(), "UPDATE");
    else       froot = new TFile(gOptROOTFilename.c_str(), "RECREATE");
    assert(froot);
  }

  //-- if the splines are to be saved in an output root file, cd into
  //   a directory with the nuclear target PDG code (create it if it
  //   does not exist)
  TDirectory * dirTgt = 0;
  if(save_in_root) {
    ostringstream dptr;
    dptr << "tgt_" << gOptTgtPdgCode;
    ostringstream dtitle;
    dtitle << "Cross section splines for target: " << gOptTgtPdgCode;

    LOG("gsplt", pINFO) 
                 << "Will store splines in ROOT TDir = " << dptr.str();

    dirTgt = (TDirectory*) froot->Get(dptr.str().c_str());
    if(!dirTgt) {
      dirTgt = new TDirectory(dptr.str().c_str(),dtitle.str().c_str());
      dirTgt->Write();
    }
  }

  //-- fill all cross section graphs for the input target
  int i=-1;
  for(kiter = keyv->begin(); kiter != keyv->end(); ++kiter) {

    i++;
    string key = *kiter;

    if(key.find(tgtmatch.str()) == string::npos) {
      gr[i] = 0;
      continue;
    }

    LOG("gsplt", pINFO) << "Drawing spline with key = \n" << key;

    //-- access the cross section spline
    const Spline * spl = splist->GetSpline(key);

    //-- set graph color/style
    int icol = TMath::Min( i % kNColors,  kNColors-1  );
    int isty = TMath::Min( i / kNMarkers, kNMarkers-1 );
    int col = colors[icol];
    int sty = markers[isty];

    LOG("gsplt", pINFO) << "color = " << col << ", marker = " << sty;

    //-- export Spline as TGraph / set color & stype
    int NP=300;
    gr[i] = spl->GetAsTGraph(NP,true,true,1.,1./units::cm2);
    gr[i]->SetLineColor(col);
    gr[i]->SetMarkerColor(col);
    gr[i]->SetMarkerStyle(sty);
    gr[i]->SetMarkerSize(0.5);

    //-- update maximum xsec in plot
    double x=0, y=0;
    for(int j=0; j<NP; j++) {
       gr[i]->GetPoint(j,x,y);
       XSmax = TMath::Max(XSmax,y);
    }

    //-- create legend entry for current graph
    ostringstream lgentry;
    lgentry << "spl-" << i;
    legend->AddEntry(gr[i], lgentry.str().c_str(),"LP");

    //-- save splines in ROOT file
    if(save_in_root) {

      dirTgt->cd(); 

      vector<string> keyv = utils::str::Split(key,"/");
      assert(keyv.size()==3);

      string algkey = keyv[0] + "/" + keyv[1];
      string intkey = keyv[2];

      TGraph * gr = spl->GetAsTGraph(4000);

      ostringstream grptrn;
      grptrn << "gr_" << i;

      // replace all ";" from the interaction key with spaces
      TString spltitle(intkey.c_str());
      spltitle = spltitle.ReplaceAll(";",1," ",1);
      spltitle.Prepend(", INT=");
      spltitle.Prepend(algkey.c_str());
      spltitle.Prepend("ALG=");
      gr->SetTitle(spltitle.Data());
      gr->Write(grptrn.str().c_str());
    }
  }
  XSmin = XSmax/300.;

  LOG("gsplt", pINFO) << "Saving plots in a postscript document";

  //-- ps output: add the 1st page with xsec spline plots
  c = new TCanvas("c","",20,20,500,500);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  c->Draw();

  LOG("gsplt", pINFO) << "Drawing frame: E    = (" << Emin  << ", " << Emax  << ")";
  LOG("gsplt", pINFO) << "Drawing frame: XSec = (" << XSmin << ", " << XSmax << ")";

  h = (TH1F*) c->DrawFrame(Emin, XSmin, Emax, XSmax);

  for(int i = 0; i < (int) keyv->size(); i++) if(gr[i]) gr[i]->Draw("LP");
  legend->Draw();

  h->GetXaxis()->SetTitle("Ev (GeV)");
  h->GetYaxis()->SetTitle("#sigma_{nuclear}/Ev (cm^{2}/GeV)");

  c->SetLogx();
  c->SetLogy();
  legend->Draw();
  c->Update();

  //-- ps output: create the 2nd page (detailed legend) with spline keys

  LOG("gsplt", pINFO) << "Creating legend page";

  ps->NewPage();

  delete c;
  c = new TCanvas("c","",20,20,500,500);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  c->Draw();

  c->Range(0,0,100,100);
  TPaveText * ptxt = new TPaveText(0,0,100,100);
  i = -1;
  for(kiter = keyv->begin(); kiter != keyv->end(); ++kiter) {
    i++;
    string key = *kiter;
    if(key.find(tgtmatch.str()) == string::npos) continue;
    ostringstream lgentry;
    lgentry << "spl-" << i << " : " << key;
    ptxt->InsertText(lgentry.str().c_str());
  }
  ptxt->Draw();
  c->Update();

  //-- close the postscript document and output ROOT file
  ps->Close();
  if(froot) froot->Close();

  //-- clean-up
  for(unsigned int j=0; j<keyv->size(); j++) { if(gr[j]) delete gr[j]; }
  delete c;
  delete ps;
  delete keyv;
  delete ptxt;
  if(froot) delete froot;

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gsplt", pINFO) << "Parsing commad line arguments";

  //input XML file name:
  try {
    LOG("gsplt", pINFO) << "Reading input XML filename";
    gOptXMLFilename = genie::utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gsplt", pFATAL) << "Unspecified input XML file!";
      PrintSyntax();
      exit(1);
    }
  }
  //output ROOT file name:
  try {
    LOG("gsplt", pINFO) << "Reading output ROOT filename";
    gOptROOTFilename = genie::utils::clap::CmdLineArgAsString(argc,argv,'o');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gsplt", pDEBUG)  
                      << "Unspecified ROOT file. Splines will not be saved.";
      gOptROOTFilename = "";
    }
  }
  //max neutrino energy
  try {
    LOG("gsplt", pINFO) << "Reading neutrino energy";
    gOptNuEnergy = genie::utils::clap::CmdLineArgAsDouble(argc,argv,'e');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gsplt", pDEBUG)  << "Unspecified Emax - Setting to 100 GeV";
      gOptNuEnergy = 100;
    }
  }
  //target PDG code:
  try {
    LOG("gsplt", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gsplt", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // print the options you got from command line arguments
  LOG("gsplt", pINFO) << "Command line arguments:";
  LOG("gsplt", pINFO) << "  Input XML file  = " << gOptXMLFilename;
  LOG("gsplt", pINFO) << "  Max neutrino E  = " << gOptNuEnergy;
  LOG("gsplt", pINFO) << "  Target PDG code = " << gOptTgtPdgCode;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gsplt", pNOTICE)
      << "\n\n" << "Syntax:" << "\n"
      << "   gsplt -f xml_file -t target_pdg [-e emax] [-o output_root_file]";
}
//____________________________________________________________________________
