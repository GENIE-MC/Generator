//____________________________________________________________________________
/*!

\program gsplt

\brief   GENIE utility program to produce cross section plots from the XML
         cross section spline data that GENIE MC drivers can load and / or
         create by themselves at job initialization.

         Syntax :
             gsplt -f xml_file -p neutrino_pdg -t target_pdg [-e emax] [-o root_file]

         Options :
           []  denotes an optional argument
           -f  the input XML file containing the cross section spline data
           -p  the neutrino pdg code
           -t  the target pdg code (format: 1aaazzz000)
           -e  the maximum energy (in generated plots -- use it to zoom at low E)
           -o  if an output ROOT file is specified the splines data will be 
               saved there in a hierarchical directory structure.
               Note that the input ROOT file, if already exists, would not be 
               recreated and the output splines would be appended to the ones
               already stored.

         Example:

           gsplt -f ~/mydata/mysplines.xml -p 14 -t 1056026000 

           will load the cross section splines from the XML file mysplines.xml,
           then will select the cross section splines that are relevant to 
           numu+Fe56 and will generate cross section plots.
           The generated cross section plots will be saved in a postscript
           document named "xsec-splines-14-1056026000.ps"

         Notes:
         - To create the cross sections splines in XML format (for some target
           list or input geometry and for some input neutrino list) run the gmkspl
           GENIE application (see $GENIE/src/stdapp/gMakeSplines.cxx)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created December 15, 2005

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
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
#include <TH1F.h>

#include "Conventions/XmlParserStatus.h"
#include "Conventions/Units.h"
#include "EVGCore/InteractionList.h"
#include "EVGDrivers/GEVGDriver.h"
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
int    gOptNuPdgCode;    // neutrino PDG code
int    gOptTgtPdgCode;   // target PDG code

const int NP=300;
const int PsType = 111; // ps type: portrait

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- figuring out the energy, xsec range
  double Emin  = 0.01;
  double Emax  = gOptNuEnergy;
  double XSmax = -9999;
  double XSmin =  9999;

  //-- define some marker styles / colors
  const unsigned int kNMarkers = 5;
  const unsigned int kNColors  = 6;
  unsigned int markers[kNMarkers] = {20, 28, 29, 27, 3};
  unsigned int colors [kNColors]  = {1, 2, 4, 6, 8, 28};

  //-- load the cross section spline list
  XSecSplineList * splist = XSecSplineList::Instance();
  XmlParserStatus_t ist = splist->LoadFromXml(gOptXMLFilename);
  assert(ist == kXmlOK);

  //-- create an event genartion driver configured for the
  //   specified initial state (cross section splines will be
  //   accessed through that driver as in event generation mode)
  InitialState init_state(gOptTgtPdgCode, gOptNuPdgCode);
           
  GEVGDriver evg_driver;
  evg_driver.Configure(init_state);
  evg_driver.CreateSplines();
  evg_driver.CreateXSecSumSpline (100, Emin, Emax);

  //-- create a postscript document for saving cross section plpots

  TCanvas * c = new TCanvas("c","",20,20,500,800);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  TLegend * legend = new TLegend(0.01,0.01,0.99,0.99);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);

  ostringstream filename;
  filename << "xsec-splines-" << gOptNuPdgCode << "-" << gOptTgtPdgCode << ".ps";
  TPostScript * ps = new TPostScript(filename.str().c_str(), PsType);

  //-- get the list of interactions that can be simulated by the driver
  const InteractionList * ilist = evg_driver.Interactions();
  unsigned int nspl = ilist->size();

  //-- book enough space for xsec plots (last one is the sum)
  TGraph * gr[nspl+1];

  //-- loop over all the simulated interactions & create the cross section graphs
  InteractionList::const_iterator ilistiter = ilist->begin();
  unsigned int i=0;
  for(; ilistiter != ilist->end(); ++ilistiter) {
    
    const Interaction * interaction = *ilistiter;
    LOG("gsplt", pINFO) 
       << "Current interaction: " << interaction->AsString();

    //-- access the cross section spline
    const Spline * spl = evg_driver.XSecSpline(interaction);
    if(!spl) {
      LOG("gsplt", pWARN) << "Can't get spline for: " << interaction->AsString();
      exit(2);
    }

    //-- set graph color/style
    int icol = TMath::Min( i % kNColors,  kNColors-1  );
    int isty = TMath::Min( i / kNMarkers, kNMarkers-1 );
    int col = colors[icol];
    int sty = markers[isty];

    LOG("gsplt", pINFO) << "color = " << col << ", marker = " << sty;

    //-- export Spline as TGraph / set color & stype
    gr[i] = spl->GetAsTGraph(NP,true,true,1.,1./units::cm2);
    gr[i]->SetLineColor(col);
    gr[i]->SetMarkerColor(col);
    gr[i]->SetMarkerStyle(sty);
    gr[i]->SetMarkerSize(0.5);

    i++;
  }

  //-- now, get the sum...
  const Spline * splsum = evg_driver.XSecSumSpline();
  if(!splsum) {
      LOG("gsplt", pWARN) << "Can't get the cross section sum spline";
      exit(2);
  }
  gr[nspl] = splsum->GetAsTGraph(NP,true,true,1.,1./units::cm2);

  //-- figure out the minimum / maximum xsec in plotted range
  double x=0, y=0;
  for(int j=0; j<NP; j++) {
    gr[nspl]->GetPoint(j,x,y);
    XSmax = TMath::Max(XSmax,y);
  }
  XSmin = XSmax/2000.;
  XSmax = XSmax*1.2;

  LOG("gsplt", pINFO) << "Drawing frame: E    = (" << Emin  << ", " << Emax  << ")";
  LOG("gsplt", pINFO) << "Drawing frame: XSec = (" << XSmin << ", " << XSmax << ")";

  //-- ps output: add the 1st page with _all_ xsec spline plots
  //
  //c->Draw();
  TH1F * h = (TH1F*) c->DrawFrame(Emin, XSmin, Emax, XSmax);
  for(unsigned int i = 0; i <= nspl; i++) if(gr[i]) gr[i]->Draw("LP");
  h->GetXaxis()->SetTitle("Ev (GeV)");
  h->GetYaxis()->SetTitle("#sigma_{nuclear}/Ev (cm^{2}/GeV)");
  c->SetLogx();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  c->Update();

  //-- plot QEL xsecs only
  //
  h = (TH1F*) c->DrawFrame(Emin, XSmin, Emax, XSmax);
  i=0;
  for(ilistiter = ilist->begin(); ilistiter != ilist->end(); ++ilistiter) {    
    const Interaction * interaction = *ilistiter;
    if(interaction->ProcInfo().IsQuasiElastic()) {
        gr[i]->Draw("LP");
        TString spltitle(interaction->AsString());
        spltitle = spltitle.ReplaceAll(";",1," ",1);
        legend->AddEntry(gr[i], spltitle.Data(),"LP");
    }
    i++;
  }
  gr[nspl]->Draw("LP");
  legend->AddEntry(gr[nspl], "sum","LP");
  h->GetXaxis()->SetTitle("Ev (GeV)");
  h->GetYaxis()->SetTitle("#sigma_{nuclear}/Ev (cm^{2}/GeV)");
  c->SetLogx();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  c->Update();
  c->Clear();
  c->Range(0,0,1,1);
  legend->Draw();
  c->Update();

  //-- plot RES xsecs only
  //
  h = (TH1F*) c->DrawFrame(Emin, XSmin, Emax, XSmax);
  legend->Clear();
  i=0;
  for(ilistiter = ilist->begin(); ilistiter != ilist->end(); ++ilistiter) {    
    const Interaction * interaction = *ilistiter;
    if(interaction->ProcInfo().IsResonant()) {
        gr[i]->Draw("LP");
        TString spltitle(interaction->AsString());
        spltitle = spltitle.ReplaceAll(";",1," ",1);
        legend->AddEntry(gr[i], spltitle.Data(),"LP");
    }
    i++;
  }
  gr[nspl]->Draw("LP");
  legend->AddEntry(gr[nspl], "sum","LP");
  h->GetXaxis()->SetTitle("Ev (GeV)");
  h->GetYaxis()->SetTitle("#sigma_{nuclear}/Ev (cm^{2}/GeV)");
  c->SetLogx();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  c->Update();
  c->Clear();
  c->Range(0,0,1,1);
  legend->Draw();
  c->Update();

  //-- plot DIS xsecs only
  //
  h = (TH1F*) c->DrawFrame(Emin, XSmin, Emax, XSmax);
  legend->Clear();
  i=0;
  for(ilistiter = ilist->begin(); ilistiter != ilist->end(); ++ilistiter) {    
    const Interaction * interaction = *ilistiter;
    if(interaction->ProcInfo().IsDeepInelastic()) {
        gr[i]->Draw("LP");
        TString spltitle(interaction->AsString());
        spltitle = spltitle.ReplaceAll(";",1," ",1);
        legend->AddEntry(gr[i], spltitle.Data(),"LP");
    }
    i++;
  }
  gr[nspl]->Draw("LP");
  legend->AddEntry(gr[nspl], "sum","LP");
  h->GetXaxis()->SetTitle("Ev (GeV)");
  h->GetYaxis()->SetTitle("#sigma_{nuclear}/Ev (cm^{2}/GeV)");
  c->SetLogx();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  c->Update();
  c->Clear();
  c->Range(0,0,1,1);
  legend->Draw();
  c->Update();

  //-- plot COH xsecs only
  //
  h = (TH1F*) c->DrawFrame(Emin, XSmin, Emax, XSmax);
  legend->Clear();
  i=0;
  for(ilistiter = ilist->begin(); ilistiter != ilist->end(); ++ilistiter) {    
    const Interaction * interaction = *ilistiter;
    if(interaction->ProcInfo().IsCoherent()) {
        gr[i]->Draw("LP");
        TString spltitle(interaction->AsString());
        spltitle = spltitle.ReplaceAll(";",1," ",1);
        legend->AddEntry(gr[i], spltitle.Data(),"LP");
    }
    i++;
  }
  gr[nspl]->Draw("LP");
  legend->AddEntry(gr[nspl], "sum","LP");
  h->GetXaxis()->SetTitle("Ev (GeV)");
  h->GetYaxis()->SetTitle("#sigma_{nuclear}/Ev (cm^{2}/GeV)");
  c->SetLogx();
  c->SetLogy();
  c->SetGridx();
  c->SetGridy();
  c->Update();
  c->Clear();
  c->Range(0,0,1,1);
  legend->Draw();
  c->Update();

  //-- close the postscript document 
  ps->Close();

  //-- clean-up
  for(unsigned int j=0; j<=nspl; j++) { if(gr[j]) delete gr[j]; }
  delete c;
  delete ps;

  //-- check whether the splines will be saved in a ROOT file - if not, exit now
  bool save_in_root = gOptROOTFilename.size()>0;
  if(!save_in_root) return 0;
  
  //-- check whether the requested filename exists
  //   if yes, then open the file in 'update' mode 
  bool exists = !(gSystem->AccessPathName(gOptROOTFilename.c_str()));

  TFile * froot = 0;
  if(exists) froot = new TFile(gOptROOTFilename.c_str(), "UPDATE");
  else       froot = new TFile(gOptROOTFilename.c_str(), "RECREATE");
  assert(froot);

  //-- create directory structure
  ostringstream dptr;
  dptr << "xs_" << gOptNuPdgCode << "_" << gOptTgtPdgCode;
  ostringstream dtitle;
  dtitle << "Cross section splines for: "
         << gOptNuPdgCode << "+" << gOptTgtPdgCode;
  LOG("gsplt", pINFO) 
      << "Will store splines in ROOT TDir = " << dptr.str();

  TDirectory * topdir = 
         dynamic_cast<TDirectory *> (froot->Get(dptr.str().c_str()));
  if(topdir) {
     LOG("gsplt", pINFO) 
       << "Directory: " << dptr.str() << " already exists!! Exiting";
     froot->Close();
     delete froot;
     return 1;
  }

  topdir = froot->mkdir(dptr.str().c_str(),dtitle.str().c_str());
  topdir->cd();

  TDirectory * qeldir = topdir->mkdir("qel","qel");
  qeldir->cd();
  TDirectory * qelccdir = qeldir->mkdir("cc","cc");
  TDirectory * qelncdir = qeldir->mkdir("nc","nc");  
  for(ilistiter = ilist->begin(); ilistiter != ilist->end(); ++ilistiter) {    
    const Interaction * interaction = *ilistiter;
    if(!interaction->ProcInfo().IsQuasiElastic()) continue;
    const Spline * spl = evg_driver.XSecSpline(interaction);
    if(interaction->ProcInfo().IsWeakCC()) {}
    else {}

  }
  qelccdir->Write();
  qelncdir->Write();
  qeldir->Write();

  TDirectory * resdir = topdir->mkdir("res","res");
  resdir->Write();

  TDirectory * disdir = topdir->mkdir("dis","dis");
  disdir->Write();

  TDirectory * cohdir = topdir->mkdir("coh","coh");
  cohdir->Write();

  topdir->Write(dptr.str().c_str());

  if(froot) {
    froot->Close();
    delete froot;
  }

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
  // neutrino PDG code:
  try {
    LOG("gevgen", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified neutrino PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
  // target PDG code:
  try {
    LOG("gevgen", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
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

  // print the options you got from command line arguments
  LOG("gsplt", pINFO) << "Command line arguments:";
  LOG("gsplt", pINFO) << "  Input XML file    = " << gOptXMLFilename;
  LOG("gsplt", pINFO) << "  Neutrino PDG code = " << gOptNuPdgCode;
  LOG("gsplt", pINFO) << "  Target PDG code   = " << gOptTgtPdgCode;
  LOG("gsplt", pINFO) << "  Max neutrino E    = " << gOptNuEnergy;
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gsplt", pNOTICE)
      << "\n\n" << "Syntax:" << "\n"
      << "   gsplt -f xml_file -p neutrino_pdg -t target_pdg"
      << " [-e emax] [-o output_root_file]";
}
//____________________________________________________________________________
