//____________________________________________________________________________
/*!

\program gMAIDValidation

\brief   Generates plots MAID Form Factors for all fitted resonances
         Transverse and Longitudinal cross sections as a function of W are also 
	 reported 

\syntax  gMAIDValidation [-o output_name.root]

         -o :
          Specifies a name to be used in the output files.
          Default: maid_validation.root

        --W-min-XSec: minimum W value for cross section calculation
	--W-max-XSec: maximum W value for cross section calculation
	--W-binning-XSec: binning on W
        --Q2-min-XSec: maximum Q2 value for cross section calculation
        --Q2-max-XSec: minimum Q2 value for cross section calculation
	--Q2-binning-XSec: binning on Q2

\author  Julia Tena Vidal <jtenavidal \at tauex.tau.ac.il>
Tel Aviv University

\created 13 Mar 2023

\cpright Copyright (c) 2023-2033, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         
*/
//____________________________________________________________________________

#include <vector>
#include <string>

#include <TFile.h>
#include <TMath.h>
#include <TPostScript.h>
#include <TPavesText.h>
#include <TText.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TPaletteAxis.h>

#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/Style.h"
#include "Framework/Utils/RunOpt.h"
#include <map>
#include "Framework/Algorithm/AlgFactory.h"
#include "Physics/Resonance/XSection/RESVectFormFactorsI.h"
#include "Physics/Resonance/XSection/RESVectFormFactorsI.h"

using namespace std;
using namespace genie;
using namespace genie::utils;

// globals
string gOptOutFile = "maid_validation.root"; // -o argument
double gQ2_min_XSec = 0 ;
double gQ2_max_XSec = 4 ;
double gQ2_binning_XSec = 100 ;
double gW_min_XSec = 1.100 ;
double gW_max_XSec = 1.800 ;
double gW_bin_width = 50 ; 

TFile * goutput_file ; 

map<Resonance_t,TGraph*> gA12_p ;
map<Resonance_t,TGraph*> gA32_p ;
map<Resonance_t,TGraph*> gS12_p ;
map<Resonance_t,TGraph*> gA12_n ;
map<Resonance_t,TGraph*> gA32_n ;
map<Resonance_t,TGraph*> gS12_n ;
 
TCanvas * gC = 0;

RESVectFormFactorsI * vffmodel_p = 0 ; 
RESVectFormFactorsI * vffmodel_n = 0 ; 

// number of resonances and GENIE resonance IDs
const int kNRes = 18 ;  
Resonance_t kResId[kNRes] = {
  kP33_1232, kS11_1535, kD13_1520, kS11_1650,
  kD13_1700, kD15_1675, kS31_1620, kD33_1700,
  kP11_1440, kP33_1600, kP13_1720, kF15_1680,
  kP31_1910, kP33_1920, kF35_1905, kF37_1950,
  kP11_1710, kF17_1970 
};

// function prototypes
void GetCommandLineArgs (int argc, char ** argv);
bool InitializeMembers(void) ; 
void MakePlots     (void);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  utils::style::SetDefaultStyle();

  GetCommandLineArgs (argc,argv);   // Get command line arguments
  if( ! InitializeMembers() ) return 0 ; 
  MakePlots();   // Produce all output plots and fill output n-tuple

  return 0;
}
//_________________________________________________________________________________
void MakePlots (void)
{
  double Q2_width = ( gQ2_max_XSec - gQ2_min_XSec ) / gQ2_binning_XSec ; 
  vector<double> Q2_binning ; 
  vector<double> A12_p , A12_n, A32_p, A32_n, S12_p, S12_n ; 
  int lepton_pdg = kPdgElectron ; 

  // Proton case : 
  Interaction * inter_p = Interaction::RESEM(1000010010, kPdgProton, lepton_pdg);   
  Interaction * inter_n = Interaction::RESEM(1000000010, kPdgNeutron, lepton_pdg);   
  double W = 0 ; 
  double Q2 = 0;

  for ( int ires = 0 ; ires < kNRes ; ++ires ) { 
    Resonance_t ResID = kResId[ires] ; 
    inter_p->ExclTagPtr()->SetResonance(ResID) ; 
    inter_n->ExclTagPtr()->SetResonance(ResID) ; 
    A12_p.clear();
    A32_p.clear();
    S12_p.clear();
    A12_n.clear();
    A32_n.clear();
    S12_n.clear();
    Q2_binning.clear();

    for( unsigned int i = 0 ; i < gQ2_binning_XSec + 1 ; ++i ) { 
      Q2 = gQ2_min_XSec + i * Q2_width ; 
      inter_p->KinePtr()->SetQ2(Q2) ;
      inter_n->KinePtr()->SetQ2(Q2) ;
      inter_p->KinePtr()->SetW(W) ;
      inter_n->KinePtr()->SetW(W) ;

      RESVectFFAmplitude vffampl_p = vffmodel_p->Compute(*inter_p);
      A12_p.push_back( vffampl_p.AmplA12() ) ;
      A32_p.push_back( vffampl_p.AmplA32() ) ;
      S12_p.push_back( vffampl_p.AmplS12() ) ;

      RESVectFFAmplitude vffampl_n = vffmodel_n->Compute(*inter_n);
      A12_n.push_back( vffampl_n.AmplA12() ) ;
      A32_n.push_back( vffampl_n.AmplA32() ) ;
      S12_n.push_back( vffampl_n.AmplS12() ) ;
      Q2_binning.push_back( Q2 ) ; 
    }

    // Create TGraph 
    gA12_p[ResID] = new TGraph( Q2_binning.size(), &Q2_binning[0], &A12_p[0] ) ;
    gA32_p[ResID] = new TGraph( Q2_binning.size(), &Q2_binning[0], &A32_p[0] ) ;
    gS12_p[ResID] = new TGraph( Q2_binning.size(), &Q2_binning[0], &S12_p[0] ) ;

    gA12_n[ResID] = new TGraph( Q2_binning.size(), &Q2_binning[0], &A12_n[0] ) ;
    gA32_n[ResID] = new TGraph( Q2_binning.size(), &Q2_binning[0], &A32_n[0] ) ;
    gS12_n[ResID] = new TGraph( Q2_binning.size(), &Q2_binning[0], &S12_n[0] ) ;
    
    string res_name = utils::res::AsString(ResID) ;
    gA12_p[ResID]->SetTitle((res_name).c_str());
    gA12_p[ResID]->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
    gA12_p[ResID]->GetYaxis()->SetTitle("A_{1/2}^{p} #left[10^{-3} GeV^{-2}#right]");
    gA32_p[ResID]->SetTitle((res_name).c_str());
    gA32_p[ResID]->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
    gA32_p[ResID]->GetYaxis()->SetTitle("A_{3/2}^{p} #left[10^{-3} GeV^{-2}#right]");
    gS12_p[ResID]->SetTitle((res_name).c_str());
    gS12_p[ResID]->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
    gS12_p[ResID]->GetYaxis()->SetTitle("S_{1/2}^{p} #left[10^{-3} GeV^{-2}#right]");
    
    gA12_n[ResID]->SetTitle((res_name).c_str());
    gA12_n[ResID]->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
    gA12_n[ResID]->GetYaxis()->SetTitle("A_{1/2}^{n} #left[10^{-3} GeV^{-2}#right]");
    gA32_n[ResID]->SetTitle((res_name).c_str());
    gA32_n[ResID]->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
    gA32_n[ResID]->GetYaxis()->SetTitle("A_{3/2}^{n} #left[10^{-3} GeV^{-2}#right]");
    gS12_n[ResID]->SetTitle((res_name).c_str());
    gS12_n[ResID]->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
    gS12_n[ResID]->GetYaxis()->SetTitle("S_{1/2}^{n} #left[10^{-3} GeV^{-2}#right]");

    gA12_p[ResID]->Write(("A12_p_"+res_name).c_str()) ; 
    gA32_p[ResID]->Write(("A32_p_"+res_name).c_str()) ; 
    gS12_p[ResID]->Write(("S12_p_"+res_name).c_str()) ; 

    gA12_n[ResID]->Write(("A12_n_"+res_name).c_str()) ; 
    gA32_n[ResID]->Write(("A32_n_"+res_name).c_str()) ; 
    gS12_n[ResID]->Write(("S12_n_"+res_name).c_str()) ; 
    
    // Clean vectors
    A12_p.clear() ; 
    A32_p.clear() ; 
    S12_p.clear() ; 
    A12_n.clear() ; 
    A32_n.clear() ; 
    S12_n.clear() ; 
    
  }
  
  goutput_file->Close();

  return ; 
}
//_________________________________________________________________________________
bool InitializeMembers(void) {

  AlgFactory * algf = AlgFactory::Instance();

  vffmodel_p = dynamic_cast<RESVectFormFactorsI*> ( algf->AdoptAlgorithm("genie::MAIDRESVectFormFactorsEMp","Default") );
  vffmodel_n = dynamic_cast<RESVectFormFactorsI*> ( algf->AdoptAlgorithm("genie::MAIDRESVectFormFactorsEMn","Default") );

  if( ! vffmodel_p ) return false ; 
  if( ! vffmodel_n ) return false ; 

  goutput_file = TFile::Open(gOptOutFile.c_str(), "RECREATE") ; 
  return true ; 
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  // necessary for setting from whence it gets ModelConfiguration.xml
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  CmdLnArgParser parser(argc,argv);

  if(parser.OptionExists('o')){
    gOptOutFile = parser.Arg('o');
  }

  if(parser.OptionExists("W-min-XSec")) {
    string Wmin = parser.Arg("W-min_XSec");
    gW_min_XSec = atof(Wmin.c_str());
  }

  if(parser.OptionExists("Q2-min-XSec")) {
    string Q2min = parser.Arg("Q2-min_XSec");
    gQ2_min_XSec = atof(Q2min.c_str());
  }

  if(parser.OptionExists("Q2-max-XSec")) {
    string Q2max = parser.Arg("Q2-max_XSec");
    gQ2_max_XSec = atof(Q2max.c_str());
  }

  if(parser.OptionExists("Q2-binning-XSec")) {
    string Q2binning = parser.Arg("Q2-binning_XSec");
    gQ2_binning_XSec = atof(Q2binning.c_str());
  }

  if(parser.OptionExists("W-max-XSec")) {
    string Wmax = parser.Arg("W-max_XSec");
    gW_max_XSec = atof(Wmax.c_str());
  }

  if(parser.OptionExists("W-binning-XSec")) {
    string BinW = parser.Arg("W-binning_XSec");
    gW_bin_width = atof(BinW.c_str());
  }
  return ; 
}
//_________________________________________________________________________________

