//____________________________________________________________________________
/*!

\program gvld_e_res_xsec

\brief   Compares GENIE with resonance electro-production cross-section data.

         The data come from the JLab Hall C archive maintained at:
         https://hallcweb.jlab.org/resdata/database/
         A local copy may be found in:
         $GENIE/data/validation/eA/xsec/differential/res/
         The archive contains ~12k data points.

         Syntax:
           gvld_e_res_xsec 
                -r model_for_resonance_xsec
                -c model_for_dis_continuum_xsec
               [-d data_archive]
               [-s datasets_to_plot]

         Options:

           [] Denotes an optional argument.

           -r Specify GENIE resonance electro-production cross-section model.

           -c Specify GENIE electron DIS cross-section model.

           -d Full path to the electron cross-section data archive.
              By default, will pick the one at:
              $GENIE/data/validation/eA/xsec/differential/res/eRES.root

           -s Info on which datasets to plot. 
              By default, will pick the one at:
              $GENIE/data/validation/eA/xsec/differential/res/datasets.txt

         Example:

            % gvld_e_res_xsec 
                  -r genie::ReinSeghalRESPXSec/Default
                  -c genie::QPMDISPXSec/Default

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created Oct 16, 2009 

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TLegend.h>
#include <TBox.h>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/KineUtils.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"

using std::ostringstream;
using std::ifstream;
using std::string;
using std::vector;

using namespace genie;
using namespace genie::constants;

//.................................................................................
// Utility class to hold info on plotted datasets
//
class eResDataSetDescription
{
public:
  eResDataSetDescription(
     int tgtpdg, string expt, double E, double theta) :
    fTgtPdg (tgtpdg), 
    fExpt   (expt), 
    fE      (E), 
    fTheta  (theta)
  {   
  }
  eResDataSetDescription() 
  {   
  }
  int    TgtPdg   (void) const { return fTgtPdg; }
  int    TgtZ     (void) const { return pdg::IonPdgCodeToZ(fTgtPdg); }
  int    TgtA     (void) const { return pdg::IonPdgCodeToA(fTgtPdg); }
  string TgtName  (void) const { 
    // all data are either in Hydrogen or Deuterium
    if(fTgtPdg == 1000010010) return "Hydrogen";
    if(fTgtPdg == 1000010020) return "Deuterium";
    return "Other";
  }
  string Expt     (void) const { return fExpt;  }
  double E        (void) const { return fE;     }
  double Theta    (void) const { return fTheta; }
  string LabelTeX (void) const { 
    ostringstream label;
    label << fExpt << " (" << this->TgtName() << "), ";
    label << "E = " << fE << " GeV, ";
    label << "#theta = " << fTheta << "^{o}";
    return label.str();
  }
private:
  int    fTgtPdg;  //
  string fExpt;    //
  double fE;       //
  double fTheta;   //
};
//.................................................................................

//
// constants
//

// defaults
const char * kDefDataArchiveFilename = "data/validation/eA/xsec/differential/res/eRES.root";  
const char * kDefDataSetsFilename    = "data/validation/eA/xsec/differential/res/datasets.txt";  

// plot config
const int kNCx = 2; // number of columns in TCanvas::Divide()
const int kNCy = 2; // number of rows    in TCanvas::Divide()

// number of resonances and GENIE resonance IDs
const int kNRes=18;  
Resonance_t kResId[kNRes] = {
   kP33_1232, kS11_1535, kD13_1520, kS11_1650,
   kD13_1700, kD15_1675, kS31_1620, kD33_1700,
   kP11_1440, kP33_1600, kP13_1720, kF15_1680,
   kP31_1910, kP33_1920, kF35_1905, kF37_1950,
   kP11_1710, kF17_1970 
};

// current program draws predictions only for the explicit resonance-production
// model at W<Wcut
const bool kDrawHatchcedScalingRegion = true; 

const double kWcut = 1.7; // Wcut from UserPhysicsOptions.xml

//
// globals
//

string gOptDataArchiveFilename = ""; // -d command-line argument
string gOptDataSetsFilename    = ""; // -s command-line argument
string gOptRESModelName        = ""; // -r command-line argument
string gOptDISModelName        = ""; // -c command-line argument

TFile *        gResDataFile  = 0;
TTree *        gResDataTree  = 0;
TPostScript *  gPS           = 0;
TCanvas *      gC            = 0;

const XSecAlgorithmI * gRESXSecModel = 0;
const XSecAlgorithmI * gDISXSecModel = 0;

vector<eResDataSetDescription *> gDataSets; // list of plotted datasets

//
// function prototypes
//

void             Init               (void);
void             End                (void);
TGraphErrors *   Data               (unsigned int iset);
vector<TGraph *> Model              (unsigned int iset, unsigned int imodel);
void             Draw               (unsigned int iset);
void             GetCommandLineArgs (int argc, char ** argv);
void             PrintSyntax        (void);


//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);

  Init();

  // loop over data sets and plot data and corresponding GENIE predictions
  for(unsigned int iset = 0; iset < gDataSets.size(); iset++) 
  {
    LOG("gvldtest", pNOTICE) 
      << "Producing plots for: " << gDataSets[iset]->LabelTeX();
    Draw(iset);
  }

  End();

  LOG("gvldtest", pINFO)  << "Done!";
  return 0;
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("gvldtest", pNOTICE) << "Initializing...";

  // Set GENIE style
  utils::style::SetDefaultStyle();

  //
  // Get TTree with electron-scattering data
  //
  if( ! utils::system::FileExists(gOptDataArchiveFilename) ) {
      LOG("gvldtest", pFATAL) 
         << "Can not find file: " << gOptDataArchiveFilename;
      gAbortingInErr = true;
      exit(1);
  }
  gResDataFile = new TFile(gOptDataArchiveFilename.c_str(),"read");  
  gResDataTree = (TTree *) gResDataFile->Get("resnt");
  if(!gResDataTree) {
      LOG("gvldtest", pFATAL) 
         << "Can not find TTree `resnt' in file: " << gOptDataArchiveFilename;
      gAbortingInErr = true;
      exit(1);
  }

  //
  // Read information on which data-sets to plot
  //
  LOG("gvldtest", pDEBUG) 
    << "Reading dataset summary info from: " << gOptDataSetsFilename;
  ifstream summary_file(gOptDataSetsFilename.c_str());
  if (!summary_file.good() ) {
      LOG("gvldtest", pFATAL) 
         << "Can't open data summary file: " << gOptDataSetsFilename;
      gAbortingInErr = true;
      exit(1);
  }
  while(1) {
      // skip header lines staring with #
      if(summary_file.peek() == '#') {
         summary_file.ignore(1000, '\n');
       } else {
         int    target = 0;
         string expt   = "";
	 double E      = 0;
         double theta  = 0;
         summary_file >> target >> expt >> E >> theta;
         summary_file.ignore(1000, '\n');
         if(summary_file.eof()) break;            

         LOG("gvldtest", pDEBUG) 
            << "target: " << target << ", experiment: " << expt 
            << ", E = " << E << " GeV, theta = " << theta << " deg";

         eResDataSetDescription * dataset = 
                 new eResDataSetDescription(target,expt,E,theta);
         gDataSets.push_back(dataset);
       }
  }
  summary_file.close();
  LOG("gvldtest", pNOTICE)
     << "Read "  << gDataSets.size() << " datasets";


  //
  // Get cross-section models
  //
  AlgFactory * algf = AlgFactory::Instance();
  gRESXSecModel = 0;
  if(gOptRESModelName != "none"){
     vector<string> modelv = utils::str::Split(gOptRESModelName,"/");
     assert(modelv.size()==2);
     string model_name = modelv[0];
     string model_conf = modelv[1];
     gRESXSecModel =
        dynamic_cast<const XSecAlgorithmI *> (
            algf->GetAlgorithm(model_name, model_conf));
  }
  gDISXSecModel = 0;
  if(gOptDISModelName != "none"){
     vector<string> modelv = utils::str::Split(gOptDISModelName,"/");
     assert(modelv.size()==2);
     string model_name = modelv[0];
     string model_conf = modelv[1];
     gDISXSecModel =
        dynamic_cast<const XSecAlgorithmI *> (
            algf->GetAlgorithm(model_name, model_conf));
  }

  // Create plot canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // Get local time to tag outputs
  string lt_for_filename   = utils::system::LocalTimeAsString("%02d.%02d.%02d_%02d.%02d.%02d");
  string lt_for_cover_page = utils::system::LocalTimeAsString("%02d/%02d/%02d %02d:%02d:%02d");


  // Create output postscript file
  string filename  = Form("genie-e_res_data_comp-%s.ps",lt_for_filename.c_str());
  gPS = new TPostScript(filename.c_str(), 111);

  // Add cover page
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText("GENIE comparison with (e,e') data on 1H and 2D in the resonance region");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(lt_for_cover_page.c_str());
  hdr.AddText(" ");
  hdr.Draw();
  gC->Update();
}
//_________________________________________________________________________________
void End(void)
{
  LOG("gvldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;

  gResDataFile->Close();
}
//_________________________________________________________________________________
vector<TGraph *> Model(unsigned int iset, unsigned int imodel)
{
// Get GENIE predictions for the `iset' dataset 

  vector<TGraph *> model;

  LOG("gvldtest", pNOTICE) 
    << "Getting GENIE prediction (model ID = " 
    << imodel << ", data set ID = " << iset << ")";

  bool calc_res = (gRESXSecModel != 0);
  bool calc_dis = (gDISXSecModel != 0);
  bool calc_tot = calc_res && calc_dis;

  double M  = kNucleonMass;
  double M2 = M*M;

  double E     = gDataSets[iset]->E();
  double theta = gDataSets[iset]->Theta();
  double costh = TMath::Cos(2*kPi*theta/360.);

  LOG("gvldtest", pINFO) 
     << " E = " << E 
     << ", theta = " << theta << " (cos(theta) = " << costh << ")";

  int Z = gDataSets[iset]->TgtZ();
  int A = gDataSets[iset]->TgtA();
  int N = A-Z;
  bool   tgt_has_p = (Z>0);
  bool   tgt_has_n = (N>0);
  double frac_p = (double) Z / (double) A;
  double frac_n = (double) N / (double) A;

  const int n = 500;

  double W2_array[n];
  double d2sigRES_dEpdOmega_array[n]; 
  double d2sigDIS_dEpdOmega_array[n];
  double d2sigTOT_dEpdOmega_array[n];

  double Epmin = 0.01;
  double Epmax = E;
  double dEp = (Epmax-Epmin)/(n-1);

  for(int i=0; i<n; i++) {
     double Ep = Epmin + i*dEp;
     double Q2 = 2*E*Ep*(1-costh);
     double W2 = M2 + 2*M*(E-Ep)-Q2;
     double W  = TMath::Sqrt( TMath::Max(0.,W2) );

     double x  = 0;
     double y  = 0;
     utils::kinematics::WQ2toXY(Ep,M,W,Q2,x,y);          

     d2sigRES_dEpdOmega_array[i] = 0; 
     d2sigDIS_dEpdOmega_array[i] = 0;
     d2sigTOT_dEpdOmega_array[i] = 0;

     if(W2 <= 0) {
        LOG("gvldtest", pDEBUG) 
          << "Ep = " << Ep << ", Q2 = " << Q2 << ", W2 = " << W2
          << "... Skipping point";
        W2_array[i] = 0.;
        continue;
     }
     W2_array[i] = W2;

     // Will calculate  d^2 sigma / dW dQ^2 and then convert to d^2sigma / dE' dOmega
     double jacobian = (E*Ep)*(M+2*E*(1-costh))/(kPi*W);

     //
     // Calculate resonance cross-section
     //

     if(calc_res) {

        double d2sigRES_dWdQ2   = 0.;
        double d2sigRES_dWdQ2_p = 0.;
        double d2sigRES_dWdQ2_n = 0.;

        //
        // e+p -> e+Resonance
        //
        if(tgt_has_p) {
           Interaction * ep_res = 
                Interaction::RESEM(1000010010, kPdgProton,  kPdgElectron, E);
           ep_res->KinePtr()->SetW (W);
           ep_res->KinePtr()->SetQ2(Q2);
           // loop over resonances
           for(int ires=0; ires<kNRes; ires++) {
              ep_res->ExclTagPtr()->SetResonance(kResId[ires]);
              double xsec = gRESXSecModel->XSec(ep_res,kPSWQ2fE) / units::nb;
              xsec = TMath::Max(0., xsec);
              d2sigRES_dWdQ2_p += xsec;
              LOG("gvldtest", pINFO) 
                  << "d2xsec_dWdQ2(ep; " << utils::res::AsString(kResId[ires])
                  << "; E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
                  << "; Ep = " << Ep << " GeV, theta = " << theta << " deg) = " 
                  << xsec << " nbarn/GeV^3";
           }//res
           delete ep_res;
        }//has_p

        //
        // e+n -> e+Resonance
        //
        if(tgt_has_n) {
           Interaction * en_res = 
                Interaction::RESEM(1000000010, kPdgNeutron, kPdgElectron, E);
           en_res->KinePtr()->SetW (W);
           en_res->KinePtr()->SetQ2(Q2);
           // loop over resonances
           for(int ires=0; ires<kNRes; ires++) {
              en_res->ExclTagPtr()->SetResonance(kResId[ires]);
              double xsec = gRESXSecModel->XSec(en_res,kPSWQ2fE) / units::nb;
              xsec = TMath::Max(0., xsec);
              d2sigRES_dWdQ2_n += xsec;
              LOG("gvldtest", pINFO) 
                  << "d2xsec_dWdQ2(en; " << utils::res::AsString(kResId[ires])
                  << "; E = " << E << " GeV, W = " << W << " GeV, Q2 = " << Q2 << " GeV^2"
                  << "; Ep = " << Ep << " GeV, theta = " << theta << " deg) = " 
                  << xsec << " nbarn/GeV^3";
           }//res
           delete en_res;
        }//has_n

        // sum over resonances and averaged over p,n
        d2sigRES_dWdQ2 = (frac_p * d2sigRES_dWdQ2_p + frac_n * d2sigRES_dWdQ2_n);

        // convert to d^2sigma / dE' dOmega and store
        double d2sigRES_dEpdOmega = jacobian * d2sigRES_dWdQ2;
        if(TMath::IsNaN(d2sigRES_dEpdOmega)) {
           LOG("gvldtest", pWARN) << "Got a NaN!";
           d2sigRES_dEpdOmega = 0;
        }
        d2sigRES_dEpdOmega_array[i] = TMath::Max(0., d2sigRES_dEpdOmega);

     }//calc_res

     //
     // Calculate DIS cross-section
     //

     if(calc_dis) {
        double d2sigDIS_dWdQ2   = 0.;
        double d2sigDIS_dWdQ2_p = 0.;
        double d2sigDIS_dWdQ2_n = 0.;
        if(tgt_has_p) {
           // Note: Not setting quark ID
           // If the quark ID is set, the code returns (eg for neutrino CC) the vq->lq' cross-section.
           // But if a quark ID is not specified then the code loops over all relevant
           // valence and sea quark species and returns (eg for neutrino CC) the vN->lX cross-section.
           Interaction * ep_dis = 
               Interaction::DISEM(1000010010, kPdgProton,  kPdgElectron, E); 
           ep_dis->KinePtr()->SetW (W);
           ep_dis->KinePtr()->SetQ2(Q2);
           ep_dis->KinePtr()->Setx (x);
           ep_dis->KinePtr()->Sety (y);
           // Note:
           d2sigDIS_dWdQ2_p = gDISXSecModel->XSec(ep_dis,kPSWQ2fE) / units::nb;
           d2sigDIS_dWdQ2_p = TMath::Max(0., d2sigDIS_dWdQ2_p);
           delete ep_dis;
        }
        if(tgt_has_n) {
           Interaction * en_dis = 
               Interaction::DISEM(1000000010, kPdgNeutron, kPdgElectron, E);
           en_dis->KinePtr()->SetW (W);
           en_dis->KinePtr()->SetQ2(Q2);
           en_dis->KinePtr()->Setx (x);
           en_dis->KinePtr()->Sety (y);
           d2sigDIS_dWdQ2_n = gDISXSecModel->XSec(en_dis,kPSWQ2fE) / units::nb;
           d2sigDIS_dWdQ2_n = TMath::Max(0., d2sigDIS_dWdQ2_n);
           delete en_dis;
        }
        d2sigDIS_dWdQ2 = (frac_p * d2sigDIS_dWdQ2_p + frac_n * d2sigDIS_dWdQ2_n);

        // convert to d^2sigma / dE' dOmega and store
        double d2sigDIS_dEpdOmega = jacobian * d2sigDIS_dWdQ2;
        if(TMath::IsNaN(d2sigDIS_dEpdOmega)) {
           LOG("gvldtest", pWARN) << "Got a NaN!";
           d2sigDIS_dEpdOmega = 0;
        }
        d2sigDIS_dEpdOmega_array[i] = TMath::Max(0., d2sigDIS_dEpdOmega);

     }//calc_dis

     //
     // Calculate total cross-section
     //
     if(calc_tot) {
        d2sigTOT_dEpdOmega_array[i] = 
           d2sigRES_dEpdOmega_array[i] +
           d2sigDIS_dEpdOmega_array[i]; 
     }
   
  }//i

  // Create graphs & store them in array
  if(calc_tot) {
     TGraph * gr = new TGraph(n,W2_array,d2sigTOT_dEpdOmega_array);
     gr->SetLineColor(kBlack);
     gr->SetTitle("Total");
     model.push_back(gr);
  }
  if(calc_res) {
     TGraph * gr = new TGraph(n,W2_array,d2sigRES_dEpdOmega_array);
     gr->SetLineColor(kRed);
     gr->SetTitle("Resonance");
     model.push_back(gr);
  }
  if(calc_dis) {
     TGraph * gr = new TGraph(n,W2_array,d2sigDIS_dEpdOmega_array);
     gr->SetLineColor(kBlue);
     gr->SetTitle("Non-resonance bkg / DIS");
     model.push_back(gr);
  }
  
  return model;
}
//_________________________________________________________________________________
TGraphErrors * Data(unsigned int iset)
{
  const double dE      = 1.0E-3;
  const double dtheta  = 2.5E-2;

  double E     = gDataSets[iset]->E();
  double theta = gDataSets[iset]->Theta();

  int Z = gDataSets[iset]->TgtZ();
  int A = gDataSets[iset]->TgtA();

  const char * selection = 
    Form("E > %f && E < %f && theta > %f && theta < %f && Z == %d && A == %d",
         E     - dE,
         E     + dE,
         theta - dtheta,
         theta + dtheta,
         Z,A);

  gResDataTree->Draw("W2:xsec:xsec_err", selection, "goff");
  
  int n = gResDataTree->GetSelectedRows();

  LOG("gvldtest", pNOTICE) 
    << "Found " << n << " data points in the xsec archive";

  if(n == 0) return 0; // return null graph

  // Data returned by TTree::Draw() are not necessarily ordered in W
  // Do the ordering here before building the graph
  int    *  idx = new int   [n];
  double *  xv  = new double[n];
  double *  yv  = new double[n];
  double *  dyv = new double[n];

  TMath::Sort(n,gResDataTree->GetV1(),idx,false);

  for(int i=0; i<n; i++) {
     int ii = idx[i];
     xv [i] = (gResDataTree->GetV1())[ii];
     yv [i] = (gResDataTree->GetV2())[ii];
     dyv[i] = (gResDataTree->GetV3())[ii];
  }

  TGraphErrors * gr = new TGraphErrors(n,xv,yv,0,dyv);
  utils::style::Format(gr, 1,1,1,1,8,0.4);

  delete [] idx;
  delete [] xv;
  delete [] yv;
  delete [] dyv;

  return gr;
}
//_________________________________________________________________________________
void Draw(unsigned int iset)
{
  // get all measurements for the current channel from the NuValidator MySQL dbase
  TGraphErrors * data = Data(iset);
  if(!data) return;

  // get the corresponding GENIE model prediction
  vector<TGraph *> model = Model(iset,0);

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

  TH1F * hframe = 0;
  double xmin =  9999999999, scale_xmin = 0.50;
  double xmax = -9999999999, scale_xmax = 1.05;
  double ymin =  9999999999, scale_ymin = 0.40;
  double ymax = -9999999999, scale_ymax = 1.30;
  // get W2 and d2sigma/dEdOmega range in the the data
  xmin  = ( data->GetX() )[TMath::LocMin(data->GetN(),data->GetX())];
  xmax  = ( data->GetX() )[TMath::LocMax(data->GetN(),data->GetX())];
  ymin  = ( data->GetY() )[TMath::LocMin(data->GetN(),data->GetY())];
  ymax  = ( data->GetY() )[TMath::LocMax(data->GetN(),data->GetY())];
  xmin  = TMath::Max(xmin, 0.5); // some data go very low 
  // take also into account the d2sigma/dEdOmega range in the the models
  // (for the W2 range of the dataset) just in case data and MC are not that similar...
  for(unsigned int imode=0; imode < model.size(); imode++) {
    TGraph * mm = model[imode];
    if(mm) {
      for(int k=0; k<mm->GetN(); k++) {
         double x = (mm->GetX())[k];
         if(x < xmin || x > xmax) continue;
         ymin = TMath::Min(ymin, (mm->GetY())[k]);
         ymax = TMath::Max(ymax, (mm->GetY())[k]);
      }//k
    }//mm
  }//imode
  LOG("gvldtest", pDEBUG) 
    << "Plot range:" 
    << "W^{2} = [" << xmin << ", " << xmax << "] GeV^{2}, "
    << "d^{2}#sigma / d#Omega dE = [" << ymin << ", " << ymax << "]  nb/sr/GeV";

  hframe = (TH1F*) gC->GetPad(iplot)->DrawFrame(
        scale_xmin*xmin, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
  hframe->Draw();
  hframe->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
  hframe->GetYaxis()->SetTitle("d^{2}#sigma / d#Omega dE (nb/sr/GeV)");

  // draw data and GENIE models
  data->Draw("P");
  for(unsigned int imode=0; imode < model.size(); imode++) {
    TGraph * mm = model[imode];
    if(!mm) continue;
    mm->Draw("L");
  }

  // add legend
  TLegend * legend = new TLegend(0.20, 0.75, 0.50, 0.85);
  legend->SetLineStyle(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.025);
  legend->SetHeader("GENIE");
  for(unsigned int imode=0; imode < model.size(); imode++) {
    TGraph * mm = model[imode];
    if(!mm) continue;
    legend->AddEntry(mm,mm->GetTitle(),"L");
  }
  legend->Draw();

  // scaling region
  TBox * scaling_region = 0;
  if(kDrawHatchcedScalingRegion) {
    double W2c = kWcut*kWcut;
    if(W2c > scale_xmin*xmin && W2c < scale_xmax*xmax) {
       scaling_region = new TBox(
           W2c, scale_ymin*ymin, scale_xmax*xmax, scale_ymax*ymax);
//       scaling_region->SetFillColor(kRed);
//       scaling_region->SetFillStyle(3905);
       scaling_region->SetFillStyle(0);
       scaling_region->SetLineColor(kRed);
       scaling_region->SetLineStyle(kDashed);
       scaling_region->Draw();
    }
  }

  // some data show the elastic peak - mark the are to avoid confusion
  if(xmin < 1) {
    double Wm2 = 1.21; // between the QE and Delta peaks
    TBox * qe_peak = new TBox(
       scale_xmin*xmin, scale_ymin*ymin, Wm2, scale_ymax*ymax);
     qe_peak->SetFillColor(kBlue);
     qe_peak->SetFillStyle(3005);
     qe_peak->Draw();
  }

  // title
  TLatex * title = new TLatex(
     scale_xmin*xmin + 0.2*(scale_xmax*xmax-scale_xmin*xmin),
    1.01*scale_ymax*ymax,gDataSets[iset]->LabelTeX().c_str());
  title->SetTextSize(0.027);
  title->Draw();

  gC->GetPad(iplot)->Update();
  gC->Update();
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  if(parser.OptionExists('d')){
     string filename = parser.ArgAsString('d');
     gOptDataArchiveFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataArchiveFilename;
        gOptDataArchiveFilename = filename;
     } else { 
        LOG("gvldtest", pFATAL) 
          << "\n Please make sure that $GENIE is defined, or use the -d option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  if(parser.OptionExists('s')){
     string filename = parser.ArgAsString('s');
     gOptDataSetsFilename = filename;
  } else {
     if(gSystem->Getenv("GENIE")) {
        string base_dir = string( gSystem->Getenv("GENIE") );
        string filename = base_dir + "/" + kDefDataSetsFilename;
        gOptDataSetsFilename = filename;
     } else { 
        LOG("gvldtest", pFATAL) 
          << "\n Please make sure that $GENIE is defined, or use the -s option"
          << "\n You didn't specify a data file and I can not pick the default one either";
        gAbortingInErr = true;
        exit(1);
     }
  }

  // Get GENIE model names to be used
  if(parser.OptionExists('r')){
     gOptRESModelName = parser.ArgAsString('r');
  }   
  if(parser.OptionExists('c')){
     gOptDISModelName = parser.ArgAsString('c');
  }

}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" 
    << "\n"
    << " gvld_e_res_xsec \n"
    << "       -r model_for_resonance_xsec \n"
    << "       -c model_for_dis_continuum_xsec \n"
    << "      [-d data_archive] \n"
    << "      [-s datasets_to_plot] \n";
}
//_________________________________________________________________________________
