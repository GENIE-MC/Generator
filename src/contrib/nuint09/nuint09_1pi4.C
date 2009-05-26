//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// 1PI.4:
// 1pi+ cross section at E_nu= 1.0 and 1.5 GeV as a function of the hadronic invarant mass.
//
// Inputs:
// - sample id: 0 (nu_mu+C12), 1 (nu_mu+O16), 2 (nu_mu+Fe56)
// - single pion source: 0 (all), 1 (P33(1232) resonance only), 2 (resonances only)
// - stage: 0 -> primary Xpi+, 1 -> final state Xpi+    
//
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
//
#include <iomanip>

//
// consts
//
const int kNSamples       = 3;
const int kNWCur          = 1;
const int kNEnergies      = 2;
const int kNRunsPerCase   = 5;
const int kNEvtPerRun     = 100000;

const char * kLabel[kNSamples] = 
{  
 // available samples
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};
const int kRunNu1PI4[kNSamples][kNWCur][kNEnergies][kNRunsPerCase] =
{
 /* indices : sample ; cc/nc ; energy */
 {
  /* 0,0,0 (nu_mu C12,  CC, 1.0 GeV) */ { {900200, 900201, 900202, 900203, 900204},
  /* 0,0,1 (nu_mu C12,  CC, 1.5 GeV) */   {900300, 900301, 900302, 900303, 900304} }
 },
 {
  /* 1,0,0 (nu_mu O16,  CC, 1.0 GeV) */ { {910200, 910201, 910202, 910203, 910204},
  /* 1,0,1 (nu_mu O16,  CC, 1.5 GeV) */   {910300, 910301, 910302, 910303, 910304} }
 },
 {
  /* 2,0,0 (nu_mu Fe56, CC, 1.0 GeV) */ { {920200, 920201, 920202, 920203, 920204},
  /* 2,0,1 (nu_mu Fe56, CC, 1.5 GeV) */   {920300, 920301, 920302, 920303, 920304} }
 }
};

void nuint09_1pi4(int isample, int single_pion_sources=0, int stage=1)
{
 cout << " ***** running: 1PI.4" << endl;

  if(isample<0 || isample >= kNSamples) return;
  if(single_pion_sources<0 || single_pion_sources>2) return;

  const char * label = kLabel[isample];

  // get cross section graph  

  TFile fsig("./cc1pip_tmp.root","read"); // generated at [1PI.1]
  TGraph * sig_graph_cc1pip = (TGraph*) fsig.Get("CC1pip");

  // range & spacing

  const int    nW   = 60;
  const double Wmin =  0.80;
  const double Wmax =  2.50;

  // create output stream

  ostringstream out_filename;
  out_filename << label;

  if      (single_pion_sources==0) out_filename << ".1pi_4a.";
  else if (single_pion_sources==1) out_filename << ".1pi_4b.";
  else if (single_pion_sources==2) out_filename << ".1pi_4c.";

  if(stage==0) out_filename << "no_FSI.";

  out_filename << label << "dsig1pi_dW.data";

  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [1PI.4]:" << endl;
  out_stream << "#  1pi+ cross section at E_nu= 1.0 and 1.5 GeV as a function of the hadronic invariant mass W" << endl;
  if(stage==0) { 
    out_stream << "#  ***** NO FSI: The {X pi+} state is a primary hadronic state" << endl;
  }
  if(single_pion_sources==0) {
     out_stream << "#  1pi sources: All" << endl;
  }
  else if(single_pion_sources==1) {
     out_stream << "#  1pi sources: P33(1232) resonance only" << endl;
  }
  else if(single_pion_sources==2) {
     out_stream << "#  1pi sources: All resonances only" << endl;
  }
  out_stream << "#  " << endl;
  out_stream << "#  Note:" << endl;
  out_stream << "#   - W in GeV, linear spacing between Wmin = " << Wmin << " GeV, Wmax = " << Wmax << " GeV "  << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 / GeV " << endl;
  out_stream << "#   - quoted cross section is nuclear cross section divided with number of nucleons A" << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  W  |   dsig(numu A -> mu- 1pi+ X; Enu = 1.0 GeV)   |  dsig(numu A -> mu- 1pi+ X; Enu = 1.5 GeV)  | "  << endl;

  out_stream << setiosflags(ios::fixed) << setprecision(6);

  //
  // load event data
  // 
    
  TChain * chain = new TChain("gst");
  
  // loop over CC/NC cases
  for(int iwkcur=0; iwkcur<kNWCur; iwkcur++) {
    // loop over energies
    for(int ie=0; ie<kNEnergies; ie++) {
       // loop over runs for current case
       for(int ir=0; ir<kNRunsPerCase; ir++) {
          // build filename
          ostringstream filename;
          int run_number = kRunNu1PI4[isample][iwkcur][ie][ir];
          filename << "../gst/gntp." << run_number << ".gst.root";
          // add to chain
          cout << "Adding " << filename.str() << " to event chain" << endl;
          chain->Add(filename.str().c_str());
       }
    }
  } 

  // 
  // get CC1pi+ cross sections at given energies for normalization purposes
  //
  double sig_cc1pip_1000MeV = sig_graph_cc1pip -> Eval (1.0);
  double sig_cc1pip_1500MeV = sig_graph_cc1pip -> Eval (1.5);
  
  //
  // book histograms
  //
  TH1D * hst_dsig_dW_1000MeV = new TH1D("hst_dsig_dW_1000MeV","dsig/dW, numu A -> mu- 1pi+ X, Enu=1.0 GeV", nW, Wmin, Wmax);
  TH1D * hst_dsig_dW_1500MeV = new TH1D("hst_dsig_dW_1500MeV","dsig/dW, numu A -> mu- 1pi+ X, Enu=1.5 GeV", nW, Wmin, Wmax);

  //
  // fill histograms
  //
  if(stage==1) {
     if(single_pion_sources==0) {
       // all sources
       chain->Draw("Ws>>hst_dsig_dW_1000MeV","cc&&Ev>0.99&&Ev<1.01&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
       chain->Draw("Ws>>hst_dsig_dW_1500MeV","cc&&Ev>1.49&&Ev<1.51&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
     }
     else if(single_pion_sources==1) {
       // P33(1232) only
       chain->Draw("Ws>>hst_dsig_dW_1000MeV","cc&&resid==0&&res&&Ev>0.99&&Ev<1.01&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
       chain->Draw("Ws>>hst_dsig_dW_1500MeV","cc&&resid==0&&res&&Ev>1.49&&Ev<1.51&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
     }
     else if(single_pion_sources==2) {
       // all resonances only
       chain->Draw("Ws>>hst_dsig_dW_1000MeV","cc&&res&&Ev>0.99&&Ev<1.01&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
       chain->Draw("Ws>>hst_dsig_dW_1500MeV","cc&&res&&Ev>1.49&&Ev<1.51&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
     }
  }
  else if(stage==0) {
     if(single_pion_sources==0) {
       // all sources
       chain->Draw("Ws>>hst_dsig_dW_1000MeV","cc&&Ev>0.99&&Ev<1.01&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
       chain->Draw("Ws>>hst_dsig_dW_1500MeV","cc&&Ev>1.49&&Ev<1.51&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
     }
     else if(single_pion_sources==1) {
       // P33(1232) only
       chain->Draw("Ws>>hst_dsig_dW_1000MeV","cc&&resid==0&&res&&Ev>0.99&&Ev<1.01&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
       chain->Draw("Ws>>hst_dsig_dW_1500MeV","cc&&resid==0&&res&&Ev>1.49&&Ev<1.51&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
     }
     else if(single_pion_sources==2) {
       // all resonances only
       chain->Draw("Ws>>hst_dsig_dW_1000MeV","cc&&res&&Ev>0.99&&Ev<1.01&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
       chain->Draw("Ws>>hst_dsig_dW_1500MeV","cc&&res&&Ev>1.49&&Ev<1.51&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
     }
  }
                
  //
  // normalize
  //
  double norm_cc1pip_1000MeV = hst_dsig_dW_1000MeV -> Integral("width") / sig_cc1pip_1000MeV;
  double norm_cc1pip_1500MeV = hst_dsig_dW_1500MeV -> Integral("width") / sig_cc1pip_1500MeV;
  
  if (norm_cc1pip_1000MeV > 0) hst_dsig_dW_1000MeV -> Scale(1./norm_cc1pip_1000MeV);
  if (norm_cc1pip_1500MeV > 0) hst_dsig_dW_1500MeV -> Scale(1./norm_cc1pip_1500MeV);

  for(int i=1; i <= hst_dsig_dW_1000MeV->GetNbinsX(); i++) {

     double W = hst_dsig_dW_1000MeV -> GetBinCenter(i);

     double dsig_dW_1000MeV = hst_dsig_dW_1000MeV -> GetBinContent(i);
     double dsig_dW_1500MeV = hst_dsig_dW_1500MeV -> GetBinContent(i);

     dsig_dW_1000MeV = TMath::Max(0., dsig_dW_1000MeV);
     dsig_dW_1500MeV = TMath::Max(0., dsig_dW_1500MeV);

     out_stream << setw(15) << W
                << setw(15) << dsig_dW_1000MeV
                << setw(15) << dsig_dW_1500MeV
                << endl;
  }

  out_stream.close();

  // visual inspection
  //TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  //...
  //c1->Update();

  delete chain;
  fsig.Close();
}
