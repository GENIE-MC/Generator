//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// COH.3:
// dSigma / dOmega_pi dKE_pi at E_nu= 0.5, 1.0, and 1.5 GeV
//
// Inputs:
// - sample id: 0 (nu_mu+C12), 1 (nu_mu+O16), 2 (nu_mu+Fe56)
//
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
//
#include <iomanip>

//
// consts
//
const int kNSamples       = 3;
const int kNWCur          = 2;
const int kNEnergies      = 3;
const int kNRunsPerCase   = 5;

const char * kLabel[kNSamples] = 
{  
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};

const int kRunNu[kNSamples][kNWCur][kNEnergies][kNRunsPerCase] = 
{  
 /* indices : sample ; cc/nc ; energy */
 {
  /* 0,0,0 (nu_mu C12,  CC, 0.5 GeV) */ { {100100, 100101, 100102, 100103, 100104},
  /* 0,0,1 (nu_mu C12,  CC, 1.0 GeV) */   {100200, 100201, 100202, 100203, 100204},
  /* 0,0,2 (nu_mu C12,  CC, 1.5 GeV) */   {100300, 100301, 100302, 100303, 100304} },
  /* 0,1,0 (nu_mu C12,  NC, 0.5 GeV) */ { {101100, 101101, 101102, 101103, 101104},
  /* 0,1,1 (nu_mu C12,  NC, 1.0 GeV) */   {101200, 101201, 101202, 101203, 101204},
  /* 0,1,2 (nu_mu C12,  NC, 1.5 GeV) */   {101300, 101301, 101302, 101303, 101304} }
 },
 {
  /* 1,0,0 (nu_mu O16,  CC, 0.5 GeV) */ { {110100, 110101, 110102, 110103, 110104},
  /* 1,0,1 (nu_mu O16,  CC, 1.0 GeV) */   {110200, 110201, 110202, 110203, 110204},
  /* 1,0,2 (nu_mu O16,  CC, 1.5 GeV) */   {110300, 110301, 110302, 110303, 110304} },
  /* 1,1,0 (nu_mu O16,  NC, 0.5 GeV) */ { {111100, 111101, 111102, 111103, 111104},
  /* 1,1,1 (nu_mu O16,  NC, 1.0 GeV) */   {111200, 111201, 111202, 111203, 111204},
  /* 1,1,2 (nu_mu O16,  NC, 1.5 GeV) */   {111300, 111301, 111302, 111303, 111304} }
 },
 {
  /* 2,0,0 (nu_mu Fe56, CC, 0.5 GeV) */ { {120100, 120101, 120102, 120103, 120104},
  /* 2,0,1 (nu_mu Fe56, CC, 1.0 GeV) */   {120200, 120201, 120202, 120203, 120204},
  /* 2,0,2 (nu_mu Fe56, CC, 1.5 GeV) */   {120300, 120301, 120302, 120303, 120304} },
  /* 2,1,0 (nu_mu Fe56, NC, 0.5 GeV) */ { {121100, 121101, 121102, 121103, 121104},
  /* 2,1,1 (nu_mu Fe56, NC, 1.0 GeV) */   {121200, 121201, 121202, 121203, 121204},
  /* 2,1,2 (nu_mu Fe56, NC, 1.5 GeV) */   {121300, 121301, 121302, 121303, 121304} }
 }
};


void nuint09_coh3(int isample)
{
  cout << " ***** running: COH.3" << endl;

  if(isample<0 || isample >= kNSamples) return;

  const char * label = kLabel[isample];

  // get cross section graphs

  TFile fsig("../sig/splines.root","read");
  TDirectory * sig_dir = (TDirectory *) fsig.Get(label);  

  TGraph * sig_graph_cohcc = (TGraph*) sig_dir->Get("coh_cc");
  TGraph * sig_graph_cohnc = (TGraph*) sig_dir->Get("coh_nc");

  // range & spacing

  const int    nKEpi    = 60;
  const double KEpimin  =  0.01;
  const double KEpimax  =  1.50;

  const int    ncosth   = 30;
  const double costhmin = -1;
  const double costhmax = +1;

  // create output stream

  ostringstream out_filename;
  out_filename << label << ".coh_3.d2sig_dKEpidOmega.data";
  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [COH.3]:" << endl;
  out_stream << "#  dSigma / dOmega_pi dKE_pi at E_nu= 0.5, 1.0, and 1.5 GeV " << endl;
  out_stream << "#  " << endl;
  out_stream << "#  Note:" << endl;
  out_stream << "#   - pion energies are _kinetic_ energies " << endl;
  out_stream << "#   - pion kinetic energy KE in GeV, between KEmin = " << KEpimin << " GeV, KEmax = " << KEpimax << " GeV "  << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 /GeV /sterad" << endl;
  out_stream << "#   - for coherent scattering we quote _nuclear_ cross section " << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  KE(pi)  |  cos(theta_pi)  |  sig(coh; CC; Ev=0.5GeV)  | sig(coh; NC; Ev=0.5GeV)  |  sig(coh; CC; Ev=1.0GeV)  | sig(coh; NC; Ev=1.0GeV)  |  sig(coh; CC; Ev=1.5GeV)  | sig(coh; NC; Ev=1.5GeV) " << endl;

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
          int run_number = kRunNu[isample][iwkcur][ie][ir];
          filename << "../gst/gntp." << run_number << ".gst.root";
          // add to chain
          cout << "Adding " << filename.str() << " to event chain" << endl;
          chain->Add(filename.str().c_str());
       }
    }
  }


  //
  // get COH pi cross sections at given energies for normalization purposes
  //
  double sig_cohcc_0500MeV = sig_graph_cohcc->Eval(0.5); 
  double sig_cohnc_0500MeV = sig_graph_cohnc->Eval(0.5); 
  double sig_cohcc_1000MeV = sig_graph_cohcc->Eval(1.0); 
  double sig_cohnc_1000MeV = sig_graph_cohnc->Eval(1.0); 
  double sig_cohcc_1500MeV = sig_graph_cohcc->Eval(1.5); 
  double sig_cohnc_1500MeV = sig_graph_cohnc->Eval(1.5); 

  //
  // book histograms (nevents in pion kinetic energy and cos(theta) bins)
  //
  TH2D * hst_d2sig_dKEpidOmg_cohcc_0500MeV = new TH2D("hst_d2sig_dKEpidOmg_cohcc_0500MeV","d^2sig / dKEpi dOmega, COH CC, Enu=0.5 GeV", nKEpi, KEpimin, KEpimax, ncosth, costhmin, costhmax);
  TH2D * hst_d2sig_dKEpidOmg_cohnc_0500MeV = new TH2D("hst_d2sig_dKEpidOmg_cohnc_0500MeV","d^2sig / dKEpi dOmega, COH NC, Enu=0.5 GeV", nKEpi, KEpimin, KEpimax, ncosth, costhmin, costhmax);
  TH2D * hst_d2sig_dKEpidOmg_cohcc_1000MeV = new TH2D("hst_d2sig_dKEpidOmg_cohcc_1000MeV","d^2sig / dKEpi dOmega, COH CC, Enu=1.0 GeV", nKEpi, KEpimin, KEpimax, ncosth, costhmin, costhmax);
  TH2D * hst_d2sig_dKEpidOmg_cohnc_1000MeV = new TH2D("hst_d2sig_dKEpidOmg_cohnc_1000MeV","d^2sig / dKEpi dOmega, COH NC, Enu=1.0 GeV", nKEpi, KEpimin, KEpimax, ncosth, costhmin, costhmax);
  TH2D * hst_d2sig_dKEpidOmg_cohcc_1500MeV = new TH2D("hst_d2sig_dKEpidOmg_cohcc_1500MeV","d^2sig / dKEpi dOmega, COH CC, Enu=1.5 GeV", nKEpi, KEpimin, KEpimax, ncosth, costhmin, costhmax);
  TH2D * hst_d2sig_dKEpidOmg_cohnc_1500MeV = new TH2D("hst_d2sig_dKEpidOmg_cohnc_1500MeV","d^2sig / dKEpi dOmega, COH NC, Enu=1.5 GeV", nKEpi, KEpimin, KEpimax, ncosth, costhmin, costhmax);

  //
  // fill histograms
  //
  chain->Draw("pzf/sqrt(pzf*pzf+pyf*pyf+pxf*pxf):(Ef-0.13957)>>hst_d2sig_dKEpidOmg_cohcc_0500MeV","coh&&cc&&Ev>0.49&&Ev<0.51&&pdgf==211","GOFF");
  chain->Draw("pzf/sqrt(pzf*pzf+pyf*pyf+pxf*pxf):(Ef-0.13498)>>hst_d2sig_dKEpidOmg_cohnc_0500MeV","coh&&nc&&Ev>0.49&&Ev<0.51&&pdgf==111","GOFF");
  chain->Draw("pzf/sqrt(pzf*pzf+pyf*pyf+pxf*pxf):(Ef-0.13957)>>hst_d2sig_dKEpidOmg_cohcc_1000MeV","coh&&cc&&Ev>0.99&&Ev<1.01&&pdgf==211","GOFF");
  chain->Draw("pzf/sqrt(pzf*pzf+pyf*pyf+pxf*pxf):(Ef-0.13498)>>hst_d2sig_dKEpidOmg_cohnc_1000MeV","coh&&nc&&Ev>0.99&&Ev<1.01&&pdgf==111","GOFF");
  chain->Draw("pzf/sqrt(pzf*pzf+pyf*pyf+pxf*pxf):(Ef-0.13957)>>hst_d2sig_dKEpidOmg_cohcc_1500MeV","coh&&cc&&Ev>1.49&&Ev<1.51&&pdgf==211","GOFF");
  chain->Draw("pzf/sqrt(pzf*pzf+pyf*pyf+pxf*pxf):(Ef-0.13498)>>hst_d2sig_dKEpidOmg_cohnc_1500MeV","coh&&nc&&Ev>1.49&&Ev<1.51&&pdgf==111","GOFF");

  //
  // normalize
  //
  double norm_cohcc_0500MeV = hst_d2sig_dKEpidOmg_cohcc_0500MeV -> Integral("width") * 2*TMath::Pi() / sig_cohcc_0500MeV;
  double norm_cohnc_0500MeV = hst_d2sig_dKEpidOmg_cohnc_0500MeV -> Integral("width") * 2*TMath::Pi() / sig_cohnc_0500MeV;
  double norm_cohcc_1000MeV = hst_d2sig_dKEpidOmg_cohcc_1000MeV -> Integral("width") * 2*TMath::Pi() / sig_cohcc_1000MeV;
  double norm_cohnc_1000MeV = hst_d2sig_dKEpidOmg_cohnc_1000MeV -> Integral("width") * 2*TMath::Pi() / sig_cohnc_1000MeV;
  double norm_cohcc_1500MeV = hst_d2sig_dKEpidOmg_cohcc_1500MeV -> Integral("width") * 2*TMath::Pi() / sig_cohcc_1500MeV;
  double norm_cohnc_1500MeV = hst_d2sig_dKEpidOmg_cohnc_1500MeV -> Integral("width") * 2*TMath::Pi() / sig_cohnc_1500MeV;

  if (norm_cohcc_0500MeV > 0) hst_d2sig_dKEpidOmg_cohcc_0500MeV -> Scale(1./norm_cohcc_0500MeV);
  if (norm_cohnc_0500MeV > 0) hst_d2sig_dKEpidOmg_cohnc_0500MeV -> Scale(1./norm_cohnc_0500MeV);
  if (norm_cohcc_1000MeV > 0) hst_d2sig_dKEpidOmg_cohcc_1000MeV -> Scale(1./norm_cohcc_1000MeV);
  if (norm_cohnc_1000MeV > 0) hst_d2sig_dKEpidOmg_cohnc_1000MeV -> Scale(1./norm_cohnc_1000MeV);
  if (norm_cohcc_1500MeV > 0) hst_d2sig_dKEpidOmg_cohcc_1500MeV -> Scale(1./norm_cohcc_1500MeV);
  if (norm_cohnc_1500MeV > 0) hst_d2sig_dKEpidOmg_cohnc_1500MeV -> Scale(1./norm_cohnc_1500MeV);

  for(int i = 1; i <= hst_d2sig_dKEpidOmg_cohnc_1500MeV->GetNbinsX(); i++) {
    for(int j = 1; j <= hst_d2sig_dKEpidOmg_cohnc_1500MeV->GetNbinsY(); j++) {

     double KEpi    = hst_d2sig_dKEpidOmg_cohnc_1500MeV -> GetXaxis() -> GetBinCenter(i);
     double costhpi = hst_d2sig_dKEpidOmg_cohnc_1500MeV -> GetYaxis() -> GetBinCenter(j);

     double d2sig_dKEpidOmg_cohcc_0500MeV = hst_d2sig_dKEpidOmg_cohcc_0500MeV -> GetBinContent(i,j);
     double d2sig_dKEpidOmg_cohnc_0500MeV = hst_d2sig_dKEpidOmg_cohnc_0500MeV -> GetBinContent(i,j);
     double d2sig_dKEpidOmg_cohcc_1000MeV = hst_d2sig_dKEpidOmg_cohcc_1000MeV -> GetBinContent(i,j);
     double d2sig_dKEpidOmg_cohnc_1000MeV = hst_d2sig_dKEpidOmg_cohnc_1000MeV -> GetBinContent(i,j);
     double d2sig_dKEpidOmg_cohcc_1500MeV = hst_d2sig_dKEpidOmg_cohcc_1500MeV -> GetBinContent(i,j);
     double d2sig_dKEpidOmg_cohnc_1500MeV = hst_d2sig_dKEpidOmg_cohnc_1500MeV -> GetBinContent(i,j);

     d2sig_dKEpidOmg_cohcc_0500MeV = TMath::Max(0., d2sig_dKEpidOmg_cohcc_0500MeV);
     d2sig_dKEpidOmg_cohnc_0500MeV = TMath::Max(0., d2sig_dKEpidOmg_cohnc_0500MeV);
     d2sig_dKEpidOmg_cohcc_1000MeV = TMath::Max(0., d2sig_dKEpidOmg_cohcc_1000MeV);
     d2sig_dKEpidOmg_cohnc_1000MeV = TMath::Max(0., d2sig_dKEpidOmg_cohnc_1000MeV);
     d2sig_dKEpidOmg_cohcc_1500MeV = TMath::Max(0., d2sig_dKEpidOmg_cohcc_1500MeV);
     d2sig_dKEpidOmg_cohnc_1500MeV = TMath::Max(0., d2sig_dKEpidOmg_cohnc_1500MeV);

     out_stream << setw(15) << KEpi 
                << setw(15) << costhpi
                << setw(15) << d2sig_dKEpidOmg_cohcc_0500MeV
                << setw(15) << d2sig_dKEpidOmg_cohnc_0500MeV
                << setw(15) << d2sig_dKEpidOmg_cohcc_1000MeV
                << setw(15) << d2sig_dKEpidOmg_cohnc_1000MeV
                << setw(15) << d2sig_dKEpidOmg_cohcc_1500MeV
                << setw(15) << d2sig_dKEpidOmg_cohnc_1500MeV
                << endl;
   }//y
  }//x

  out_stream.close();

  // visual inspection

  //TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  //c1->Update();
}
