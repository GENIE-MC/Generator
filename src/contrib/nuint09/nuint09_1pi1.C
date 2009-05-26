//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// 1PI.1:
// Total cross section for 1 pion (pi+ only) production as a function of neutrino energy
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
const int kNRunsPerCase   = 10;
const int kNEvtPerRun     = 100000;

const char * kLabel[kNSamples] = 
{  
 // available samples
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};
const int kRunNu1PI1[kNSamples][kNWCur][kNRunsPerCase] =
{
 /* indices : sample ; cc/nc */
 {
  /* 0,0 (nu_mu C12,  CC) */ {909900, 909901, 909902, 909903, 909904, 909905, 909906, 909907, 909908, 909909}
 },
 {
  /* 1,0 (nu_mu O16,  CC) */ {919900, 919901, 919902, 919903, 919904, 919905, 919906, 919907, 919908, 919909}
 },
 {
  /* 2,0 (nu_mu Fe56, CC) */ {929900, 929901, 929902, 929903, 929904, 929905, 929906, 929907, 929908, 929909}
 }
};
int kA[kNSamples] = 
{  
 // A for nuclear target at each sample
 /* 0 */  12,
 /* 1 */  16,
 /* 2 */  56
};

void nuint09_1pi1(int isample, int single_pion_sources=0, int stage=1)
{
  cout << " ***** running: 1PI.1" << endl;

  if(isample<0 || isample >= kNSamples) return;
  if(single_pion_sources<0 || single_pion_sources>2) return;

  const char * label = kLabel[isample];

  int A = kA[isample];

  // get cross section graphs

  TFile fsig("../sig/splines.root","read");
  TDirectory * sig_dir = (TDirectory *) fsig.Get(label);  

  TGraph * sig_graph_totcc = (TGraph*) sig_dir->Get("tot_cc");

  // range & spacing

  const int    nEnu   = 60;
  const double Enumin =  0.05;
  const double Enumax = 30.00;

  // create output stream

  ostringstream out_filename;
  out_filename << label;

  if      (single_pion_sources==0) out_filename << ".1pi_1a.";
  else if (single_pion_sources==1) out_filename << ".1pi_1b.";
  else if (single_pion_sources==2) out_filename << ".1pi_1c.";

  if(stage==0) out_filename << "no_FSI.";

  out_filename << "sig1pi_vs_Enu.data";
  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [1PI.1]:" << endl;
  out_stream << "#  Total cross section for 1 pion (pi+ only) production as a function of energy" << endl;
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
  out_stream << "#   - neutrino energy E in GeV, linear spacing between Emin = " << Enumin << " GeV, Emax = " << Enumax << " GeV "  << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 " << endl;
  out_stream << "#   - quoted cross section is nuclear cross section divided with number of nucleons A" << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  Energy     |   sig(nu_mu + A -> mu- 1pi+ X)   |  "  << endl;

  out_stream << setiosflags(ios::fixed) << setprecision(6);

  //
  // load event data
  // 
    
  TChain * chain = new TChain("gst");
  
  // loop over CC/NC cases
  for(int iwkcur=0; iwkcur<kNWCur; iwkcur++) {
       // loop over runs for current case
       for(int ir=0; ir<kNRunsPerCase; ir++) {
          // build filename
          ostringstream filename;
          int run_number = kRunNu1PI1[isample][iwkcur][ir];

          cout << "isample = " << isample << ", iwkcur = " << iwkcur << ", ir = " << ir << ", run = " << run_number << endl;

          filename << "../gst/gntp." << run_number << ".gst.root";
          // add to chain
          cout << "Adding " << filename.str() << " to event chain" << endl;
          chain->Add(filename.str().c_str());
       }
  } 

  //
  // book histograms
  //
  TH1D * hst_cc1pip = new TH1D("hst_cc1pip","CC1pi+ events, Enu spectrum", nEnu, Enumin, Enumax);
  TH1D * hst_allcc  = new TH1D("hst_allcc", "all CC events, Enu spectrum", nEnu, Enumin, Enumax);

  //
  // fill histograms
  //

  chain->Draw("Ev>>hst_allcc", "cc","GOFF");

  if(stage==1) {
    if(single_pion_sources==0) {
      //all sources
      chain->Draw("Ev>>hst_cc1pip","cc&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
    }
    else if(single_pion_sources==1) {
      //P33(1232) only
      chain->Draw("Ev>>hst_cc1pip","cc&&resid==0&&res&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
    }
    else if(single_pion_sources==2) {
      //all resonances only
      chain->Draw("Ev>>hst_cc1pip","cc&&res&&pdgf==211&&nfpip==1&&nfpim==0&&nfpi0==0","GOFF");
    }
  }
  else if(stage==0) {
    if(single_pion_sources==0) {
      //all sources
      chain->Draw("Ev>>hst_cc1pip","cc&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
    }
    else if(single_pion_sources==1) {
      //P33(1232) only
      chain->Draw("Ev>>hst_cc1pip","cc&&resid==0&&res&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
    }
    else if(single_pion_sources==2) {
      //all resonances only
      chain->Draw("Ev>>hst_cc1pip","cc&&res&&pdgi==211&&nipip==1&&nipim==0&&nipi0==0","GOFF");
    }
  }

  cout << "CC1pi+ nevents: " << hst_cc1pip->GetEntries() << endl;
  cout << "ALLCC  nevents: " << hst_allcc->GetEntries() << endl;

  //
  // compute sig(CC1pi+) and write out in txt file expected by NuINT organizers
  // Also write out root graph in temp file to be re-used at later benchmarking calc
  //

  double e[nEnu], sig[nEnu];

  for(int i=1; i <= hst_cc1pip->GetNbinsX(); i++) {

     double Enu = hst_cc1pip->GetBinCenter(i);

     double Ncc1pip = hst_cc1pip -> GetBinContent(i);
     double Nallcc  = hst_allcc  -> GetBinContent(i);

     double sig_cc1pip=0;
     if(Nallcc>0) { sig_cc1pip = (Ncc1pip/Nallcc) * sig_graph_totcc->Eval(Enu) / A; }

     out_stream << setw(15) << Enu << setw(15) << sig_cc1pip << endl;

     e  [i-1] = Enu;
     sig[i-1] = sig_cc1pip;
  }

  out_stream.close();

  TFile ftmp("./cc1pip_tmp.root","recreate");
  TGraph * grCC1pip = new TGraph(nEnu,e,sig);  
  grCC1pip->Write("CC1pip");
  hst_allcc->Write();
  hst_cc1pip->Write();

  // visual inspection
  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  grCC1pip->Draw("alp");
  c1->Update();

  ftmp.Close();

  delete chain;
}
