//
// Chain JNUBEAM SuperK flux ntuples & write out flux histogram files
//
// % root
// root[0] .L generate_sk_flux_histograms.C++
// root[1] generate_sk_flux_histograms(flux_dir, file_prefix, file_suffix, min_run, max_run, Emin, Emax, Ebinsize);
//
// C.Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
//

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

#include <TFile.h>
#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>

using namespace std;

void generate_sk_flux_histograms(

   string flux_dir    = "./", 
   string file_prefix = "nu.sk_flukain.", 
   string file_suffix = ".root", 
   int    min_run     = 0, 
   int    max_run     = 499,
   double Emin        = 0.,    /* GeV */
   double Emax        = 30.,   /* GeV */
   double Ebinsz      = 0.05   /* GeV */  )

{
  //
  // chain
  //
  TChain chain("h2000", "JNUBEAM SuperK flux ntuples");
  for(int irun = min_run; irun <= max_run; irun++) {
    ostringstream file;
    file << flux_dir << "/" << file_prefix << irun << file_suffix;
    cout << "Adding file..... : " << file.str() << endl;
    chain.AddFile(file.str().c_str());
  }

  //
  // book histograms
  //
  assert ( Emin<Emax );
  assert ( Emin>=0.  );
  assert ( Ebinsz>0. );
  int nbins = TMath::FloorNint((Emax-Emin)/Ebinsz);
  TH1D * numu_flux    = new TH1D("numu_flux",    "numu flux at SuperK",    nbins, Emin, Emax);
  TH1D * numubar_flux = new TH1D("numubar_flux", "numubar flux at SuperK", nbins, Emin, Emax);
  TH1D * nue_flux     = new TH1D("nue_flux",     "nue flux at SuperK",     nbins, Emin, Emax);
  TH1D * nuebar_flux  = new TH1D("nuebar_flux",  "nuebar flux at SuperK",  nbins, Emin, Emax);

  //
  // fill histograms
  //
  chain.Draw("Enu>>numu_flux",    "norm*(mode>=11 && mode<=19)", "goff");
  chain.Draw("Enu>>numubar_flux", "norm*(mode>=21 && mode<=29)", "goff");
  chain.Draw("Enu>>nue_flux",     "norm*(mode>=31 && mode<=39)", "goff");
  chain.Draw("Enu>>nuebar_flux",  "norm*(mode>=41 && mode<=49)", "goff");

  //
  // save histograms
  //
  TFile f("./sk_flux_histograms.root", "recreate");
  numu_flux    -> Write ("numu_flux"   );
  numubar_flux -> Write ("numubar_flux");
  nue_flux     -> Write ("nue_flux"    );
  nuebar_flux  -> Write ("nuebar_flux" );
  f.Close();
}
