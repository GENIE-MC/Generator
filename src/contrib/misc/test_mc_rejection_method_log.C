//
// MC rejection method: Sanity check
//
// Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
// University of Liverpool & STFC Rutherford Appleton Lab
//

#include <iostream>

#include <TRandom.h>
#include <TF1.h>
#include <TH1D.h>

void test_mc_rejection_method_log()
{
  const int    N        = 100000;
  const int    nbins    = 200;
  const double ymax     = 1.1;
  const double xmin     = 1E-2;
  const double xmax     = 10;
  const double logxmin  = TMath::Log(xmin);
  const double logxmax  = TMath::Log(xmax);
  const double dlogx    = logxmax-logxmin;

  TRandom rg;

  TF1 *  func = new TF1  ("func","1/(x+1)",0,10);
  TH1D * hgen = new TH1D ("hgen","generated",nbins,xmin,xmax);

  for(int i=0; i<N; i++) {
     cout << "..................." << i << endl;
     bool selected=false;
     while(1) {
       double xg = TMath::Exp(logxmin + rg.Uniform() * dlogx);
       double yg = ymax * rg.Uniform();
       double yc = xg * func->Eval(xg);
       selected = (yg<yc);
       if(selected) {
           hgen->Fill(xg);
           break;
       }
     }
  }

  double IF = func->Integral(xmin,xmax);
  double IH = hgen->Integral("width");
  double sc = IF/IH;
  hgen->Scale(sc);

  hgen->Draw();
  func->Draw("same");
}
