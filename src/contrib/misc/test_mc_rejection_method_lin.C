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

void test_mc_rejection_method_lin()
{
  const int    N     = 500000;
  const int    nbins = 300;
  const double ymax  = 1.1;
  const double xmin  = 0;
  const double xmax  = 10;
  const double dx    = xmax-xmin;

  TRandom rg;

  TF1 *  func = new TF1  ("func","1/(x+1)",0,10);
  TH1D * hgen = new TH1D ("hgen","generated",nbins,xmin,xmax);

  for(int i=0; i<N; i++) {
     cout << "..................." << i << endl;
     bool selected=false;
     while(1) {
       double xg = xmin + rg.Uniform() * dx;
       double yg = ymax * rg.Uniform();
       double yc = func->Eval(xg);
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
