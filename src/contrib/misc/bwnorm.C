//
// Test Breit-Wigner normalization
//
// To run:
//    root[0] .L bwnorm.C++
//    root[1] bwnorm()
//
// Costas Andreopoulos
//

#include <iostream>

#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>

//
// constants
//

const int kNRes = 18;

const char * kResName[kNRes] =  {
  "P33(1232)", "S11(1535)", "D13(1520)", "S11(1650)", "D13(1700)", 
  "D15(1675)", "S31(1620)", "D33(1700)", "P11(1440)", "P33(1600)",
  "P13(1720)", "F15(1680)", "P31(1910)", "P33(1920)", "F35(1905)",
  "F37(1950)", "P11(1710)", "F17(1970)"
};

const int kResN[kNRes] = {
  0, 1, 1, 1, 1, 
  1, 1, 1, 1, 9,
  2, 2, 2, 2, 2,
  2, 2, 0
};

const int kResL[kNRes] = {
  1, 0, 2, 0, 2, 
  2, 0, 2, 1, 1,
  1, 3, 1, 1, 3,
  3, 1, 0
};

const double kResMass[kNRes] = {
  1.232, 1.535, 1.520, 1.650, 1.700, 
  1.675, 1.620, 1.700, 1.440, 1.600,
  1.720, 1.680, 1.910, 1.920, 1.905,
  1.950, 1.710, 1.970
};

const double kResWidth[kNRes] = {
  0.120, 0.150, 0.120, 0.150, 0.100, 
  0.150, 0.150, 0.300, 0.350, 0.350,
  0.150, 0.130, 0.250, 0.200, 0.350,
  0.300, 0.100, 0.325
};

const int kResColor[kNRes] = {
  1, 8, 2, 8, 2, 
  2, 8, 2, 1, 1,
  1, 5, 1, 1, 5,
  5, 1, 8
};

const int kResStyle[kNRes] = {
  1, 1, 1, 2, 2, 
  3, 3, 4, 2, 3,
  4, 1, 2, 2, 2,
  3, 3, 4
};

//
// func prototypes
//

void   bwnorm             (void);
void   print_norm         (int method, double Wmax = 3.0);
double bwfunc             (double mR, double gRo, int L, double W);
double bwintegrate        (double mR, double gRo, int L, double Wmin, double Wmax);
double bwintegrate_neugen (double mR, double gRo, int L, int n);


void bwnorm(void)
{
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  c->SetFillColor(0);
  c->SetBorderMode(0);

  TLegend * lg = new TLegend(0.15, 0.30, 0.35, 0.85);
  lg->SetFillColor(0);
  
  const int kNW = 150;

  double norm[kNRes][kNW];
  double wmax[kNRes][kNW];

  TGraph * grnorm[kNRes];

  for(int ires=0; ires<kNRes; ires++) {

    double m = kResMass  [ires];
    double w = kResWidth [ires];
    int    l = kResL     [ires];

    double Wmin = 0.01;
    double Wmax = m+w;
    for(int iw=0; iw<kNW; iw++) {
      Wmax+=w;
      wmax[ires][iw] = Wmax;
      norm[ires][iw] = bwintegrate(m,w,l,Wmin,Wmax);
    }//w

    grnorm[ires] = new TGraph(kNW, wmax[ires], norm[ires]);

    grnorm[ires]->SetLineWidth(2);
    grnorm[ires]->SetLineStyle(kResStyle[ires]);
    grnorm[ires]->SetLineColor(kResColor[ires]);

    lg->AddEntry(grnorm[ires], kResName[ires], "L");

  }//res

  c->cd();
  TH1F * hframe = (TH1F*)c->DrawFrame(1.0, 0.6, 3.0, 1.2);
  hframe->GetXaxis()->SetTitle("W_{max} (GeV) in integration");
  hframe->GetYaxis()->SetTitle("normalization");
  hframe->Draw();
  for(int ires=0; ires<kNRes; ires++) {
    grnorm[ires]->Draw("L");
  }
  lg->Draw();
  c->Update();
}
//__________________________________________________________________________
void print_norm(int method, double Wmax)
{
  for(int ires=0; ires<kNRes; ires++) {

    double m = kResMass  [ires];
    double w = kResWidth [ires];
    int    l = kResL     [ires];
    int    n = kResN     [ires];

    double norm=0;
    if      (method==0) norm = bwintegrate_neugen(m,w,l,n);
    else if (method==1) norm = bwintegrate(m,w,l,0.01,Wmax);

    cout << kResName[ires] << " --> Norm = " << norm 
         << " (multiplied by pi -> " 
         << norm * TMath::Pi() << ")" << endl;
  }
}
//__________________________________________________________________________
double bwintegrate(
   double mR, double gRo, int L, double Wmin, double Wmax)
{
// integrate within input W range
//
  int N = 1000* TMath::Nint( (Wmax-Wmin)/gRo );
  if(N%2==0) N++;

  double dW = (Wmax-Wmin)/(N-1);

  double sum = 0.5 * (bwfunc(mR,gRo,L,Wmin) + bwfunc(mR,gRo,L,Wmax));

  for(int i=1; i<N-1; i++) {
    double W = Wmin + i*dW;
    sum += ( bwfunc(mR,gRo,L,W) * (i%2+1) );
  }
  sum *= (2.*dW/3.);

  return sum;
}
//__________________________________________________________________________
double bwintegrate_neugen(double mR, double gRo, int L, int n)
{
// integrate using std neugen W range
//
  int NW = 4;
  if(n==2) NW=2;
  if(n==0) NW=6;

  double Wmin = 0.001;
  double Wmax = mR + NW*gRo;

  return bwintegrate(mR,gRo,L,Wmin,Wmax);
}
//__________________________________________________________________________
double bwfunc(double mR, double gRo, int L, double W)
{
// breit-wigner function
//
  double pi   = 3.141;
  double mN   = 0.938;
  double mPi  = 0.140;
  double mR2  = TMath::Power(mR, 2);
  double mN2  = TMath::Power(mN, 2);
  double mPi2 = TMath::Power(mPi,2);
  double W2   = TMath::Power(W,  2);
  double qpW2 = TMath::Power(W2  - mN2 - mPi2, 2) - 4*mN2*mPi2;
  double qpM2 = TMath::Power(mR2 - mN2 - mPi2, 2) - 4*mN2*mPi2;
  double qpW  = TMath::Sqrt(TMath::Max(0.,qpW2)) / (2*W);
  double qpM  = TMath::Sqrt(TMath::Max(0.,qpM2)) / (2*mR);
  double gR   = gRo * TMath::Power( qpW/qpM, 2*L+1 );
  double gR2  = TMath::Power(gR, 2);
  double brwg = (0.5/pi)*gR / (TMath::Power(W-mR, 2) + 0.25*gR2);
  return brwg;
}
//__________________________________________________________________________
