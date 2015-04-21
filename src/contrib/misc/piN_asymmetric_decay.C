//
// test asymmetric decays of low multiplicity (2-body pi+N) hadronic systems
// using data (spherical harmonic expansion) from Nucl.Phys.B343(1990) 285-309
//
// To run:
// root[0] .L piN_asymmetric_decay.C++
// root[1] plot_input_data();  //or
// root[1] run_evgen_test(W);
//
// Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
// University of Liverpool & STFC Rutherford Appleton Lab
//

//
// Refresher on spherical harmonics:
//
// Y^{m}_{l} := Yml = N * exp(i*m*phi) * Pml(costh)
//
// Y01 =  (1/2) * sqrt(3/pi)   * costh
// Y02 =  (1/4) * sqrt(5/pi)   * (3*cos^2th - 1)
// Y11 = -(1/2) * sqrt(3/2pi)  * sinth           * exp(i*phi)
// Y12 = -(1/2) * sqrt(15/2pi) * sinth*costh     * exp(i*phi)
// Y22 =  (1/4) * sqrt(15/2pi) * sin^2th         * exp(i*2*phi)
//

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TSpline.h>
#include <TRandom3.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

// p pi+
const int kNWBins = 6;
double Wmin        [kNWBins] = {  1.100,  1.180,  1.230,  1.280,  1.500,  1.700 };
double Wmax        [kNWBins] = {  1.180,  1.230,  1.280,  1.500,  1.700,  2.000 };
double Wc          [kNWBins] = {  1.140,  1.205,  1.255,  1.390,  1.600,  1.850 };

double CoeffReY01  [kNWBins] = { -0.22, -0.01,  0.29,  0.32,  0.61,  0.84 };
double CoeffReY02  [kNWBins] = { -0.06, -0.13, -0.19, -0.08,  0.49,  0.80 };
double CoeffReY11  [kNWBins] = {  0.21, -0.02, -0.20, -0.05,  0.12, -0.03 };
double CoeffImY11  [kNWBins] = { -0.09, -0.06, -0.06,  0.01, -0.08, -0.17 };
double CoeffReY12  [kNWBins] = {  0.26,  0.05,  0.02,  0.00,  0.09, -0.05 };
double CoeffImY12  [kNWBins] = {  0.07,  0.00,  0.16, -0.01, -0.04, -0.03 };
double CoeffReY22  [kNWBins] = { -0.04,  0.08,  0.08,  0.04,  0.09, -0.08 };
double CoeffImY22  [kNWBins] = { -0.07, -0.07, -0.01,  0.06, -0.15, -0.03 };

double dCoeffReY01 [kNWBins] = {  0.11,  0.07,  0.09,  0.09,  0.13,  0.13 };
double dCoeffReY02 [kNWBins] = {  0.12,  0.07,  0.10,  0.09,  0.15,  0.14 };
double dCoeffReY11 [kNWBins] = {  0.09,  0.06,  0.08,  0.06,  0.08,  0.07 };
double dCoeffImY11 [kNWBins] = {  0.08,  0.06,  0.08,  0.06,  0.09,  0.09 };
double dCoeffReY12 [kNWBins] = {  0.08,  0.06,  0.08,  0.06,  0.07,  0.06 };
double dCoeffImY12 [kNWBins] = {  0.07,  0.05,  0.08,  0.06,  0.08,  0.07 };
double dCoeffReY22 [kNWBins] = {  0.07,  0.05,  0.07,  0.06,  0.09,  0.09 };
double dCoeffImY22 [kNWBins] = {  0.08,  0.05,  0.08,  0.06,  0.08,  0.07 };

// cubic splines (built from data above)
//
TSpline3 splCoeffReY01("splCoeffReY01",Wc,CoeffReY01,kNWBins,"",0,0);
TSpline3 splCoeffReY02("splCoeffReY02",Wc,CoeffReY02,kNWBins,"",1,1);
TSpline3 splCoeffReY11("splCoeffReY11",Wc,CoeffReY11,kNWBins,"",1,1);
TSpline3 splCoeffImY11("splCoeffImY11",Wc,CoeffImY11,kNWBins,"",1,1);
TSpline3 splCoeffReY12("splCoeffReY12",Wc,CoeffReY12,kNWBins,"",1,1);
TSpline3 splCoeffImY12("splCoeffImY12",Wc,CoeffImY12,kNWBins,"",1,1);
TSpline3 splCoeffReY22("splCoeffReY22",Wc,CoeffReY22,kNWBins,"",1,1);
TSpline3 splCoeffImY22("splCoeffImY22",Wc,CoeffImY22,kNWBins,"",1,1);

// consts
//
const double kPi = 3.14;

// prototypes
//
double angular_distribution_pdf(double W, double th, double phi);
double coeff (double W, int m, int l, bool real);
double Y     (double th, double phi, int m, int l, bool real);
double Y00   (double th, double phi, bool real);
double Y01   (double th, double phi, bool real);
double Y02   (double th, double phi, bool real);
double Y11   (double th, double phi, bool real);
double Y12   (double th, double phi, bool real);
double Y22   (double th, double phi, bool real);

//_______________________________________________________________________________
void run_evgen_test(double W)
{
// generate events with both the old & new scheme for comparison
//
  TFile f("./events.root","recreate");
  TNtuple nt("nt","","m:i:th:phi");

  TRandom rndgen(18928928);

  const int nev = 10000000;

  double th_min  = 0;
  double th_max  = kPi;
  double phi_min = 0;
  double phi_max = 2*kPi;

// scan for max
//
  const  int nth  = 180;
  const  int nphi = 360;
  double dth  = (th_max  - th_min ) / (nth -1);
  double dphi = (phi_max - phi_min) / (nphi-1);

  double pdf_max = -1;
  for(int i=0; i<nth; i++)
      for(int j=0; j<nphi; j++)
        pdf_max = TMath::Max(pdf_max, 
              angular_distribution_pdf(W,th_min+i*dth,phi_min+j*dphi));  
  pdf_max *= 1.2;

// generate N theta,phi pairs using the new method
//
  for(int iev=0; iev<nev; iev++) {
    bool generated=false;
    while(!generated) {
      double th  = rndgen.Rndm() * kPi; 
      double phi = rndgen.Rndm() * 2*kPi;
      double t   = rndgen.Rndm() * pdf_max;
      double pdf = angular_distribution_pdf(W,th,phi);
      generated = (t<pdf);
      if(generated) { nt.Fill(1,iev,th,phi); }
    }
  }

// generate N theta,phi pairs using the old method (phase space decay)
//
   const int    Np   = 2;
   const double MN   = 0.938;
   const double Mpi  = 0.139;
   TLorentzVector p4(0.0, 0.0, 0.0, W);
   Double_t masses[Np] = { MN, Mpi } ;

   TGenPhaseSpace phspg;
   Bool_t is_allowed = phspg.SetDecay(p4, Np, masses);
   assert(is_allowed);

   double max_wght = -1;
   for (Int_t iev=0; iev<200; iev++) {
      double weight = phspg.Generate();
      max_wght = TMath::Max(max_wght, weight);
   }
   max_wght *= 1.2;

  for(int iev=0; iev<nev; iev++) {
    bool generated=false;
    while(!generated) {
      double weight = phspg.Generate();
      double t = rndgen.Rndm() * max_wght;
      generated = (t<weight);
      if(generated) {
         TLorentzVector * p4pi = phspg.GetDecay(1);
         double th  = p4pi->Theta();
         double phi = p4pi->Phi() + kPi;
         nt.Fill(-1,iev,th,phi); 
      }
    }
  }

  nt.Write();
  f.Close();
}
//_______________________________________________________________________________
void plot_input_data(void)
{
// save all distributions that are input to the new 2-body decay scheme
//
  TFile f("./inputs.root","recreate");
  TNtuple nt("nt","","th:phi:W:sum:r01:r02:r11:i11:r12:i12:r22:i22");

  const  int nw   = 11;
  const  int nth  = 180;
  const  int nphi = 360;

  double th_min   = 0;
  double th_max   = kPi;
  double phi_min  = 0;
  double phi_max  = 2*kPi;
  double dth      = (th_max  - th_min ) / (nth -1);
  double dphi     = (phi_max - phi_min) / (nphi-1);

  double W[nw] = { 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 };

  for(int i=0; i<nth; i++) {
      double th = th_min + i*dth;
      for(int j=0; j<nphi; j++) {
         double phi = phi_min + j*dphi;
         for(int k=0; k<nw; k++) {
           double w = W[k];

           double sum = angular_distribution_pdf(w, th, phi);
           double r01 = coeff(w, 0,1,true ) * Y(th,phi, 0,1,true );
           double r02 = coeff(w, 0,2,true ) * Y(th,phi, 0,2,true );
           double r11 = coeff(w, 1,1,true ) * Y(th,phi, 1,1,true );
           double i11 = coeff(w, 1,1,false) * Y(th,phi, 1,1,false);
           double r12 = coeff(w, 1,2,true ) * Y(th,phi, 1,2,true );
           double i12 = coeff(w, 1,2,false) * Y(th,phi, 1,2,false);
           double r22 = coeff(w, 2,2,true ) * Y(th,phi, 2,2,true );
           double i22 = coeff(w, 2,2,false) * Y(th,phi, 2,2,false);

           nt.Fill(th,phi,w,sum,r01,r02,r11,i11,r12,i12,r22,i22);
         }//w
      }//phi
  }//theta

  double wmin = 1.1;
  double wmax = 2.0;
  double nwb = 100;
  TH1D cr01("cr01","coefficients of R{Y01} vs W",nwb,wmin,wmax);
  TH1D cr02("cr02","coefficients of R{Y02} vs W",nwb,wmin,wmax);
  TH1D cr11("cr11","coefficients of R{Y11} vs W",nwb,wmin,wmax);
  TH1D ci11("ci11","coefficients of I{Y11} vs W",nwb,wmin,wmax);
  TH1D cr12("cr12","coefficients of R{Y12} vs W",nwb,wmin,wmax);
  TH1D ci12("ci12","coefficients of I{Y12} vs W",nwb,wmin,wmax);
  TH1D cr22("cr22","coefficients of R{Y22} vs W",nwb,wmin,wmax);
  TH1D ci22("ci22","coefficients of I{Y22} vs W",nwb,wmin,wmax);

  for(int i=1; i<cr01.GetNbinsX(); i++) {
    double w = cr01.GetBinCenter(i);
    cr01.Fill( w, coeff(w,0,1,true)  );
    cr02.Fill( w, coeff(w,0,2,true)  );
    cr11.Fill( w, coeff(w,1,1,true)  );
    ci11.Fill( w, coeff(w,1,1,false) );
    cr12.Fill( w, coeff(w,1,2,true)  );
    ci12.Fill( w, coeff(w,1,2,false) );
    cr22.Fill( w, coeff(w,2,2,true)  );
    ci22.Fill( w, coeff(w,2,2,false) );
  }

  double dW[kNWBins];
  for(int i=0; i<kNWBins; i++) {dW[i] = Wmax[i]-Wc[i];}

  TGraphAsymmErrors dr01(kNWBins,Wc,CoeffReY01,dW,dW,dCoeffReY01,dCoeffReY01);
  TGraphAsymmErrors dr02(kNWBins,Wc,CoeffReY02,dW,dW,dCoeffReY02,dCoeffReY02);
  TGraphAsymmErrors dr11(kNWBins,Wc,CoeffReY11,dW,dW,dCoeffReY11,dCoeffReY11);
  TGraphAsymmErrors di11(kNWBins,Wc,CoeffImY11,dW,dW,dCoeffImY11,dCoeffImY11);
  TGraphAsymmErrors dr12(kNWBins,Wc,CoeffReY12,dW,dW,dCoeffReY12,dCoeffReY12);
  TGraphAsymmErrors di12(kNWBins,Wc,CoeffImY12,dW,dW,dCoeffImY12,dCoeffImY12);
  TGraphAsymmErrors dr22(kNWBins,Wc,CoeffReY22,dW,dW,dCoeffReY22,dCoeffReY22);
  TGraphAsymmErrors di22(kNWBins,Wc,CoeffImY22,dW,dW,dCoeffImY22,dCoeffImY22);

  dr01.Write("dr01");
  dr02.Write("dr02");
  dr11.Write("dr11");
  di11.Write("di11");
  dr12.Write("dr12");
  di12.Write("di12");
  dr22.Write("dr22");
  di22.Write("di22");
  cr01.Write();
  cr02.Write();
  cr11.Write();
  ci11.Write();
  cr12.Write();
  ci12.Write();
  cr22.Write();
  ci22.Write();
  nt.Write();
  f.Close();
}
//_______________________________________________________________________________
//_______________________________________________________________________________
// Code implementing the new 2-body decays
//_______________________________________________________________________________
//_______________________________________________________________________________
double angular_distribution_pdf(double W, double th, double phi)
{
// return Sum{m,l} { (interpolated coeff) x (spherical harmonic)}
//
   double fr = 1 +
               coeff(W, 0,1,true ) * Y(th,phi, 0,1,true ) +
               coeff(W, 0,2,true ) * Y(th,phi, 0,2,true ) +
               coeff(W, 1,1,true ) * Y(th,phi, 1,1,true ) +
               coeff(W, 1,2,true ) * Y(th,phi, 1,2,true ) +
               coeff(W, 2,2,true ) * Y(th,phi, 2,2,true );

   double fi = coeff(W, 1,1,false) * Y(th,phi, 1,1,false) +
               coeff(W, 1,2,false) * Y(th,phi, 1,2,false) +
               coeff(W, 2,2,false) * Y(th,phi, 2,2,false);

   double f = TMath::Sqrt(fr*fr+fi*fi);
   return f;
}
//_______________________________________________________________________________
double coeff(double W, int m, int l, bool real)
{
  if      (m==0 && l==0 &&  real) return 1;
  else if (m==0 && l==1 &&  real) return splCoeffReY01.Eval(W);
  else if (m==0 && l==2 &&  real) return splCoeffReY02.Eval(W);
  else if (m==1 && l==1 &&  real) return splCoeffReY11.Eval(W);
  else if (m==1 && l==1 && !real) return splCoeffImY11.Eval(W);
  else if (m==1 && l==2 &&  real) return splCoeffReY12.Eval(W);
  else if (m==1 && l==2 && !real) return splCoeffImY12.Eval(W);
  else if (m==2 && l==2 &&  real) return splCoeffReY22.Eval(W);
  else if (m==2 && l==2 && !real) return splCoeffImY22.Eval(W);
  else
     return 0;  
}
//_______________________________________________________________________________
double Y(double th, double phi, int m, int l, bool real)
{
  if      (m==0 && l==0) return Y00(th,phi,real);
  else if (m==0 && l==1) return Y01(th,phi,real);
  else if (m==0 && l==2) return Y02(th,phi,real);
  else if (m==1 && l==1) return Y11(th,phi,real);
  else if (m==1 && l==2) return Y12(th,phi,real);
  else if (m==2 && l==2) return Y22(th,phi,real);
  else
     return 0;
}
//_______________________________________________________________________________
double Y00(double /*th*/, double /*phi*/, bool real)
{
  if(!real) return 0;
  return 0.5/TMath::Sqrt(kPi);
}
//_______________________________________________________________________________
double Y01(double th, double /*phi*/, bool real)
{
  if(!real) return 0;
  return 0.5 * TMath::Sqrt(3./kPi) * TMath::Cos(th);
}
//_______________________________________________________________________________
double Y02(double th, double /*phi*/, bool real)
{
  if(!real) return 0;
  return 0.25 * TMath::Sqrt(5./kPi) * (3.* TMath::Power(TMath::Cos(th),2.) - 1);
}
//_______________________________________________________________________________
double Y11(double th, double phi, bool real)
{
  double a = -0.5 * TMath::Sqrt(1.5/kPi) * TMath::Sin(th);
  if(!real) return a * TMath::Sin(phi);
  else      return a * TMath::Cos(phi);
}
//_______________________________________________________________________________
double Y12(double th, double phi, bool real)
{
  double a = -0.5 * TMath::Sqrt(7.5/kPi) * TMath::Cos(th) * TMath::Sin(th);
  if(!real) return a * TMath::Sin(phi);
  else      return a * TMath::Cos(phi);
}
//_______________________________________________________________________________
double Y22(double th, double phi, bool real)
{
  double a = 0.25 * TMath::Sqrt(7.5/kPi) * TMath::Power( TMath::Sin(th), 2.);
  if(!real) return a * TMath::Sin(2*phi);
  else      return a * TMath::Cos(2*phi);
}
//_______________________________________________________________________________
