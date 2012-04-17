//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TH1D.h>
#include <TMath.h>

#include "Messenger/Messenger.h"
#include "validation/NuXSec/NuXSecFunc.h"

using namespace genie;
using namespace genie::mc_vs_data;


//____________________________________________________________________________
NuXSecFromEventSample::NuXSecFromEventSample(
  string xsec_dir, string incl_xsec_spline, string incl_selection, string selection) :
NuXSecFunc(),
fXSecDir          (xsec_dir),           
fInclXSecSpline   (incl_xsec_spline),   
fInclEvtSelection (incl_selection),
fEvtSelection     (selection)
{

}
//............................................................................
NuXSecFromEventSample::~NuXSecFromEventSample() 
{ 

}
//............................................................................
TGraph * NuXSecFromEventSample::operator()
    (TFile * genie_xsec_file, TChain * genie_event_tree, 
     double Emin, double Emax, int n, bool scale_with_E)
{
  if(!genie_xsec_file ) return 0;
  if(!genie_event_tree) return 0;

  TDirectory * xsec_dir  = (TDirectory *) genie_xsec_file->Get(fXSecDir.c_str());
  if(!xsec_dir) return 0;
  TGraph * incl_xsec_spline = (TGraph*) xsec_dir->Get(fInclXSecSpline.c_str());
  if(!incl_xsec_spline) return 0;

  TH1D * h     = 0;
  TH1D * hincl = 0;
  LOG("gvldtest", pNOTICE)   << "Emin = " << Emin << ", Emax = " << Emax;

  bool inlogE = false;
  if(inlogE) {
     h     = new TH1D("h",     "", n, TMath::Log10(Emin), TMath::Log10(Emax));
     hincl = new TH1D("hincl", "", n, TMath::Log10(Emin), TMath::Log10(Emax)); 
     genie_event_tree->Draw("log10(Ev)>>h",     fEvtSelection.c_str(),     "goff");
     genie_event_tree->Draw("log10(Ev)>>hincl", fInclEvtSelection.c_str(), "goff");
  } else { 
     h     = new TH1D("h",     "", n, Emin, Emax);
     hincl = new TH1D("hincl", "", n, Emin, Emax); 
     genie_event_tree->Draw("Ev>>h", fEvtSelection.c_str(), "goff");
     LOG("gvldtest", pNOTICE)  
        << "Selection: " << fEvtSelection 
        << " retrieved " << genie_event_tree->GetSelectedRows() << " entries";
     genie_event_tree->Draw("Ev>>hincl", fInclEvtSelection.c_str(), "goff");
     LOG("gvldtest", pNOTICE)  
        << "Selection: " << fInclEvtSelection 
        << " retrieved " << genie_event_tree->GetSelectedRows() << " entries";
  }
  h->Divide(hincl);

  double * energy_array = new double [n];
  double * xsec_array   = new double [n];

  for(int i = 0; i < n; i++) {
    int ibin = i+1;
    double energy = (inlogE) ? 
       TMath::Power(10., h->GetBinCenter(ibin)) : h->GetBinCenter(ibin);
    double event_fraction = h->GetBinContent(ibin);
    double incl_xsec =  incl_xsec_spline->Eval(energy);
    double xsec = event_fraction * incl_xsec;
  
    if(scale_with_E) {
      assert(energy>0);
      xsec /= energy;
    }

    LOG("gvldtest", pNOTICE)  
        << "xsec(E=" << energy << " GeV) = " << xsec 
        << " (event fraction = " << event_fraction << ", incl xsec = " << incl_xsec << ")";

    energy_array[i] = energy;
    xsec_array[i]   = xsec;
  }
  TGraph * model = new TGraph(n,energy_array,xsec_array);
  delete h;
  delete hincl;
  delete [] energy_array;
  delete [] xsec_array;
  return model;
}
//____________________________________________________________________________
NuXSecDirectlyFromXSecFile::NuXSecDirectlyFromXSecFile(
    string xsec_dir, string xsec_spline) :
NuXSecFunc(),
fXSecDir    (xsec_dir),           
fXSecSpline (xsec_spline)
{
}
//............................................................................
NuXSecDirectlyFromXSecFile::~NuXSecDirectlyFromXSecFile() 
{ 
}
//............................................................................
TGraph * NuXSecDirectlyFromXSecFile::operator()
  (TFile * genie_xsec_file, TChain * /*genie_event_tree*/, 
   double Emin, double Emax, int n, bool scale_with_E)
{
  if(!genie_xsec_file) return 0;

  TDirectory * xsec_dir  = (TDirectory *) genie_xsec_file->Get(fXSecDir.c_str());
  if(!xsec_dir) return 0;

  TGraph * xsec_spline = (TGraph*) xsec_dir->Get(fXSecSpline.c_str());
  if(!xsec_spline) return 0;

  double * energy_array = new double[n];
  double * xsec_array   = new double[n];

  bool inlogE = true;
  for(int i = 0; i < n; i++) {
    double energy = (inlogE) ? 
      TMath::Power(10., TMath::Log10(Emin) + i * TMath::Log10(Emax-Emin)/(n-1)) : 
         Emin + i*(Emax-Emin)/(n-1);
    double xsec = xsec_spline->Eval(energy);
    if(scale_with_E) {
      assert(energy>0);
      xsec /= energy;
    }
    energy_array[i] = energy;
    xsec_array[i] = xsec;
  }

  TGraph * model = new TGraph(n,energy_array,xsec_array);
  delete [] energy_array;
  delete [] xsec_array;
  return model;
}
//____________________________________________________________________________
NuXSecCombineSplinesFromXSecFile::NuXSecCombineSplinesFromXSecFile(
   double factor_1, string xsec_dir_1, string xsec_spline_1,
   double factor_2, string xsec_dir_2, string xsec_spline_2) :
NuXSecFunc(),
fFactor1     (factor_1),
fXSecDir1    (xsec_dir_1),           
fXSecSpline1 (xsec_spline_1),
fFactor2     (factor_2),
fXSecDir2    (xsec_dir_2),           
fXSecSpline2 (xsec_spline_2)
{

}
//............................................................................
NuXSecCombineSplinesFromXSecFile::~NuXSecCombineSplinesFromXSecFile() 
{ 
  
}
//............................................................................
TGraph * NuXSecCombineSplinesFromXSecFile::operator()
  (TFile * genie_xsec_file, TChain * /*genie_event_tree*/, 
   double Emin, double Emax, int n, bool scale_with_E)
{
  if(!genie_xsec_file) return 0;

  TDirectory * xsec_dir_1  = (TDirectory *) genie_xsec_file->Get(fXSecDir1.c_str());
  if(!xsec_dir_1) return 0;
  TGraph * xsec_spline_1 = (TGraph*) xsec_dir_1->Get(fXSecSpline1.c_str());
  if(!xsec_spline_1) return 0;

  TDirectory * xsec_dir_2  = (TDirectory *) genie_xsec_file->Get(fXSecDir2.c_str());
  if(!xsec_dir_2) return 0;
  TGraph * xsec_spline_2 = (TGraph*) xsec_dir_2->Get(fXSecSpline2.c_str());
  if(!xsec_spline_2) return 0;

  double f_1 = fFactor1;
  double f_2 = fFactor2;

  double * energy_array = new double[n];
  double * xsec_array   = new double[n];

  bool inlogE = true;
  for(int i = 0; i < n; i++) {
    double energy = (inlogE) ? 
       TMath::Power(10., TMath::Log10(Emin) + i * TMath::Log10(Emax-Emin)/(n-1)) : 
       Emin + i*(Emax-Emin)/(n-1);     
    double xsec_1 = xsec_spline_1->Eval(energy);
    double xsec_2 = xsec_spline_2->Eval(energy);
    double xsec = f_1*xsec_1 + f_2*xsec_2;
    if(scale_with_E) {
      assert(energy>0);
      xsec /= energy;
    }
    energy_array[i] = energy;
    xsec_array[i] = xsec;
  }

  TGraph * model = new TGraph(n,energy_array,xsec_array);
  delete [] energy_array;
  delete [] xsec_array;
  return model;
}
//____________________________________________________________________________
