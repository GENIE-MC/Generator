//____________________________________________________________________________
/*!

\class    genie::mc_vs_data::NuXSecFunc

\brief    A set of functors used in cross-section validation/tuning programs.
          Each functor contains a recipe for extracting a cross-secion out
          of ...

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Apr 17, 2012

\cpright  Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _NU_XSEC_FUNCTIONS_H_
#define _NU_XSEC_FUNCTIONS_H_

#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>

using std::string;

namespace genie {
namespace mc_vs_data {

// Base functor class
//
class NuXSecFunc
{
public:
  virtual ~NuXSecFunc() {}

  virtual TGraph * operator()
    (TFile *  /* genie_xsec_file */, 
     TChain * /* genie_event_tree */, 
     double   /* Emin */, 
     double   /* Emax */, 
     int      /* n */, 
     bool     /* scale_with_E */) 
    { 
      return 0; 
    }

protected:
  NuXSecFunc() {}
};

//............................................................................
// Extracting cross-section for an exclusive f/s from an event sample as
// the fraction of events with that f/s times the inclusive cross-section.
// This is a powerfull method that allows one to calculate any cross-section
// for any f/s state as a function of any kinematical parameter.
//
class NuXSecFromEventSample: public NuXSecFunc
{
public:
  NuXSecFromEventSample(
    string xsec_dir, string incl_xsec_spline, string incl_selection, string selection);
 ~NuXSecFromEventSample();

  TGraph * operator() (TFile * genie_xsec_file, TChain * genie_event_tree, double Emin, double Emax, int n, bool scale_with_E);

private:
  string fXSecDir;          
  string fInclXSecSpline;   
  string fInclEvtSelection; 
  string fEvtSelection;     
};

//
//............................................................................
// Get cross-section directly from the input cross-section spline file.
// This is a straightforward method but there is only a limited number of
// cases where it is applicable. For almost all exclusive inelastic reactions
// one needs to use the nuXSecFromEventSample functor.
//
class NuXSecDirectlyFromXSecFile: public NuXSecFunc
{
public:
  NuXSecDirectlyFromXSecFile(string xsec_dir, string xsec_spline);
 ~NuXSecDirectlyFromXSecFile();

  TGraph * operator() (TFile * genie_xsec_file, TChain * /*genie_event_tree*/, double Emin, double Emax, int n, bool scale_with_E);

private:
  string fXSecDir;          
  string fXSecSpline;   
};

//
//............................................................................
// Get cross-section by combining two cross-section splines taken directly 
// from the input cross-section file.
// This is usefull when, for example, we need to calculate the cross-section
// for an isoscalar target by averaging the vp and vn cross-sections
//
class NuXSecCombineSplinesFromXSecFile: public NuXSecFunc
{
public:
  NuXSecCombineSplinesFromXSecFile(
      double factor_1, string xsec_dir_1, string xsec_spline_1,
      double factor_2, string xsec_dir_2, string xsec_spline_2);
 ~NuXSecCombineSplinesFromXSecFile();

  TGraph * operator() (TFile * genie_xsec_file, TChain * /*genie_event_tree*/, double Emin, double Emax, int n, bool scale_with_E);

private:
  double fFactor1;
  string fXSecDir1;          
  string fXSecSpline1;   
  double fFactor2;
  string fXSecDir2;          
  string fXSecSpline2;   
};
//____________________________________________________________________________

} // mc_vs_data namespace
} // genie namepace

#endif  // _NU_XSEC_FUNCTIONS_H_
