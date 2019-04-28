//________________________________________________________________________________________
/*!

\macro   make_sk_xsec_table.C

\brief   Write-out the cross-section table needed by SuperK softw for MC normalization

\usage   % genie 'make_sk_xsec_table.C("/some/path/genie_t2k_splines.xml")'

\author  Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Ryan Terri <r.terri \at qmul.ac.uk>>
         Queen Mary, University of London

\created Nov 24, 2008

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_________________________________________________________________________________________

#include <fstream>
#include <iomanip>
#include <cassert>

#include "Conventions/Units.h"

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;

using namespace genie;

void make_sk_xsec_table
(
  string genie_cross_section_xml_file,
  double EvMin  =  0.025, // minimum energy in output file, GeV
  double EvMax  = 15.000, // maximum energy in output file, GeV
  double dEv   =  0.050   // energy binning, GeV
)
{
  //
  // Load cross-section splines
  //

  XSecSplineList * xspl = XSecSplineList::Instance();
  XmlParserStatus_t ist = xspl->LoadFromXml(genie_cross_section_xml_file);
  assert(ist == kXmlOK);

  //
  // Create GENIE & configure event generation drivers for all init states
  //

  cout << "Initializing event generation drivers!" << endl;

  GEVGDriver numu_H1;
  GEVGDriver numu_O16;
  GEVGDriver numubar_H1;
  GEVGDriver numubar_O16;
  GEVGDriver nue_H1;
  GEVGDriver nue_O16;
  GEVGDriver nuebar_H1;
  GEVGDriver nuebar_O16;
  numu_H1.    Configure (kPdgNuMu,     pdg::IonPdgCodeToZ(kPdgTgtFreeP), pdg::IonPdgCodeToA(kPdgTgtFreeP));
  numu_O16.   Configure (kPdgNuMu,     pdg::IonPdgCodeToZ(kPdgTgtO16),   pdg::IonPdgCodeToA(kPdgTgtO16)  );
  numubar_H1. Configure (kPdgAntiNuMu, pdg::IonPdgCodeToZ(kPdgTgtFreeP), pdg::IonPdgCodeToA(kPdgTgtFreeP));
  numubar_O16.Configure (kPdgAntiNuMu, pdg::IonPdgCodeToZ(kPdgTgtO16),   pdg::IonPdgCodeToA(kPdgTgtO16)  );
  nue_H1.     Configure (kPdgNuE,      pdg::IonPdgCodeToZ(kPdgTgtFreeP), pdg::IonPdgCodeToA(kPdgTgtFreeP));
  nue_O16.    Configure (kPdgNuE,      pdg::IonPdgCodeToZ(kPdgTgtO16),   pdg::IonPdgCodeToA(kPdgTgtO16)  );
  nuebar_H1.  Configure (kPdgAntiNuE,  pdg::IonPdgCodeToZ(kPdgTgtFreeP), pdg::IonPdgCodeToA(kPdgTgtFreeP));
  nuebar_O16. Configure (kPdgAntiNuE,  pdg::IonPdgCodeToZ(kPdgTgtO16),   pdg::IonPdgCodeToA(kPdgTgtO16)  );
  numu_H1.    UseSplines();
  numu_O16.   UseSplines();
  numubar_H1. UseSplines();
  numubar_O16.UseSplines();
  nue_H1.     UseSplines();
  nue_O16.    UseSplines();
  nuebar_H1.  UseSplines();
  nuebar_O16. UseSplines();

  //
  // Instruct the drivers to sum-up the cross-section splines for all simulated processes
  // (compute total cross-section)
  //

  cout << "Calculating total interaction cross-sections" << endl;

  int    splNumKnots = 300;    // number of knots in total cross section spline
  double splEvMin    =  0.010; // min energy in the spline
  double splEvMax    = 50.000; // max energy in the spline
  bool   splInLogE   = true;   // spline knots distributed logarithmicaly in energy
  numu_H1.    CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);
  numu_O16.   CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);
  numubar_H1. CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);
  numubar_O16.CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);
  nue_H1.     CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);
  nue_O16.    CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);
  nuebar_H1.  CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);
  nuebar_O16. CreateXSecSumSpline (splNumKnots, splEvMin, splEvMax, splInLogE);

  //
  // Write out the cross-section table in the format required by SKDETSIM
  //

  cout << "Writing out the GENIE cross-section table" << endl;

  ofstream outf("./genie_sk_xsec_table.dat",ios::out);

  double w_H1  = 2;  // # of H in H2O
  double w_O16 = 1;  // # of O in H2O

  double Ev = EvMin;
  while(1) {
    double xsec_numu_H1     = numu_H1.     XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);
    double xsec_numu_O16    = numu_O16.    XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);
    double xsec_numubar_H1  = numubar_H1.  XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);
    double xsec_numubar_O16 = numubar_O16. XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);
    double xsec_nue_H1      = nue_H1.      XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);
    double xsec_nue_O16     = nue_O16.     XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);
    double xsec_nuebar_H1   = nuebar_H1.   XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);
    double xsec_nuebar_O16  = nuebar_O16.  XSecSumSpline() -> Evaluate(Ev) / (1E-38 * units::cm2);

    double xsec_numu_H20    = w_H1 * xsec_numu_H1    + w_O16 * xsec_numu_O16     ;
    double xsec_numubar_H20 = w_H1 * xsec_numubar_H1 + w_O16 * xsec_numubar_O16  ;
    double xsec_nue_H20     = w_H1 * xsec_nue_H1     + w_O16 * xsec_nue_O16      ;
    double xsec_nuebar_H20  = w_H1 * xsec_nuebar_H1  + w_O16 * xsec_nuebar_O16   ;

    outf << setfill(' ') << setw(10) << std::fixed << setprecision(5) << Ev
         << setfill(' ') << setw(12) << std::fixed << setprecision(5) << xsec_numu_H20
         << setfill(' ') << setw(12) << std::fixed << setprecision(5) << xsec_numubar_H20
         << setfill(' ') << setw(12) << std::fixed << setprecision(5) << xsec_nue_H20
         << setfill(' ') << setw(12) << std::fixed << setprecision(5) << xsec_nuebar_H20
         << endl;

    Ev += dEv;
    if(Ev > EvMax) break;
  }

  outf.close();

  cout << "Done!" << endl;
}
