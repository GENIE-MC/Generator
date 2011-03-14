//____________________________________________________________________________
/*!

\program gvld_sf_sum_rule_test

\brief   Check sum rules

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 30, 2009

\cpright Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TGraph.h>
#include <TPostScript.h>
#include <TH1D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TStyle.h>
#include <TLegend.h>

#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"

using std::ostringstream;
using std::string;

using namespace genie;

// function prototypes

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

void TestAdlerSumRule               (void);
void TestGrossLlewellynSmithSumRule (void);
void TestGottfriedSumRule           (void);
void TestBjorkenSumRule             (void);

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);

  TestAdlerSumRule               ();
  TestGrossLlewellynSmithSumRule ();
  TestGottfriedSumRule           ();
  TestBjorkenSumRule             ();

  LOG("gvldtest", pINFO)  << "Done!";

  return 0;
}
//_________________________________________________________________________________
void TestAdlerSumRule(void)
{
//
// ***********************************************
// * Adler sum rule:                             *
// *                                             *
// *          /1                                 *
// * S_{A} =  |dx  (F2{vn} - F2{vp}) / 2x  = 1.  *
// *          /0                                 *
// ***********************************************
//

/*
  const int    kNQ2    = 10;
  const int    kNlogx  = 100;

  const double xmin    = 1E-6;
  const double xmax    = 1;
  const double logxmin = TMath::Log10(xmin);
  const double logxmax = TMath::Log10(xmax);
  const double dlogx   = (logxmax-logxmin)/(kNlogx-1);

  double I  [kNQ2];
  double Q2 [kNQ2] = { 0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0 };

  double x  [kNlogx];
  for(int ix=0; ix<kNlogx; ix++) {
    x[ix] = TMath::Power(10., logxmin + ix * dlogx);
  }

  double F1, F2, xF3;

  for(int iq2=0; iq2<kNQ2; iq2++) {
    double Q2curr = Q2[iq2];

    for(int ix=0; ix<kNlogx; ix++) {
	double xcurr = x[ix];

        //
        // ...
        //
    }//x
  }//Q2
*/

}
//_________________________________________________________________________________
void TestGrossLlewellynSmithSumRule(void)
{
//
// *******************************************
// * Gross - Llewellyn Smith sum rule:       *
// *                                         *
// *            /1                           *
// * S_{GLS} =  | dx xF3{vN} / 2x = 3.       * 
// *           /0                            *
// *******************************************
//


}
//_________________________________________________________________________________
void TestGottfriedSumRule(void)
{
//
// ************************************************
// * Gottfried sum rule:                          *
// *                                              *
// *         /1                                   *
// * S_{G} = | dx (F2{mup} - F2{mun}) / x = 1./3. * 
// *        /0                                    *
// ************************************************
//


}
//_________________________________________________________________________________
void TestBjorkenSumRule(void)
{

}
//_________________________________________________________________________________
void GetCommandLineArgs(int /*argc*/, char ** /*argv*/)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

}
//_________________________________________________________________________________
void PrintSyntax(void)
{

}
//_________________________________________________________________________________
