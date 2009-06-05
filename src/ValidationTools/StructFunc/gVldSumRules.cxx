//____________________________________________________________________________
/*!

\program gvld_structfunc_sumrules

\brief   A GENIE utility that generates structure function comparison plots

         Syntax:
           gvld_structfunc_sumrules 
             -r res_model -d dis_model -c charm_dis_model -p lepton -t target

         Options:
		      
\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 30, 2009

\cpright Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
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
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
#include "ValidationTools/StructFunc/ExtractStructFunc.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::vld_structfunc;

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
        vld_structfunc::ExtractStructFunc(xcurr,Q2curr,F1,F2,xF3);
    }//x



  }//Q2

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
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing commad line arguments";

  // get input GENIE cross section file
  try {
    gOptSfFilename_curr = utils::clap::CmdLineArgAsString(argc,argv,'f');
    bool ok = CheckRootFilename(gOptSfFilename_curr.c_str());
    if(!ok) {
      PrintSyntax();
      exit(1);
    }
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      PrintSyntax();
      exit(1);
    }
  }

  // get [reference] input GENIE cross section file
  try {
    gOptSfFilename_ref0 = utils::clap::CmdLineArgAsString(argc,argv,'r');
    bool ok = CheckRootFilename(gOptSfFilename_ref0.c_str());
    if(!ok) {
      PrintSyntax();
      exit(1);
    }
    gOptHaveRef = true;
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gvldtest", pNOTICE) << "No reference cross section file";
      gOptHaveRef = false;
    }
  }

  // check whether to compare with data
  gOptCmpWithData = genie::utils::clap::CmdLineArgAsBool(argc,argv,'d');

  if(gOptCmpWithData) {

   // get DB URL
   try {
     gOptDbURL = utils::clap::CmdLineArgAsString(argc,argv,'h');
   } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
       gOptDbURL = kDefDbURL;
     }
   }

   // get DB username
   try {
     gOptDbUser = utils::clap::CmdLineArgAsString(argc,argv,'u');
   } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
       PrintSyntax();
       exit(1);
     }
   }

   // get DB passwd
   try {
     gOptDbPasswd = utils::clap::CmdLineArgAsString(argc,argv,'p');
   } catch(exceptions::CmdLineArgParserException e) {
     if(!e.ArgumentFound()) {
       PrintSyntax();
       exit(1);
     }
   }
  } // -d enabled?
}
//_________________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gvldtest", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gNuSampleTest  -f sample.root [-n nev] [-r reference_sample.root]\n";
}
//_________________________________________________________________________________
