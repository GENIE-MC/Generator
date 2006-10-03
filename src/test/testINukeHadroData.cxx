//____________________________________________________________________________
/*!

\program testINukeHadroData

\brief   test program used for testing the hadron cross section data used in
         INTRANUKE.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 12, 2004

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include "HadronTransport/INukeHadroData.h"
#include "Numerical/Spline.h"

using namespace genie;

int main(int argc, char ** argv)
{
  //-- load SAID hadron cross section data

  INukeHadroData * inukd = INukeHadroData::Instance();

  //-- save splines to a ROOT file

  // Note: splines are beeing saved to the ROOT files as genie::Spline
  // objects. To read them back, one in a ROOT session load the GENIE
  // libraries first using the $GENIE/src/scripts/loadlibs.C script
  // eg:
  // root [0] .x $GENIE/src/scripts/loadlibs.C
  // root [1] TFile f("./hdxs.root");
  // root [2] mashnik_pife_pi0->Draw();
  // root [3] said_pip_p_elas->Draw("same");
  //
  
  // Cross sections from SAID (Arndt, Workman, Strakovsky) PWA fit.
  //  
  inukd -> XsPipPElas () -> SaveAsROOT("hdxs.root", "said_pip_p_elas", true);
  inukd -> XsPipPReac () -> SaveAsROOT("hdxs.root", "said_pip_p_reac", false);
  inukd -> XsPipDAbs  () -> SaveAsROOT("hdxs.root", "said_pip_d_abs",  false);
  inukd -> XsPimPElas () -> SaveAsROOT("hdxs.root", "said_pim_p_elas", false);
  inukd -> XsPimPReac () -> SaveAsROOT("hdxs.root", "said_pim_p_reac", false);
  inukd -> XsPimPCEx  () -> SaveAsROOT("hdxs.root", "said_pim_p_cex",  false);
  inukd -> XsPPElas   () -> SaveAsROOT("hdxs.root", "said_p_p_elas",   false);
  inukd -> XsPPReac   () -> SaveAsROOT("hdxs.root", "said_p_p_reac",   false);
  inukd -> XsNPElas   () -> SaveAsROOT("hdxs.root", "said_n_p_elas",   false);
  inukd -> XsNPReac   () -> SaveAsROOT("hdxs.root", "said_n_p_reac",   false);

  // Cross sections from Mashnik's calculations (p+Fe)
  //
  inukd -> XsPFeElas  () -> SaveAsROOT("hdxs.root", "mashnik_pfe_elas", false);
  inukd -> XsPFeReac  () -> SaveAsROOT("hdxs.root", "mashnik_pfe_reac", false);
  inukd -> XsPFeP     () -> SaveAsROOT("hdxs.root", "mashnik_pfe_p",    false);
  inukd -> XsPFePP    () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pp",   false);
  inukd -> XsPFeNPP   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_npp",  false);
  inukd -> XsPFeNNP   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_nnp",  false);
  inukd -> XsPFeNNPP  () -> SaveAsROOT("hdxs.root", "mashnik_pfe_nnpp", false);
  inukd -> XsPFePim   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pim",  false);
  inukd -> XsPFePi0   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pi0",  false);
  inukd -> XsPFePip   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pip",  false);

  // Cross sections from Mashnik's calculations (pi+Fe)
  //
  inukd -> XsPiFeP    () -> SaveAsROOT("hdxs.root", "mashnik_pife_p",    false);
  inukd -> XsPiFePP   () -> SaveAsROOT("hdxs.root", "mashnik_pife_pp",   false);
  inukd -> XsPiFePPP  () -> SaveAsROOT("hdxs.root", "mashnik_pife_ppp",  false);
  inukd -> XsPiFeN    () -> SaveAsROOT("hdxs.root", "mashnik_pife_n",    false);
  inukd -> XsPiFeNN   () -> SaveAsROOT("hdxs.root", "mashnik_pife_nn",   false);
  inukd -> XsPiFeNNN  () -> SaveAsROOT("hdxs.root", "mashnik_pife_nnn",  false);
  inukd -> XsPiFeNP   () -> SaveAsROOT("hdxs.root", "mashnik_pife_np",   false);
  inukd -> XsPiFeNPP  () -> SaveAsROOT("hdxs.root", "mashnik_pife_npp",  false);
  inukd -> XsPiFeNPPP () -> SaveAsROOT("hdxs.root", "mashnik_pife_nppp", false);
  inukd -> XsPiFeNNP  () -> SaveAsROOT("hdxs.root", "mashnik_pife_nnp",  false);
  inukd -> XsPiFeNNPP () -> SaveAsROOT("hdxs.root", "mashnik_pife_nnpp", false);
  inukd -> XsPiFePi0  () -> SaveAsROOT("hdxs.root", "mashnik_pife_pi0",  false);
}
