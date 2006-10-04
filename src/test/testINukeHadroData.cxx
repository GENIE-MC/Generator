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

  // Hadronic x-sections from experimental data (Ash & Carrol)
  //
  inukd -> XsAshPiFeAbs () -> SaveAsROOT("hdxs.root", "exp_ash_pife_abs",  false);
  inukd -> XsAshPiFeReac() -> SaveAsROOT("hdxs.root", "exp_ash_pife_reac", false);
  inukd -> XsCarPiFeTot () -> SaveAsROOT("hdxs.root", "exp_car_pife_tot",  false);

  // Total x-sections needed from mean free path
  //
  inukd -> XsPipTot() -> SaveAsROOT("hdxs.root", "tot_pip",  false);
  inukd -> XsPimTot() -> SaveAsROOT("hdxs.root", "tot_pim",  false);
  inukd -> XsPi0Tot() -> SaveAsROOT("hdxs.root", "tot_pi0",  false);
  inukd -> XsPTot  () -> SaveAsROOT("hdxs.root", "tot_p",    false);
  inukd -> XsNTot  () -> SaveAsROOT("hdxs.root", "tot_n",    false);

/*
  // Fractions of total x-sections needed for particle fates
  //
  inukd -> FrPipCEx    () -> SaveAsROOT("hdxs.root", "fr_pip_cex",    false);
  inukd -> FrPipElas   () -> SaveAsROOT("hdxs.root", "fr_pip_elas",   false);
  inukd -> FrPipReac   () -> SaveAsROOT("hdxs.root", "fr_pip_reac",   false);
  inukd -> FrPipAbs    () -> SaveAsROOT("hdxs.root", "fr_pip_abs",    false);
  inukd -> FrPimCEx    () -> SaveAsROOT("hdxs.root", "fr_pim_cex",    false);
  inukd -> FrPimElas   () -> SaveAsROOT("hdxs.root", "fr_pim_elas",   false);
  inukd -> FrPimReac   () -> SaveAsROOT("hdxs.root", "fr_pim_reac",   false);
  inukd -> FrPimAbs    () -> SaveAsROOT("hdxs.root", "fr_pim_abs",    false);
  inukd -> FrPi0CEx    () -> SaveAsROOT("hdxs.root", "fr_pi0_cex",    false);
  inukd -> FrPi0Elas   () -> SaveAsROOT("hdxs.root", "fr_pi0_elas",   false);
  inukd -> FrPi0Reac   () -> SaveAsROOT("hdxs.root", "fr_pi0_reac",   false);
  inukd -> FrPi0Abs    () -> SaveAsROOT("hdxs.root", "fr_pi0_abs",    false);
  inukd -> FrPReac     () -> SaveAsROOT("hdxs.root", "fr_p_reac",     false);
  inukd -> FrNReac     () -> SaveAsROOT("hdxs.root", "fr_n_reac",     false);
  inukd -> FrPiAElas   () -> SaveAsROOT("hdxs.root", "fr_pia_elas",   false);
  inukd -> FrPiAInel   () -> SaveAsROOT("hdxs.root", "fr_pia_inel",   false);
  inukd -> FrPiACEx    () -> SaveAsROOT("hdxs.root", "fr_pia_cex",    false);
  inukd -> FrPiAAbs    () -> SaveAsROOT("hdxs.root", "fr_pia_abs",    false);
  inukd -> FrPiAPP     () -> SaveAsROOT("hdxs.root", "fr_pia_pp",     false);
  inukd -> FrPiANPP    () -> SaveAsROOT("hdxs.root", "fr_pia_npp",    false);
  inukd -> FrPiANNP    () -> SaveAsROOT("hdxs.root", "fr_pia_nnp",    false);
  inukd -> FrPiA4N4P   () -> SaveAsROOT("hdxs.root", "fr_pia_4n4p",   false);
  inukd -> FrPiAPiProd () -> SaveAsROOT("hdxs.root", "fr_pia_piprod", false);
  inukd -> FrPAElas    () -> SaveAsROOT("hdxs.root", "fr_pa_elas",    false);
  inukd -> FrPAInel    () -> SaveAsROOT("hdxs.root", "fr_pa_inel",    false);
  inukd -> FrPAAbs     () -> SaveAsROOT("hdxs.root", "fr_pa_abs",     false);
  inukd -> FrPAPP      () -> SaveAsROOT("hdxs.root", "fr_pa_pp",      false);
  inukd -> FrPANPP     () -> SaveAsROOT("hdxs.root", "fr_pa_npp",     false);
  inukd -> FrPANNP     () -> SaveAsROOT("hdxs.root", "fr_pa_nnp",     false);
  inukd -> FrPA4N4P    () -> SaveAsROOT("hdxs.root", "fr_pa_4n4p",    false);
  inukd -> FrPAPiProd  () -> SaveAsROOT("hdxs.root", "fr_pa_piprod",  false); */
}
