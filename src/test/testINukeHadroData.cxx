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

int main(int /*argc*/, char ** /*argv*/)
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
  inukd -> XSecPipP_Elas () -> SaveAsROOT("hdxs.root", "said_pip_p_elas", true);
  inukd -> XSecPipP_Reac () -> SaveAsROOT("hdxs.root", "said_pip_p_reac", false);
  inukd -> XSecPipD_Abs  () -> SaveAsROOT("hdxs.root", "said_pip_d_abs",  false);
  inukd -> XSecPimP_Elas () -> SaveAsROOT("hdxs.root", "said_pim_p_elas", false);
  inukd -> XSecPimP_Reac () -> SaveAsROOT("hdxs.root", "said_pim_p_reac", false);
  inukd -> XSecPimP_CEx  () -> SaveAsROOT("hdxs.root", "said_pim_p_cex",  false);
  inukd -> XSecPP_Elas   () -> SaveAsROOT("hdxs.root", "said_p_p_elas",   false);
  inukd -> XSecPP_Reac   () -> SaveAsROOT("hdxs.root", "said_p_p_reac",   false);
  inukd -> XSecNP_Elas   () -> SaveAsROOT("hdxs.root", "said_n_p_elas",   false);
  inukd -> XSecNP_Reac   () -> SaveAsROOT("hdxs.root", "said_n_p_reac",   false);

  // Cross sections from Mashnik's calculations (p+Fe)
  //
  inukd -> XSecPFe_Elas  () -> SaveAsROOT("hdxs.root", "mashnik_pfe_elas", false);
  inukd -> XSecPFe_Reac  () -> SaveAsROOT("hdxs.root", "mashnik_pfe_reac", false);
  inukd -> XSecPFe_P     () -> SaveAsROOT("hdxs.root", "mashnik_pfe_p",    false);
  inukd -> XSecPFe_N     () -> SaveAsROOT("hdxs.root", "mashnik_pfe_n",    false);
  inukd -> XSecPFe_PP    () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pp",   false);
  inukd -> XSecPFe_NP    () -> SaveAsROOT("hdxs.root", "mashnik_pfe_np",   false);
  inukd -> XSecPFe_NPP   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_npp",  false);
  inukd -> XSecPFe_NNP   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_nnp",  false);
  inukd -> XSecPFe_NNPP  () -> SaveAsROOT("hdxs.root", "mashnik_pfe_nnpp", false);
  inukd -> XSecPFe_Pim   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pim",  false);
  inukd -> XSecPFe_Pi0   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pi0",  false);
  inukd -> XSecPFe_Pip   () -> SaveAsROOT("hdxs.root", "mashnik_pfe_pip",  false);

  // Cross sections from Mashnik's calculations (pi+Fe)
  //
  inukd -> XSecPiFe_P    () -> SaveAsROOT("hdxs.root", "mashnik_pife_p",    false);
  inukd -> XSecPiFe_PP   () -> SaveAsROOT("hdxs.root", "mashnik_pife_pp",   false);
  inukd -> XSecPiFe_PPP  () -> SaveAsROOT("hdxs.root", "mashnik_pife_ppp",  false);
  inukd -> XSecPiFe_N    () -> SaveAsROOT("hdxs.root", "mashnik_pife_n",    false);
  inukd -> XSecPiFe_NN   () -> SaveAsROOT("hdxs.root", "mashnik_pife_nn",   false);
  inukd -> XSecPiFe_NNN  () -> SaveAsROOT("hdxs.root", "mashnik_pife_nnn",  false);
  inukd -> XSecPiFe_NP   () -> SaveAsROOT("hdxs.root", "mashnik_pife_np",   false);
  inukd -> XSecPiFe_NPP  () -> SaveAsROOT("hdxs.root", "mashnik_pife_npp",  false);
  inukd -> XSecPiFe_NPPP () -> SaveAsROOT("hdxs.root", "mashnik_pife_nppp", false);
  inukd -> XSecPiFe_NNP  () -> SaveAsROOT("hdxs.root", "mashnik_pife_nnp",  false);
  inukd -> XSecPiFe_NNPP () -> SaveAsROOT("hdxs.root", "mashnik_pife_nnpp", false);
  inukd -> XSecPiFe_Pi0  () -> SaveAsROOT("hdxs.root", "mashnik_pife_pi0",  false);

  // Hadronic x-sections from experimental data (Ash & Carrol)
  //
  inukd -> XSecAshPiFe_Abs () -> SaveAsROOT("hdxs.root", "exp_ash_pife_abs",  false);
  inukd -> XSecAshPiFe_Reac() -> SaveAsROOT("hdxs.root", "exp_ash_pife_reac", false);
  inukd -> XSecCarPiFe_Tot () -> SaveAsROOT("hdxs.root", "exp_car_pife_tot",  false);

  // Total x-sections needed from mean free path
  //
  inukd -> XSecPip_Tot() -> SaveAsROOT("hdxs.root", "tot_pip",  false);
  inukd -> XSecPim_Tot() -> SaveAsROOT("hdxs.root", "tot_pim",  false);
  inukd -> XSecPi0_Tot() -> SaveAsROOT("hdxs.root", "tot_pi0",  false);
  inukd -> XSecP_Tot  () -> SaveAsROOT("hdxs.root", "tot_p",    false);
  inukd -> XSecN_Tot  () -> SaveAsROOT("hdxs.root", "tot_n",    false);

  // Fractions of total x-sections needed for particle fates
  //
  inukd -> FracPip_CEx    () -> SaveAsROOT("hdxs.root", "fr_pip_cex",    false);
  inukd -> FracPip_Elas   () -> SaveAsROOT("hdxs.root", "fr_pip_elas",   false);
  inukd -> FracPip_Reac   () -> SaveAsROOT("hdxs.root", "fr_pip_reac",   false);
  inukd -> FracPip_Abs    () -> SaveAsROOT("hdxs.root", "fr_pip_abs",    false);
  inukd -> FracPim_CEx    () -> SaveAsROOT("hdxs.root", "fr_pim_cex",    false);
  inukd -> FracPim_Elas   () -> SaveAsROOT("hdxs.root", "fr_pim_elas",   false);
  inukd -> FracPim_Reac   () -> SaveAsROOT("hdxs.root", "fr_pim_reac",   false);
  inukd -> FracPim_Abs    () -> SaveAsROOT("hdxs.root", "fr_pim_abs",    false);
  inukd -> FracPi0_CEx    () -> SaveAsROOT("hdxs.root", "fr_pi0_cex",    false);
  inukd -> FracPi0_Elas   () -> SaveAsROOT("hdxs.root", "fr_pi0_elas",   false);
  inukd -> FracPi0_Reac   () -> SaveAsROOT("hdxs.root", "fr_pi0_reac",   false);
  inukd -> FracPi0_Abs    () -> SaveAsROOT("hdxs.root", "fr_pi0_abs",    false);
  inukd -> FracP_Reac     () -> SaveAsROOT("hdxs.root", "fr_p_reac",     false);
  inukd -> FracN_Reac     () -> SaveAsROOT("hdxs.root", "fr_n_reac",     false);
  inukd -> FracPiA_Elas   () -> SaveAsROOT("hdxs.root", "fr_pia_elas",   false);
  inukd -> FracPiA_Inel   () -> SaveAsROOT("hdxs.root", "fr_pia_inel",   false);
  inukd -> FracPiA_CEx    () -> SaveAsROOT("hdxs.root", "fr_pia_cex",    false);
  inukd -> FracPiA_Abs    () -> SaveAsROOT("hdxs.root", "fr_pia_abs",    false);
  inukd -> FracPiA_PP     () -> SaveAsROOT("hdxs.root", "fr_pia_pp",     false);
  inukd -> FracPiA_NPP    () -> SaveAsROOT("hdxs.root", "fr_pia_npp",    false);
  inukd -> FracPiA_NNP    () -> SaveAsROOT("hdxs.root", "fr_pia_nnp",    false);
  inukd -> FracPiA_4N4P   () -> SaveAsROOT("hdxs.root", "fr_pia_4n4p",   false);
  inukd -> FracPiA_PiProd () -> SaveAsROOT("hdxs.root", "fr_pia_piprod", false);
  inukd -> FracPA_Elas    () -> SaveAsROOT("hdxs.root", "fr_pa_elas",    false);
  inukd -> FracPA_Inel    () -> SaveAsROOT("hdxs.root", "fr_pa_inel",    false);
  inukd -> FracPA_Abs     () -> SaveAsROOT("hdxs.root", "fr_pa_abs",     false);
  inukd -> FracPA_PP      () -> SaveAsROOT("hdxs.root", "fr_pa_pp",      false);
  inukd -> FracPA_NPP     () -> SaveAsROOT("hdxs.root", "fr_pa_npp",     false);
  inukd -> FracPA_NNP     () -> SaveAsROOT("hdxs.root", "fr_pa_nnp",     false);
  inukd -> FracPA_4N4P    () -> SaveAsROOT("hdxs.root", "fr_pa_4n4p",    false);
  inukd -> FracPA_PiProd  () -> SaveAsROOT("hdxs.root", "fr_pa_piprod",  false); 
}
