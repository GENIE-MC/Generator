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
  
  // Write out hN mode hadron x-section splines
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  inukd -> XSecPN_Tot       () -> SaveAsROOT("hdxs.root", "pN_tot",    true);
  inukd -> XSecPN_Elas      () -> SaveAsROOT("hdxs.root", "pN_elas",   false);
  inukd -> XSecPN_Reac      () -> SaveAsROOT("hdxs.root", "pN_reac",   false);
  inukd -> XSecNN_Tot       () -> SaveAsROOT("hdxs.root", "nN_tot",    false);
  inukd -> XSecNN_Elas      () -> SaveAsROOT("hdxs.root", "nN_elas",   false);
  inukd -> XSecNN_Reac      () -> SaveAsROOT("hdxs.root", "nN_reac",   false);
  inukd -> XSecPipN_Tot     () -> SaveAsROOT("hdxs.root", "pipN_tot",  false);
  inukd -> XSecPipN_CEx     () -> SaveAsROOT("hdxs.root", "pipN_cex",  false);
  inukd -> XSecPipN_Elas    () -> SaveAsROOT("hdxs.root", "pipN_elas", false);
  inukd -> XSecPipN_Reac    () -> SaveAsROOT("hdxs.root", "pipN_reac", false);
  inukd -> XSecPipN_Abs     () -> SaveAsROOT("hdxs.root", "pipN_abs",  false);
  inukd -> XSecPimN_Tot     () -> SaveAsROOT("hdxs.root", "pimN_tot",  false);
  inukd -> XSecPimN_CEx     () -> SaveAsROOT("hdxs.root", "pimN_cex",  false);
  inukd -> XSecPimN_Elas    () -> SaveAsROOT("hdxs.root", "pimN_elas", false);
  inukd -> XSecPimN_Reac    () -> SaveAsROOT("hdxs.root", "pimN_reac", false);
  inukd -> XSecPimN_Abs     () -> SaveAsROOT("hdxs.root", "pimN_abs",  false);
  inukd -> XSecPi0N_Tot     () -> SaveAsROOT("hdxs.root", "pi0N_tot",  false);
  inukd -> XSecPi0N_CEx     () -> SaveAsROOT("hdxs.root", "pi0N_cex",  false);
  inukd -> XSecPi0N_Elas    () -> SaveAsROOT("hdxs.root", "pi0N_elas", false);
  inukd -> XSecPi0N_Reac    () -> SaveAsROOT("hdxs.root", "pi0N_reac", false);
  inukd -> XSecPi0N_Abs     () -> SaveAsROOT("hdxs.root", "pi0N_abs",  false);
  
  // Write out hA mode hadron x-section splines
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  inukd -> XSecPA_Tot       () -> SaveAsROOT("hdxs.root", "pA_tot",       false);
  inukd -> XSecPA_Elas      () -> SaveAsROOT("hdxs.root", "pA_elas",      false);
  inukd -> XSecPA_Inel      () -> SaveAsROOT("hdxs.root", "pA_inel",      false);
  inukd -> XSecPA_CEx       () -> SaveAsROOT("hdxs.root", "pA_cex",       false);
  inukd -> XSecPA_NP        () -> SaveAsROOT("hdxs.root", "pA_np",        false);
  inukd -> XSecPA_PP        () -> SaveAsROOT("hdxs.root", "pA_pp",        false);
  inukd -> XSecPA_NPP       () -> SaveAsROOT("hdxs.root", "pA_npp",       false);
  inukd -> XSecPA_NNP       () -> SaveAsROOT("hdxs.root", "pA_nnp",       false);
  inukd -> XSecPA_NNPPP     () -> SaveAsROOT("hdxs.root", "pA_nnppp",     false);
  inukd -> XSecPA_NPip      () -> SaveAsROOT("hdxs.root", "pA_pip",       false);
  inukd -> XSecPA_NPipPi0   () -> SaveAsROOT("hdxs.root", "pA_npippi0",   false);
  inukd -> XSecNA_Tot       () -> SaveAsROOT("hdxs.root", "nA_tot",       false);
  inukd -> XSecNA_Elas      () -> SaveAsROOT("hdxs.root", "nA_elas",      false);
  inukd -> XSecNA_Inel      () -> SaveAsROOT("hdxs.root", "nA_inel",      false);
  inukd -> XSecNA_CEx       () -> SaveAsROOT("hdxs.root", "nA_cex",       false);
  inukd -> XSecNA_NP        () -> SaveAsROOT("hdxs.root", "nA_np",        false);
  inukd -> XSecNA_PP        () -> SaveAsROOT("hdxs.root", "nA_pp",        false);
  inukd -> XSecNA_NPP       () -> SaveAsROOT("hdxs.root", "nA_npp",       false); 
  inukd -> XSecNA_NNP       () -> SaveAsROOT("hdxs.root", "nA_nnp",       false);
  inukd -> XSecNA_NNPPP     () -> SaveAsROOT("hdxs.root", "nA_nnppp",     false);
  inukd -> XSecNA_NPip      () -> SaveAsROOT("hdxs.root", "nA_npip",      false);
  inukd -> XSecNA_NPipPi0   () -> SaveAsROOT("hdxs.root", "nA_npippi0",   false);
  inukd -> XSecPipA_Tot     () -> SaveAsROOT("hdxs.root", "pipA_tot",     false);
  inukd -> XSecPipA_Elas    () -> SaveAsROOT("hdxs.root", "pipA_elas",    false);
  inukd -> XSecPipA_Inel    () -> SaveAsROOT("hdxs.root", "pipA_inel",    false);
  inukd -> XSecPipA_CEx     () -> SaveAsROOT("hdxs.root", "pipA_cex",     false);
  inukd -> XSecPipA_NP      () -> SaveAsROOT("hdxs.root", "pipA_np",      false);  
  inukd -> XSecPipA_PP      () -> SaveAsROOT("hdxs.root", "pipA_pp",      false);
  inukd -> XSecPipA_NPP     () -> SaveAsROOT("hdxs.root", "pipA_npp",     false);
  inukd -> XSecPipA_NNP     () -> SaveAsROOT("hdxs.root", "pipA_nnp",     false);
  inukd -> XSecPipA_NNPP    () -> SaveAsROOT("hdxs.root", "pipA_nnpp",    false);
  inukd -> XSecPipA_NPipPi0 () -> SaveAsROOT("hdxs.root", "pipA_npippi0", false);
  inukd -> XSecPimA_Tot     () -> SaveAsROOT("hdxs.root", "pimA_tot",     false);
  inukd -> XSecPimA_Elas    () -> SaveAsROOT("hdxs.root", "pimA_elas",    false);
  inukd -> XSecPimA_Inel    () -> SaveAsROOT("hdxs.root", "pimA_inel",    false);
  inukd -> XSecPimA_CEx     () -> SaveAsROOT("hdxs.root", "pimA_cex",     false);
  inukd -> XSecPimA_NP      () -> SaveAsROOT("hdxs.root", "pimA_np",      false);
  inukd -> XSecPimA_PP      () -> SaveAsROOT("hdxs.root", "pimA_pp",      false);
  inukd -> XSecPimA_NPP     () -> SaveAsROOT("hdxs.root", "pimA_npp",     false);
  inukd -> XSecPimA_NNP     () -> SaveAsROOT("hdxs.root", "pimA_nnp",     false);
  inukd -> XSecPimA_NNPP    () -> SaveAsROOT("hdxs.root", "pimA_nnpp",    false);
  inukd -> XSecPimA_NPipPi0 () -> SaveAsROOT("hdxs.root", "pimA_npippi0", false);
  inukd -> XSecPi0A_Tot     () -> SaveAsROOT("hdxs.root", "pi0A_tot",     false);
  inukd -> XSecPi0A_Elas    () -> SaveAsROOT("hdxs.root", "pi0A_elas",    false);
  inukd -> XSecPi0A_Inel    () -> SaveAsROOT("hdxs.root", "pi0A_inel",    false);
  inukd -> XSecPi0A_CEx     () -> SaveAsROOT("hdxs.root", "pi0A_cex",     false);
  inukd -> XSecPi0A_NP      () -> SaveAsROOT("hdxs.root", "pi0A_np",      false);
  inukd -> XSecPi0A_PP      () -> SaveAsROOT("hdxs.root", "pi0A_pp",      false);
  inukd -> XSecPi0A_NPP     () -> SaveAsROOT("hdxs.root", "pi0A_npp",     false);
  inukd -> XSecPi0A_NNP     () -> SaveAsROOT("hdxs.root", "pi0A_nnp",     false);
  inukd -> XSecPi0A_NNPP    () -> SaveAsROOT("hdxs.root", "pi0A_nnpp",    false);
  inukd -> XSecPi0A_NPipPi0 () -> SaveAsROOT("hdxs.root", "pi0A_npippi0", false);
}
