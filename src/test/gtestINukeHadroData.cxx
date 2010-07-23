//____________________________________________________________________________
/*!

\program testINukeHadroData

\brief   test program used for testing the hadron cross section data used in
         INTRANUKE. 
         The program will save all the hadron cross section splines at an output 
         ROOT file (hadxs.root).
         If a hadron kinetic energy is passed as command line argument then,
         rather than saving the cross section splines at an output file, the
         program will print out the cross section values at the specified
         kinetic energy.

\syntax  gtestINukeHadroData [-k energy(MeV)]

         []  denotes an optionan argument
         -k  can be used to pass a kinetic energy for which all the hadron
             cross sections will be printed out

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created May 12, 2004

\cpright Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <math.h>

#include "HadronTransport/INukeHadroData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "Utils/CmdLnArgParser.h"

using namespace genie;

void  SaveSplinesToRootFile (void);
void  PrintOutForInputKE    (double ke);

int main(int argc, char ** argv)
{
  double ke = -1; // input hadron kinetic energy

  CmdLnArgParser parser(argc,argv);
  if( parser.OptionExists('k') ) {
    ke = parser.ArgAsDouble('k');
  }

  if(ke < 0) {
    SaveSplinesToRootFile();
  }
  else {     
    PrintOutForInputKE(ke);
  }
  
  return 0;
}
//____________________________________________________________________________
void SaveSplinesToRootFile(void)
{
  string filename = "hadxs.root";

  LOG("testINukeHadroData", pNOTICE) 
        << "Saving INTRANUKE hadron x-section splines to file: " << filename;

  //-- load SAID hadron cross section data
  INukeHadroData * inukd = INukeHadroData::Instance();

  //-- save splines to a ROOT file
  //
  // Note: splines are beeing saved to the ROOT files as genie::Spline
  // objects. To read them back, one in a ROOT session load the GENIE
  // libraries first using the $GENIE/src/scripts/loadlibs.C script
  // eg:
  // root [0] .x $GENIE/src/scripts/loadlibs.C
  // root [1] TFile f("./hdxs.root");
  // root [2] pN_tot->Draw();
  // root [3] pi0N_tot->Draw("same");
  //
  
  // Write out hN mode hadron x-section splines
  inukd -> XSecPN_Tot       () -> SaveAsROOT(filename, "pN_tot",    true);
  inukd -> XSecPN_Elas      () -> SaveAsROOT(filename, "pN_elas",   false);
  inukd -> XSecPN_Reac      () -> SaveAsROOT(filename, "pN_reac",   false);
  inukd -> XSecNN_Tot       () -> SaveAsROOT(filename, "nN_tot",    false);
  inukd -> XSecNN_Elas      () -> SaveAsROOT(filename, "nN_elas",   false);
  inukd -> XSecNN_Reac      () -> SaveAsROOT(filename, "nN_reac",   false);
  inukd -> XSecPipN_Tot     () -> SaveAsROOT(filename, "pipN_tot",  false);
  inukd -> XSecPipN_CEx     () -> SaveAsROOT(filename, "pipN_cex",  false);
  inukd -> XSecPipN_Elas    () -> SaveAsROOT(filename, "pipN_elas", false);
  inukd -> XSecPipN_Reac    () -> SaveAsROOT(filename, "pipN_reac", false);
  inukd -> XSecPipN_Abs     () -> SaveAsROOT(filename, "pipN_abs",  false);
  inukd -> XSecPimN_Tot     () -> SaveAsROOT(filename, "pimN_tot",  false);
  inukd -> XSecPimN_CEx     () -> SaveAsROOT(filename, "pimN_cex",  false);
  inukd -> XSecPimN_Elas    () -> SaveAsROOT(filename, "pimN_elas", false);
  inukd -> XSecPimN_Reac    () -> SaveAsROOT(filename, "pimN_reac", false);
  inukd -> XSecPimN_Abs     () -> SaveAsROOT(filename, "pimN_abs",  false);
  inukd -> XSecPi0N_Tot     () -> SaveAsROOT(filename, "pi0N_tot",  false);
  inukd -> XSecPi0N_CEx     () -> SaveAsROOT(filename, "pi0N_cex",  false);
  inukd -> XSecPi0N_Elas    () -> SaveAsROOT(filename, "pi0N_elas", false);
  inukd -> XSecPi0N_Reac    () -> SaveAsROOT(filename, "pi0N_reac", false);
  inukd -> XSecPi0N_Abs     () -> SaveAsROOT(filename, "pi0N_abs",  false);
  
  // Write out hA mode hadron x-section splines
  inukd -> XSecPA_Tot       () -> SaveAsROOT(filename, "pA_tot",       false);
  inukd -> XSecPA_Elas      () -> SaveAsROOT(filename, "pA_elas",      false);
  inukd -> XSecPA_Inel      () -> SaveAsROOT(filename, "pA_inel",      false);
  inukd -> XSecPA_CEx       () -> SaveAsROOT(filename, "pA_cex",       false);
  inukd -> XSecPA_NP        () -> SaveAsROOT(filename, "pA_np",        false);
  inukd -> XSecPA_PP        () -> SaveAsROOT(filename, "pA_pp",        false);
  inukd -> XSecPA_NPP       () -> SaveAsROOT(filename, "pA_npp",       false);
  inukd -> XSecPA_NNP       () -> SaveAsROOT(filename, "pA_nnp",       false);
  inukd -> XSecPA_NNPPP     () -> SaveAsROOT(filename, "pA_nnppp",     false);
  inukd -> XSecPA_NPip      () -> SaveAsROOT(filename, "pA_pip",       false);
  inukd -> XSecPA_NPipPi0   () -> SaveAsROOT(filename, "pA_npippi0",   false);
  inukd -> XSecNA_Tot       () -> SaveAsROOT(filename, "nA_tot",       false);
  inukd -> XSecNA_Elas      () -> SaveAsROOT(filename, "nA_elas",      false);
  inukd -> XSecNA_Inel      () -> SaveAsROOT(filename, "nA_inel",      false);
  inukd -> XSecNA_CEx       () -> SaveAsROOT(filename, "nA_cex",       false);
  inukd -> XSecNA_NP        () -> SaveAsROOT(filename, "nA_np",        false);
  inukd -> XSecNA_PP        () -> SaveAsROOT(filename, "nA_pp",        false);
  inukd -> XSecNA_NPP       () -> SaveAsROOT(filename, "nA_npp",       false); 
  inukd -> XSecNA_NNP       () -> SaveAsROOT(filename, "nA_nnp",       false);
  inukd -> XSecNA_NNPPP     () -> SaveAsROOT(filename, "nA_nnppp",     false);
  inukd -> XSecNA_NPip      () -> SaveAsROOT(filename, "nA_npip",      false);
  inukd -> XSecNA_NPipPi0   () -> SaveAsROOT(filename, "nA_npippi0",   false);
  inukd -> XSecPipA_Tot     () -> SaveAsROOT(filename, "pipA_tot",     false);
  inukd -> XSecPipA_Elas    () -> SaveAsROOT(filename, "pipA_elas",    false);
  inukd -> XSecPipA_Inel    () -> SaveAsROOT(filename, "pipA_inel",    false);
  inukd -> XSecPipA_CEx     () -> SaveAsROOT(filename, "pipA_cex",     false);
  inukd -> XSecPipA_NP      () -> SaveAsROOT(filename, "pipA_np",      false);  
  inukd -> XSecPipA_PP      () -> SaveAsROOT(filename, "pipA_pp",      false);
  inukd -> XSecPipA_NPP     () -> SaveAsROOT(filename, "pipA_npp",     false);
  inukd -> XSecPipA_NNP     () -> SaveAsROOT(filename, "pipA_nnp",     false);
  inukd -> XSecPipA_NNPP    () -> SaveAsROOT(filename, "pipA_nnpp",    false);
  inukd -> XSecPipA_NPipPi0 () -> SaveAsROOT(filename, "pipA_npippi0", false);
  inukd -> XSecPimA_Tot     () -> SaveAsROOT(filename, "pimA_tot",     false);
  inukd -> XSecPimA_Elas    () -> SaveAsROOT(filename, "pimA_elas",    false);
  inukd -> XSecPimA_Inel    () -> SaveAsROOT(filename, "pimA_inel",    false);
  inukd -> XSecPimA_CEx     () -> SaveAsROOT(filename, "pimA_cex",     false);
  inukd -> XSecPimA_NP      () -> SaveAsROOT(filename, "pimA_np",      false);
  inukd -> XSecPimA_PP      () -> SaveAsROOT(filename, "pimA_pp",      false);
  inukd -> XSecPimA_NPP     () -> SaveAsROOT(filename, "pimA_npp",     false);
  inukd -> XSecPimA_NNP     () -> SaveAsROOT(filename, "pimA_nnp",     false);
  inukd -> XSecPimA_NNPP    () -> SaveAsROOT(filename, "pimA_nnpp",    false);
  inukd -> XSecPimA_NPipPi0 () -> SaveAsROOT(filename, "pimA_npippi0", false);
  inukd -> XSecPi0A_Tot     () -> SaveAsROOT(filename, "pi0A_tot",     false);
  inukd -> XSecPi0A_Elas    () -> SaveAsROOT(filename, "pi0A_elas",    false);
  inukd -> XSecPi0A_Inel    () -> SaveAsROOT(filename, "pi0A_inel",    false);
  inukd -> XSecPi0A_CEx     () -> SaveAsROOT(filename, "pi0A_cex",     false);
  inukd -> XSecPi0A_NP      () -> SaveAsROOT(filename, "pi0A_np",      false);
  inukd -> XSecPi0A_PP      () -> SaveAsROOT(filename, "pi0A_pp",      false);
  inukd -> XSecPi0A_NPP     () -> SaveAsROOT(filename, "pi0A_npp",     false);
  inukd -> XSecPi0A_NNP     () -> SaveAsROOT(filename, "pi0A_nnp",     false);
  inukd -> XSecPi0A_NNPP    () -> SaveAsROOT(filename, "pi0A_nnpp",    false);
  inukd -> XSecPi0A_NPipPi0 () -> SaveAsROOT(filename, "pi0A_npippi0", false);
}
//____________________________________________________________________________
void PrintOutForInputKE(double ke)
{
  LOG("testINukeHadroData", pNOTICE) 
          << "Printing out INTRANUKE hadron x-sections";
 
  //-- load SAID hadron cross section data
  INukeHadroData * inukd = INukeHadroData::Instance();

  //-- Print out the hN mode hadron x-section 
  LOG("testINukeHadroData", pNOTICE)  
     << "\n hN mode x-sections:"
     << "\n XSec[pN/tot]      (K=" << ke << " MeV) = " << inukd -> XSecPN_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pN/elas]     (K=" << ke << " MeV) = " << inukd -> XSecPN_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pN/reac]     (K=" << ke << " MeV) = " << inukd -> XSecPN_Reac      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nN/tot]      (K=" << ke << " MeV) = " << inukd -> XSecNN_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nN/elas]     (K=" << ke << " MeV) = " << inukd -> XSecNN_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nN/reac]     (K=" << ke << " MeV) = " << inukd -> XSecNN_Reac      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipN/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPipN_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipN/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPipN_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipN/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPipN_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipN/reac]   (K=" << ke << " MeV) = " << inukd -> XSecPipN_Reac    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipN/abs]    (K=" << ke << " MeV) = " << inukd -> XSecPipN_Abs     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimN/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPimN_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimN/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPimN_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimN/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPimN_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimN/reac]   (K=" << ke << " MeV) = " << inukd -> XSecPimN_Reac    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimN/abs]    (K=" << ke << " MeV) = " << inukd -> XSecPimN_Abs     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0N/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPi0N_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0N/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPi0N_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0N/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPi0N_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0N/reac]   (K=" << ke << " MeV) = " << inukd -> XSecPi0N_Reac    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0N/abs]    (K=" << ke << " MeV) = " << inukd -> XSecPi0N_Abs     () -> Evaluate(ke) << " mbarn"
     << "\n"
     << "\n hA mode x-sections:"
     << "\n XSec[pA/tot]      (K=" << ke << " MeV) = " << inukd -> XSecPA_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/elas]     (K=" << ke << " MeV) = " << inukd -> XSecPA_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/inel]     (K=" << ke << " MeV) = " << inukd -> XSecPA_Inel      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/cex]      (K=" << ke << " MeV) = " << inukd -> XSecPA_CEx       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/np]       (K=" << ke << " MeV) = " << inukd -> XSecPA_NP        () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/pp]       (K=" << ke << " MeV) = " << inukd -> XSecPA_PP        () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/npp]      (K=" << ke << " MeV) = " << inukd -> XSecPA_NPP       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/nnp]      (K=" << ke << " MeV) = " << inukd -> XSecPA_NNP       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/nnppp]    (K=" << ke << " MeV) = " << inukd -> XSecPA_NNPPP     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/npip]     (K=" << ke << " MeV) = " << inukd -> XSecPA_NPip      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pA/npippi0]  (K=" << ke << " MeV) = " << inukd -> XSecPA_NPipPi0   () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/tot]      (K=" << ke << " MeV) = " << inukd -> XSecNA_Tot       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/elas]     (K=" << ke << " MeV) = " << inukd -> XSecNA_Elas      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/inel]     (K=" << ke << " MeV) = " << inukd -> XSecNA_Inel      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/reac]     (K=" << ke << " MeV) = " << inukd -> XSecNA_CEx       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/np]       (K=" << ke << " MeV) = " << inukd -> XSecNA_NP        () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/pp]       (K=" << ke << " MeV) = " << inukd -> XSecNA_PP        () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/npp]      (K=" << ke << " MeV) = " << inukd -> XSecNA_NPP       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/nnp]      (K=" << ke << " MeV) = " << inukd -> XSecNA_NNP       () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/nnppp]    (K=" << ke << " MeV) = " << inukd -> XSecNA_NNPPP     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/npip]     (K=" << ke << " MeV) = " << inukd -> XSecNA_NPip      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[nA/npippi0]  (K=" << ke << " MeV) = " << inukd -> XSecNA_NPipPi0   () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPipA_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPipA_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/inel]   (K=" << ke << " MeV) = " << inukd -> XSecPipA_Inel    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPipA_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/np]     (K=" << ke << " MeV) = " << inukd -> XSecPipA_NP      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/pp]     (K=" << ke << " MeV) = " << inukd -> XSecPipA_PP      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/npp]    (K=" << ke << " MeV) = " << inukd -> XSecPipA_NPP     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/nnp]    (K=" << ke << " MeV) = " << inukd -> XSecPipA_NNP     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pipA/nnpp]   (K=" << ke << " MeV) = " << inukd -> XSecPipA_NNPP    () -> Evaluate(ke) << " mbarn" 
     << "\n XSec[pipA/npippi0](K=" << ke << " MeV) = " << inukd -> XSecPipA_NPipPi0 () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPimA_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPimA_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/inel]   (K=" << ke << " MeV) = " << inukd -> XSecPimA_Inel    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPimA_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/np]     (K=" << ke << " MeV) = " << inukd -> XSecPimA_NP      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/pp]     (K=" << ke << " MeV) = " << inukd -> XSecPimA_PP      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/npp]    (K=" << ke << " MeV) = " << inukd -> XSecPimA_NPP     () -> Evaluate(ke) << " mbarn" 
     << "\n XSec[pimA/nnp]    (K=" << ke << " MeV) = " << inukd -> XSecPimA_NNP     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/nnpp]   (K=" << ke << " MeV) = " << inukd -> XSecPimA_NNPP    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pimA/npippi0](K=" << ke << " MeV) = " << inukd -> XSecPimA_NPipPi0 () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/tot]    (K=" << ke << " MeV) = " << inukd -> XSecPi0A_Tot     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/elas]   (K=" << ke << " MeV) = " << inukd -> XSecPi0A_Elas    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/inel]   (K=" << ke << " MeV) = " << inukd -> XSecPi0A_Inel    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/cex]    (K=" << ke << " MeV) = " << inukd -> XSecPi0A_CEx     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/np]     (K=" << ke << " MeV) = " << inukd -> XSecPi0A_NP      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/pp]     (K=" << ke << " MeV) = " << inukd -> XSecPi0A_PP      () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/npp]    (K=" << ke << " MeV) = " << inukd -> XSecPi0A_NPP     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/nnp]    (K=" << ke << " MeV) = " << inukd -> XSecPi0A_NNP     () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/nnpp]   (K=" << ke << " MeV) = " << inukd -> XSecPi0A_NNPP    () -> Evaluate(ke) << " mbarn"
     << "\n XSec[pi0A/npippi0](K=" << ke << " MeV) = " << inukd -> XSecPi0A_NPipPi0 () -> Evaluate(ke) << " mbarn";
}
//____________________________________________________________________________

