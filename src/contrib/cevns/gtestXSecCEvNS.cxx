//____________________________________________________________________________
/*!

\program gtestXSecCEvNS

\brief   test program used for testing the CEvNS cross-section calculation

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created July 15, 2019

\cpright Copyright (c) 2003-2019, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include <TFile.h>
#include <TGraph.h>

#include "Framework/Algorithm/Algorithm.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Utils/RunOpt.h"

using namespace genie;
//__________________________________________________________________________
int main(int argc, char ** argv)
{
  // <temp/>
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);
  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("test", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();
  // </temp>

  LOG("test",pINFO) << "Testing the CEvNS cross-section calculation";

  // Instantiate the coherent elastic cross-section model
  AlgFactory * algf = AlgFactory::Instance();
  AlgId id("genie::PattonCEvNSPXSec","Default");
  const Algorithm * alg = algf->GetAlgorithm(id);
  LOG("test", pINFO) << *alg;
  const XSecAlgorithmI * xsecalg = dynamic_cast<const XSecAlgorithmI*>(alg);

  const int target = 1000180400;
  const int probe  = 12;
  double E         = 0.010;  // GeV
  double Q2        = 0.0001; //GeV^2

  Interaction * interaction = Interaction::CEvNS(target, probe, E);
  interaction->KinePtr()->SetQ2(Q2);

/*
  // test differential cross-section calculation at a single point
  double dxsec = xsecalg->XSec(interaction,kPSQ2fE);

  LOG("test", pNOTICE)
   << "dxsec[CEvNS; target = " << target << "]/dQ2"
   << "(E = " << E << " GeV, Q2 = " << Q2 << " GeV^2) = "
   << dxsec/(units::cm2) << " cm^2/GeV^2";
*/

/*
  // test integrated cross-section calculation at a single point
  double xsec = xsecalg->Integral(interaction);

  LOG("test", pNOTICE)
   << "xsec[CEvNS; target = " << target << "]"
   << "(E = " << E << " GeV) = " << xsec/(units::cm2) << " cm^2";
*/

  // test integrated cross-section calculation at a range of energies
  // and compare with published predictions
  const int n = 100;
  const double Emin = 0.005;
  const double Emax = 0.055;
  const double dE   = (Emax-Emin)/(n-1);

  double e_array[n] = {0};
  double genie_xsec_array[n] = {0};

  for(int i=0; i<n; i++) {
    double E_current = Emin + i*dE;
    interaction->InitStatePtr()->SetProbeE(E_current);
    double xsec_current = xsecalg->Integral(interaction)/(units::cm2);
    e_array[i]          = E_current;
    genie_xsec_array[i] = xsec_current;
    LOG("test", pNOTICE)
     << "xsec[CEvNS; target = " << target << "]"
     << "(E = " << E_current << " GeV) = " << xsec_current << " cm^2";
  }
  TGraph * genie_xsec = new TGraph(n, e_array, genie_xsec_array);
  TGraph * published_xsec = new TGraph("$GENIE/src/contrib/cevns/cevns_arXiv180309183_Ar40_prediction.data");
  published_xsec->SetMarkerStyle(20);
  published_xsec->SetMarkerColor(kRed);
  genie_xsec->SetLineColor(kBlue);
  TFile f("cevns.root","recreate");
  genie_xsec->Write("genie_xsec_nuAr40");
  published_xsec->Write("published_xsec_nuAr40");
  f.Close();
  return 0;
}
//__________________________________________________________________________
