//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - August 17, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

//#include <TSystem.h>
#include <TLorentzVector.h>
//#include <TClonesArray.h>
//#include <TMCParticle6.h>
#include <TH1D.h>
#include <TMath.h>
//#include <TF1.h>

//#include "Algorithm/AlgConfigPool.h"
//#include "Algorithm/AlgFactory.h"
#include "Conventions/Constants.h"
//#include "Conventions/Controls.h"
//#include "Decay/DecayModelI.h"
#include "Fragmentation/HadronizationModelBase.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
//#include "Numerical/RandomGen.h"
//#include "Numerical/Spline.h"
//#include "PDG/PDGLibrary.h"
//#include "PDG/PDGCodeList.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
//#include "Utils/KineUtils.h"
//#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
//using namespace genie::controls;
//using namespace genie::utils::print;

//____________________________________________________________________________
HadronizationModelBase::HadronizationModelBase(void) :
HadronizationModelI()
{

}
//____________________________________________________________________________
HadronizationModelBase::HadronizationModelBase(string name) :
HadronizationModelI(name)
{

}
//____________________________________________________________________________
HadronizationModelBase::HadronizationModelBase(string name, string config) :
HadronizationModelI(name, config)
{

}
//____________________________________________________________________________
HadronizationModelBase::~HadronizationModelBase()
{

}
//____________________________________________________________________________
bool HadronizationModelBase::AssertWMin(const Interaction * interaction) const
{
  double W = interaction->GetKinematics().W();

  if(W <= kNucleonMass+kPionMass) {
     LOG("KNOHad", pWARN)  << "Low invariant mass, W = " << W << " GeV!!";
     return false;
  }
  return true;
}
//____________________________________________________________________________
double HadronizationModelBase::MaxMult(const Interaction * interaction) const
{
  double W = interaction->GetKinematics().W();

  double maxmult = TMath::Floor(1 + (W-kNeutronMass)/kPionMass);
  return maxmult;
}
//____________________________________________________________________________
TH1D * HadronizationModelBase::CreateMultProbHist(double maxmult) const
{
  double minmult = 2;
  int    nbins   = TMath::Nint(maxmult-minmult+1);

  TH1D * mult_prob = new TH1D("mult_prob", 
      "hadronic multiplicity distribution", nbins, minmult-0.5, maxmult+0.5);
  mult_prob->SetDirectory(0);

  return mult_prob;
}
//____________________________________________________________________________
void HadronizationModelBase::ApplyRijk(
		  const Interaction * interaction, bool norm, TH1D * mp) const
{
// Apply the NEUGEN multiplicity probability scaling factors
//
  if(!mp) return;

  const InitialState & init_state = interaction->GetInitialState();
  int nu_pdg  = init_state.GetProbePDGCode();
  int nuc_pdg = init_state.GetTarget().StruckNucleonPDGCode();

  double R2=1., R3=1.;

  const ProcessInfo & proc_info = interaction->GetProcessInfo();
  bool isCC = proc_info.IsWeakCC();

  if(pdg::IsNeutrino(nu_pdg) && pdg::IsProton(nuc_pdg))  {
    R2 = (isCC) ? fRvpCCm2 : fRvpNCm2;
    R3 = (isCC) ? fRvpCCm3 : fRvpNCm3;
  } else 
  if(pdg::IsNeutrino(nu_pdg) && pdg::IsNeutron(nuc_pdg)) {
    R2 = (isCC) ? fRvnCCm2 : fRvnNCm2;
    R3 = (isCC) ? fRvnCCm3 : fRvnNCm3;
  } else 
  if(pdg::IsAntiNeutrino(nu_pdg) && pdg::IsProton(nuc_pdg))  {
    R2 = (isCC) ? fRvbpCCm2 :   fRvbpNCm2;
    R3 = (isCC) ? fRvbpCCm3 :   fRvbpNCm3;
  } else 
  if(pdg::IsAntiNeutrino(nu_pdg) && pdg::IsNeutron(nuc_pdg)) {
    R2 = (isCC) ? fRvbnCCm2 : fRvbnNCm2;
    R3 = (isCC) ? fRvbnCCm3 : fRvbnNCm3;
  } else {
    LOG("KNOHad", pERROR) << "Invalid initial state: " << init_state;
  }

  int nbins = mp->GetNbinsX();
  for(int i = 1; i <= nbins; i++) {
     int n = TMath::Nint( mp->GetBinCenter(i) ); 

     double R=1;
     if      (n==2) R=R2;
     else if (n==3) R=R3;

     if(n==2 || n==3) {
        double P   = mp->GetBinContent(i);
        double Psc = R*P;
        LOG("KNOHad", pDEBUG) 
          << "n=" << n << "/ Scaling factor R = " 
                              << R << "/ P " << P << " --> " << Psc;
        mp->SetBinContent(i, Psc);
     }
     if(n>3) break;
  }

  // renormalize the histogram?
  if(norm) {
     mp->Scale( 1.0 / mp->Integral("width") );
  }
}
//____________________________________________________________________________
