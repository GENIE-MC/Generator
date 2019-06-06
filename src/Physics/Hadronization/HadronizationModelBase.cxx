//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TLorentzVector.h>
#include <TH1D.h>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Physics/Hadronization/HadronizationModelBase.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

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
double HadronizationModelBase::Wmin(void) const
{
  return (kNucleonMass+kPionMass);
}
//____________________________________________________________________________
double HadronizationModelBase::MaxMult(const Interaction * interaction) const
{
  double W = interaction->Kine().W();

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

  const InitialState & init_state = interaction->InitState();
  int probe_pdg = init_state.ProbePdg();
  int nuc_pdg   = init_state.Tgt().HitNucPdg();

  const ProcessInfo & proc_info = interaction->ProcInfo();
  bool is_CC = proc_info.IsWeakCC();
  bool is_NC = proc_info.IsWeakNC();
  bool is_EM = proc_info.IsEM();
  // EDIT
  bool is_dm = proc_info.IsDarkMatter();

  //
  // get the R2, R3 factors
  //

  double R2=1., R3=1.;

  // weak CC or NC case
  // EDIT
  if(is_CC || is_NC || is_dm) {
     bool is_nu    = pdg::IsNeutrino     (probe_pdg); 
     bool is_nubar = pdg::IsAntiNeutrino (probe_pdg);
     bool is_p     = pdg::IsProton       (nuc_pdg);
     bool is_n     = pdg::IsNeutron      (nuc_pdg);
     bool is_dmi   = pdg::IsDarkMatter   (probe_pdg);  // EDIT
     if((is_nu && is_p) || (is_dmi && is_p))  {
         R2 = (is_CC) ? fRvpCCm2 : fRvpNCm2;
         R3 = (is_CC) ? fRvpCCm3 : fRvpNCm3;
     } else 
      if((is_nu && is_n) || (is_dmi && is_n)) {
         R2 = (is_CC) ? fRvnCCm2 : fRvnNCm2;
         R3 = (is_CC) ? fRvnCCm3 : fRvnNCm3;
      } else 
      if(is_nubar && is_p)  {
         R2 = (is_CC) ? fRvbpCCm2 :   fRvbpNCm2;
         R3 = (is_CC) ? fRvbpCCm3 :   fRvbpNCm3;
      } else 
      if(is_nubar && is_n) {
         R2 = (is_CC) ? fRvbnCCm2 : fRvbnNCm2;
         R3 = (is_CC) ? fRvbnCCm3 : fRvbnNCm3;
      } else {
         LOG("BaseHad", pERROR) 
            << "Invalid initial state: " << init_state;
     }
  }//cc||nc?

  // EM case (apply the NC tuning factors)

  if(is_EM) {
     bool is_l     = pdg::IsNegChargedLepton (probe_pdg); 
     bool is_lbar  = pdg::IsPosChargedLepton (probe_pdg);
     bool is_p     = pdg::IsProton           (nuc_pdg);
     bool is_n     = pdg::IsNeutron          (nuc_pdg);
     if(is_l && is_p)  {
         R2 = fRvpNCm2;
         R3 = fRvpNCm3;
      } else 
      if(is_l && is_n) {
         R2 = fRvnNCm2;
         R3 = fRvnNCm3;
      } else 
      if(is_lbar && is_p)  {
         R2 = fRvbpNCm2;
         R3 = fRvbpNCm3;
      } else 
      if(is_lbar && is_n) {
         R2 = fRvbnNCm2;
         R3 = fRvbnNCm3;
      } else {
         LOG("BaseHad", pERROR) 
            << "Invalid initial state: " << init_state;
     }
  }//em?

  //
  // Apply to the multiplicity probability distribution
  //

  int nbins = mp->GetNbinsX();
  for(int i = 1; i <= nbins; i++) {
     int n = TMath::Nint( mp->GetBinCenter(i) ); 

     double R=1;
     if      (n==2) R=R2;
     else if (n==3) R=R3;

     if(n==2 || n==3) {
        double P   = mp->GetBinContent(i);
        double Psc = R*P;
        LOG("BaseHad", pDEBUG) 
          << "n=" << n << "/ Scaling factor R = " 
                              << R << "/ P " << P << " --> " << Psc;
        mp->SetBinContent(i, Psc);
     }
     if(n>3) break;
  }

  // renormalize the histogram?
  if(norm) {
     double histo_norm = mp->Integral("width");
     if(histo_norm>0) mp->Scale(1.0/histo_norm);
  }
}
//____________________________________________________________________________
