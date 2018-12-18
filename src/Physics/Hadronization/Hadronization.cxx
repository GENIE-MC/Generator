//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab
         Shivesh Mandalia <s.p.mandalia \at qmul.ac.uk>
         Queen Mary University of London
 */
//____________________________________________________________________________


#include <algorithm>

#include <TH1D.h>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Hadronization/Hadronization.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
Hadronization::Hadronization() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
Hadronization::Hadronization(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
Hadronization::Hadronization(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
Hadronization::~Hadronization()
{

}
//___________________________________________________________________________
double Hadronization::Wmin(void) const
{
  return (kNucleonMass+kPionMass);
}
//____________________________________________________________________________
double Hadronization::MaxMult(const Interaction * interaction) const
{
  double W = interaction->Kine().W();

  double maxmult = TMath::Floor(1 + (W-kNeutronMass)/kPionMass);
  return maxmult;
}
//____________________________________________________________________________
TH1D * Hadronization::CreateMultProbHist(double maxmult) const
{
  double minmult = 2;
  int    nbins   = TMath::Nint(maxmult-minmult+1);

  TH1D * mult_prob = new TH1D("mult_prob",
      "hadronic multiplicity distribution", nbins, minmult-0.5, maxmult+0.5);
  mult_prob->SetDirectory(0);

  return mult_prob;
}
//____________________________________________________________________________
void Hadronization::ApplyRijk(
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
         LOG("Hadronization", pERROR)
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
         LOG("Hadronization", pERROR)
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
        LOG("Hadronization", pDEBUG)
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
bool Hadronization::AssertValidity(const Interaction * interaction) const
{
  // check that there is no charm production
  // (GENIE uses a special model for these cases)
  if(interaction->ExclTag().IsCharmEvent()) {
     LOG("Hadronization", pWARN) << "Can't hadronize charm events";
     return false;
  }
  // check the available mass
  double W = utils::kinematics::W(interaction);
  if(W < this->Wmin()) {
     LOG("Hadronization", pWARN) << "Low invariant mass, W = "
         << W << " GeV!!";
     return false;
  }

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  if( ! target.HitQrkIsSet() ) {
     LOG("Hadronization", pWARN) << "Hit quark was not set!";
     return false;
  }

  int  probe       = init_state.ProbePdg();
  int  hit_nucleon = target.HitNucPdg();
  int  hit_quark   = target.HitQrkPdg();
//bool from_sea    = target.HitSeaQrk();

  // check hit-nucleon assignment, input neutrino & weak current
  bool isp  = pdg::IsProton           (hit_nucleon);
  bool isn  = pdg::IsNeutron          (hit_nucleon);
  bool isv  = pdg::IsNeutrino         (probe);
  bool isvb = pdg::IsAntiNeutrino     (probe);
  bool isdm = pdg::IsDarkMatter         (probe);
  bool isl  = pdg::IsNegChargedLepton (probe);
  bool islb = pdg::IsPosChargedLepton (probe);
  bool iscc = proc_info.IsWeakCC      ();
  bool isnc = proc_info.IsWeakNC      ();
  bool isdmi = proc_info.IsDarkMatter  ();
  bool isem = proc_info.IsEM          ();
  if( !(iscc||isnc||isem||isdmi) ) {
    LOG("Hadronization", pWARN)
       << "Can only handle electro-weak interactions";
    return false;
  }
  if( !(isp||isn) || !(isv||isvb||isl||islb||isdm) ) {
    LOG("Hadronization", pWARN)
      << "Invalid initial state: probe = "
      << probe << ", hit_nucleon = " << hit_nucleon;
    return false;
  }

  // assert that the interaction mode is allowed
  bool isu  = pdg::IsUQuark     (hit_quark);
  bool isd  = pdg::IsDQuark     (hit_quark);
  bool iss  = pdg::IsSQuark     (hit_quark);
  bool isub = pdg::IsAntiUQuark (hit_quark);
  bool isdb = pdg::IsAntiDQuark (hit_quark);
  bool issb = pdg::IsAntiSQuark (hit_quark);

  bool allowed = (iscc && isv  && (isd||isub||iss))  ||
                 (iscc && isvb && (isu||isdb||issb)) ||
                 (isnc && (isv||isvb) && (isu||isd||isub||isdb||iss||issb)) ||
                 (isdmi && isdm && (isu||isd||isub||isdb||iss||issb)) ||
                 (isem && (isl||islb) && (isu||isd||isub||isdb||iss||issb));
  if(!allowed) {
    LOG("Hadronization", pWARN)
      << "Impossible interaction type / probe / hit quark combination!";
    return false;
  }

  return true;
}
//____________________________________________________________________________
void Hadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigRemoveMe();

  fAllowReconfig = false;
}
//___________________________________________________________________________
void Hadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigRemoveMe();

  fAllowReconfig = false;
}
//___________________________________________________________________________
void Hadronization::LoadConfigRemoveMe(void)
{
  // Check whether to generate weighted or unweighted particle decays
  fGenerateWeighted = false ;
  //this->GetParam("GenerateWeighted", fGenerateWeighted, false);{

  // decayer
  fDecayer = 0;
  if( GetConfig().Exists("Decayer") ) {
     fDecayer = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("Decayer"));
     assert(fDecayer);
  }

  // Load Wcut determining the phase space area where the multiplicity prob.
  // scaling factors would be applied -if requested-
  this->GetParam( "Wcut", fWcut ) ;

  // Load NEUGEN multiplicity probability scaling parameters Rijk
  // neutrinos
  this->GetParam( "DIS-HMultWgt-vp-CC-m2",  fRvpCCm2  ) ;
  this->GetParam( "DIS-HMultWgt-vp-CC-m3",  fRvpCCm3  ) ;
  this->GetParam( "DIS-HMultWgt-vp-NC-m2",  fRvpNCm2  ) ;
  this->GetParam( "DIS-HMultWgt-vp-NC-m3",  fRvpNCm3  ) ;
  this->GetParam( "DIS-HMultWgt-vn-CC-m2",  fRvnCCm2  ) ;
  this->GetParam( "DIS-HMultWgt-vn-CC-m3",  fRvnCCm3  ) ;
  this->GetParam( "DIS-HMultWgt-vn-NC-m2",  fRvnNCm2  ) ;
  this->GetParam( "DIS-HMultWgt-vn-NC-m3",  fRvnNCm3  ) ;
  //Anti-neutrinos
  this->GetParam( "DIS-HMultWgt-vbp-CC-m2", fRvbpCCm2 ) ;
  this->GetParam( "DIS-HMultWgt-vbp-CC-m3", fRvbpCCm3 ) ;
  this->GetParam( "DIS-HMultWgt-vbp-NC-m2", fRvbpNCm2 ) ;
  this->GetParam( "DIS-HMultWgt-vbp-NC-m3", fRvbpNCm3 ) ;
  this->GetParam( "DIS-HMultWgt-vbn-CC-m2", fRvbnCCm2 ) ;
  this->GetParam( "DIS-HMultWgt-vbn-CC-m3", fRvbnCCm3 ) ;
  this->GetParam( "DIS-HMultWgt-vbn-NC-m2", fRvbnNCm2 ) ;
  this->GetParam( "DIS-HMultWgt-vbn-NC-m3", fRvbnNCm3 ) ;
}
//___________________________________________________________________________
