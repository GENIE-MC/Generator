//____________________________________________________________________________
/*!

\class    genie::BodekRitchie

\brief    The Bodek-Ritchie model for the probability distribution of nucleon
          momenta within a nucleus.
          Implements the NuclMomentumModelI interface.

\ref      Bodek and Ritchie, Phys.Rev.D23:1070 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Messenger/Messenger.h"
#include "Nuclear/BodekRitchie.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BodekRitchie::BodekRitchie() :
NuclMomentumModelI("genie::BodekRitchie")
{

}
//____________________________________________________________________________
BodekRitchie::BodekRitchie(string config) :
NuclMomentumModelI("genie::BodekRitchie", config)
{

}
//____________________________________________________________________________
BodekRitchie::~BodekRitchie()
{

}
//____________________________________________________________________________
TH1D * BodekRitchie::ProbabilityDistribution(const Target & target) const
{
  LOG("BodekRitchie", pNOTICE)
             << "Computing P = f(p_nucleon) for: " << target.AsString();
  LOG("BodekRitchie", pNOTICE)
               << "P(cut-off) = " << fPCutOff << ", P(max) = " << fPMax;

  //-- get information for the nuclear target
  int target_pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  int nucleon_pdgc = target.HitNucPdg();
  assert( pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc) );
  double Z = (double) target.Z();
  double N = (double) target.N();
  double A = (double) target.A();

  //-- look up the Fermi momentum for this Target
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable("Default");
  double KF = kft->FindClosestKF(target_pdgc, nucleon_pdgc);
  LOG("BodekRitchie", pNOTICE) << "KF = " << KF;

  //-- correct the Fermi momentum for the struck nucleon
  assert(target.HitNucIsSet());
  bool is_p = pdg::IsProton(nucleon_pdgc);
  if(is_p) KF *= TMath::Power( 2*Z/A, 1./3.);
  else     KF *= TMath::Power( 2*N/A, 1./3.);
  LOG("BodekRitchie", pINFO) << "Corrected KF = " << KF;

  double a  = 2.0;
  double C  = 4. * kPi * TMath::Power(KF,3) / 3.;
  double R  = 1. / (1.- KF/4.);
  LOG("BodekRitchie", pDEBUG) << "a  = " << a;
  LOG("BodekRitchie", pDEBUG) << "C  = " << C;
  LOG("BodekRitchie", pDEBUG) << "R  = " << R;

  //-- create the probability distribution

  TH1D * prob = new TH1D("", "", fNPBins, 0, fPMax);

  double dp = fPMax / (fNPBins-1);
  double iC = (C>0) ? 1./C : 0.;
  double kfa_pi_2 = TMath::Power(KF*a/kPi,2);

  for(int i = 0; i < fNPBins; i++) {
     double p  = i * dp;
     double p2 = TMath::Power(p,2);

     // calculate |phi(p)|^2
     double phi2 = 0;
     if (p <= KF)
        phi2 = iC * (1. - 6.*kfa_pi_2);
     else if ( p > KF && p < fPCutOff)
        phi2 = iC * (2*R*kfa_pi_2*TMath::Power(KF/p,4.));

     // calculate probability density : dProbability/dp
     double dP_dp = 4*kPi * p2 * phi2;
     LOG("BodekRitchie", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
     prob->Fill(p, dP_dp);
  }

  //-- normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  return prob; // Note: The calling function adopts the object
}
//____________________________________________________________________________
void BodekRitchie::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//____________________________________________________________________________
void BodekRitchie::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
}
//____________________________________________________________________________
void BodekRitchie::LoadConfigData(void)
{
  fNPBins  = fConfig->GetIntDef("n-momentum-bins", 250);
  fPMax    = fConfig->GetDoubleDef("p-max",     10.0);
  fPCutOff = fConfig->GetDoubleDef("p-cut-off",  0.5);

  assert(fNPBins > 1 && fPMax > 0 && fPCutOff > 0 && fPCutOff < fPMax);
}
//____________________________________________________________________________

