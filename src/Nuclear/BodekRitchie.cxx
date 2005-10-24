//____________________________________________________________________________
/*!

\class    genie::BodekRitchie

\brief    The Bodek-Ritchie model for the probability distribution of nucleon
          momenta within a nucleus.

          Implements the NuclearPDistributionModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

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
NuclearPDistributionModelI("genie::BodekRitchie")
{

}
//____________________________________________________________________________
BodekRitchie::BodekRitchie(string config) :
NuclearPDistributionModelI("genie::BodekRitchie", config)
{

}
//____________________________________________________________________________
BodekRitchie::~BodekRitchie()
{

}
//____________________________________________________________________________
TH1D * BodekRitchie::ProbabilityDistribution(const Target & target) const
{
  //-- get information for the nuclear target
  int target_pdgc  = pdg::IonPdgCode(target.A(), target.Z());
  int nucleon_pdgc = target.StruckNucleonPDGCode();
  assert( pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc) );
  double Z = (double) target.Z();
  double N = (double) target.N();
  double A = (double) target.A();

  //-- look up the Fermi momentum for this Target
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable("Default");
  double KF = kft->FindClosestKF(target_pdgc, nucleon_pdgc);
  LOG("BodekRitchie", pDEBUG) << "KF = " << KF;

  //-- correct the Fermi momentum for the struck nucleon
  bool is_p = pdg::IsProton(nucleon_pdgc);
  if(is_p) KF *= pow( 2 * Z/A, 1./3.);
  else     KF *= pow( 2 * N/A, 1./3.);
  LOG("BodekRitchie", pDEBUG) << "Corrected KF = " << KF;

  double a  = 2.0;
  double C  = 4. * kPi * pow(KF,3) / 3.;
  double R  = 1. / ( 1. - KF/4.);
  LOG("BodekRitchie", pDEBUG) << "a  = " << a;
  LOG("BodekRitchie", pDEBUG) << "C  = " << C;
  LOG("BodekRitchie", pDEBUG) << "R  = " << R;

  //-- create the probability distribution
  int    n       = fConfig->GetIntDef("n-momentum-bins", 250);
  double pmax    = fConfig->GetDoubleDef("p-max", 10);
  double pcutoff = fConfig->GetDoubleDef("p-cut-off", 4);

  LOG("BodekRitchie", pDEBUG)
      << "P(cut-off) = " << pcutoff << ", P(max) = " << pmax;

  assert(n > 1 && pmax > 0 && pcutoff > 0 && pcutoff < pmax);

  TH1D * prob = new TH1D("", "", n, 0, pmax);

  double dp = pmax / (n-1);
  for(int i = 0; i < n; i++) {
     double p = i * dp;

     //-- calculate |phi(p)|^2
     double phi2 = 0;
     if (p <= KF)
        phi2 = (1./C) * (1. - 6.*pow(KF*a/kPi,2));
     else if ( p > KF && p < pcutoff)
        phi2 = (1./C) * (2*R*pow(KF*a/kPi,2.)*pow(KF/p,4.));

     //-- calculate probability density : dProbability/dp
     double dP_dp = 4 * kPi * p*p * phi2;
     LOG("BodekRitchie", pDEBUG) << "p = " << p << ", dP/dp = " << dP_dp;
     prob->Fill(p, dP_dp);
  }

  //----- Normalize the probability distribution
  prob->Scale( 1.0 / prob->Integral("width") );

  return prob; // Note: The calling function adopts the object
}
//____________________________________________________________________________
