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
#include "Nuclear/BodekRitchie.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
BodekRitchie::BodekRitchie() :
NuclearPDistributionModelI()
{
  fName = "genie::BodekRitchie";
}
//____________________________________________________________________________
BodekRitchie::BodekRitchie(const char * param_set) :
NuclearPDistributionModelI(param_set)
{
  fName = "genie::BodekRitchie";

  FindConfig();
}
//____________________________________________________________________________
BodekRitchie::~BodekRitchie()
{

}
//____________________________________________________________________________
TH1D * BodekRitchie::ProbabilityDistribution(const Target & target) const
{
  //-- look up the Fermi momentum for this Target
  
  double KF = 0.225; // GeV - just use a constant value for now...

  //-- correct the Fermi momentum for the struck nucleon

  int nucleon_pdgc = target.StruckNucleonPDGCode();
  
  assert( pdg::IsProton(nucleon_pdgc) || pdg::IsNeutron(nucleon_pdgc) );

  double Z = (double) target.Z();
  double N = (double) target.N();
  double A = (double) target.A();
  
  if( pdg::IsProton(nucleon_pdgc) )
     KF *= pow( 2 * Z/A, 1./3.);
  else
     KF *= pow( 2 * N/A, 1./3.);

  double a  = 2.0;
  double C  = 4. * kPi * pow(KF,3) / 3.;
  double R  = 1. / ( 1. - KF/4.);

  LOG("BodekRitchie", pDEBUG) << "KF = " << KF;
  LOG("BodekRitchie", pDEBUG) << "a  = " << a;
  LOG("BodekRitchie", pDEBUG) << "C  = " << C;
  LOG("BodekRitchie", pDEBUG) << "R  = " << R;
        
  //-- create the probability distribution

  int    n = (fConfig->Exists("n-momentum-bins")) ?
                                 fConfig->GetInt("n-momentum-bins") : 100;
  double pmax = (fConfig->Exists("p-max")) ?
                                         fConfig->GetDouble("p-max") : 10;  
  double pcutoff = (fConfig->Exists("p-cut-off")) ?
                                      fConfig->GetDouble("p-cut-off") : 4;

  LOG("BodekRitchie", pDEBUG) << "P-cut-off = " << pcutoff;

  assert(n > 1 && pmax > 0 && pcutoff > 0 && pcutoff < pmax);
                                        
  TH1D * prob = new TH1D("", "", n, 0, pmax);

  double dp = pmax / (n-1);
  
  for(int i = 0; i < n; i++) {

     double p = i * dp;

     //-- calculate |phi(p)|^2
     
     double phi2 = 0;

     if (p <= KF)
     {
        phi2 = (1./C) * (1. - 6.*pow(KF*a/kPi,2));
     }
     else if ( p > KF && p < pcutoff)
     {
        phi2 = (1./C) * (2*R*pow(KF*a/kPi,2.)*pow(KF/p,4.));
     }
     
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
