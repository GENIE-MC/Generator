//____________________________________________________________________________
/*!

\namespace  genie::utils::nuclear

\brief      Simple nuclear physics empirical formulas

\author     Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
            CCLRC, Rutherford Appleton Laboratory

\created    May 06, 2004

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "PDG/PDGUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
double genie::utils::nuclear::BindEnergy(const Target & target)
{
// Compute the average binding energy (in GeV) using the semi-empirical
// formula from Wapstra (Handbuch der Physik, XXXVIII/1)

  if(!target.IsNucleus()) return 0;

  double a = 15.835;
  double b =  18.33;
  double s =  23.20;
  double d =   0.714;

  double delta = 0;                       /*E-O*/
  if (target.IsOddOdd()  ) delta =  11.2; /*O-O*/
  if (target.IsEvenEven()) delta = -11.2; /*E-E*/

  double N = (double) target.N();
  double Z = (double) target.Z();
  double A = (double) target.A();

  double BE =  a * A
             - b * TMath::Power(A,0.667)
             - s * TMath::Power(N-Z,2.0)/A
             - d * TMath::Power(Z,2.0)/TMath::Power(A,0.333)
             - delta / TMath::Sqrt(A); // MeV

  return (1e-3 * BE); // GeV
}
//___________________________________________________________________________
double genie::utils::nuclear::BindEnergyPerNucleon(const Target & target)
{
// Compute the average binding energy per nucleon (in GeV)

  if(!target.IsNucleus()) return 0;

  return (utils::nuclear::BindEnergy(target) / target.A());
}
//___________________________________________________________________________
double genie::utils::nuclear::BindEnergyLastNucleon(const Target & target)
{
// Compute the binding for the most loose nucleon (in GeV)

  if(!target.IsNucleus()) return 0;

  //-- temporarily, return the binding energy per nucleon rather than the
  //   separation energy of the last nucleon
   return (utils::nuclear::BindEnergy(target) / target.A());
}
//___________________________________________________________________________
double genie::utils::nuclear::Radius(const Target & target)
{
// Compute the nuclear radius (in GeV^-1)
// (to get it in F, for example, use: Radius(target)/genie::units::fermi)

  return Radius(target.A());
}
//___________________________________________________________________________
double genie::utils::nuclear::Radius(int A)
{
// Compute the nuclear radius (in GeV^-1)
// (to get it in F, for example, use: Radius(target)/genie::units::fermi)

  if(A<1) return 0;

  double Rn = kNucRo * TMath::Power(A, 0.3333333); // in GeV^-1
  return Rn;
}
//___________________________________________________________________________
double genie::utils::nuclear::NuclQELXSecSuppression(
                string kftable, double pmax, const Interaction * interaction)
{
// Computes the suppression of the EL/QEL differential cross section due to
// nuclear effects.
//
// A direct adaptation of NeuGEN's qelnuc() and r_factor()
//
// Hugh's comments from the original code: 
// "This routine is based on an analytic calculation of the rejection factor 
//  in the Fermi Gas model using the form for the fermi momentum distribution 
//  given in the  Bodek and Ritchie paper.  [P.R.  D23 (1981) 1070]
//  R is the ratio of the differential cross section from the nuclear material 
//  specified by (kf,pf) to the differential cross section for a free nucleon".
//  (kf,pf = Fermi Gas model Fermi momentum for initial,final nucleons)
//
  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();
  const Target &       target     = init_state.GetTarget();

  if (!target.IsNucleus()) return 1.0; // no suppression for free nucleon tgt

  double R = 1; // computed suppression factor

  int target_pdgc         = target.PDGCode();
  int struck_nucleon_pdgc = target.StruckNucleonPDGCode();
  int final_nucleon_pdgc  = struck_nucleon_pdgc;

  if(proc_info.IsWeakCC()) 
          final_nucleon_pdgc = pdg::SwitchProtonNeutron(struck_nucleon_pdgc);

  // get the requested Fermi momentum table 
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable(kftable);

  // Fermi Gas model Fermi momentum for initial,final nucleons

  double kFi = kft->FindClosestKF(target_pdgc, struck_nucleon_pdgc);
  double kFf = (struck_nucleon_pdgc==final_nucleon_pdgc) ? kFi : 
               kft->FindClosestKF(target_pdgc, final_nucleon_pdgc );

  // Compute magnitude of the 3-momentum transfer to the nucleon

  double Mn2 = target.StruckNucleonP4()->M2(); // can be off m/shell

  const Kinematics & kine = interaction->GetKinematics();
  double q2     = kine.q2();
  double magq2  = q2 * (0.25*q2/Mn2 - 1.);
  double q      = TMath::Sqrt(TMath::Max(0.,magq2));

  double kfa   = kFi * 2./kPi;
  double kfa2  = TMath::Power(kfa,2);
  double kFi4  = TMath::Power(kFi,4);
  double rkf   = 1./(1. - kFi/4.);
  double alpha = 1. - 6.*kfa2;
  double beta  = 2. * rkf * kfa2 * kFi4;

  double fm_area = FmArea(alpha,beta,kFi,pmax);

  if (q <= kFf) {

     double p1   = kFf - q;
     double p2   = kFf + q;
     double fmi1 = FmI1(alpha,beta,p1,p2,  kFi,kFf,q);
     double fmi2 = FmI2(alpha,beta,p2,pmax,kFi,kFf,q);

     R = 2*kPi * (fmi1 + 2*fmi2) / fm_area;

  } else if (q > kFf && q <= (pmax-kFf)) {

     double p1    = q - kFf;
     double p2    = q + kFf;
     double fmi1  = FmI1(alpha,beta, p1, p2,   kFi, kFf,q);
     double fmi2a = FmI2(alpha,beta, 0., p1,   kFi, kFf,q);
     double fmi2b = FmI2(alpha,beta, p2, pmax, kFi, kFf,q);

     R = 2*kPi * (fmi1 + 2*(fmi2a+fmi2b)) / fm_area;

  } else if (q > (pmax-kFf) && q <= (pmax+kFf)) {

     double p1   = q - kFf;
     double fmi2 = FmI2(alpha,beta, 0.,p1,  kFi,kFf,q);
     double fmi1 = FmI1(alpha,beta, p1,pmax,kFi,kFf,q);

     R = 2*kPi * (2.*fmi2 + fmi1) / fm_area;

  } else if (q > (pmax+kFf)) {
     R = 1.; 

  } else {
     LOG("Nuclear", pFATAL) << "Illegal input q = " << q;
     exit(1);
  }
  return R;
}
//___________________________________________________________________________
double genie::utils::nuclear::FmI1(double alpha, double beta, 
                       double a, double b,  double kFi, double kFf, double q)
{
// Adapted from NeuGEN's fm_integral1() used in r_factor()

  double f=0;

  double q2   = TMath::Power(q,  2);
  double a2   = TMath::Power(a,  2);
  double b2   = TMath::Power(b,  2);
  double kFi2 = TMath::Power(kFi,2);
  double kFf2 = TMath::Power(kFf,2);

  if(kFi < a) {
     double lg = TMath::Log(b/a);

     f = -beta * (kFf2-q2)/(4.*q) * (1./a2 - 1./b2) + beta/(2.*q) * lg;

  } else if (kFi > b) {
     double a4 = TMath::Power(a2,2);
     double b4 = TMath::Power(b2,2);

     f = - (kFf2-q2) * alpha/(4.*q) * (b2-a2) + alpha/(8.*q) * (b4-a4);

  } else {
     double a4   = TMath::Power(a2,2);
     double kFi4 = TMath::Power(kFi2,2);
     double lg   = TMath::Log(b/kFi);

     f = - alpha*(kFf2-q2)/(4.*q)*(kFi2-a2) + alpha/(8.*q)*(kFi4 - a4)
         - beta*(kFf2-q2)/(4.*q)*(1./kFi2 - 1./b2) + beta/(2.*q)*lg;
  }

  double integral2 = FmI2(alpha,beta,a,b,kFi,kFf,q);
  double integral1 = integral2 + f;
  
  return integral1;
}
//___________________________________________________________________________
double genie::utils::nuclear::FmI2(double alpha, double beta, 
                       double a, double b,  double kFi, double kFf, double q)
{
// Adapted from NeuGEN's fm_integral2() used in r_factor()

  double integral2 = 0;

  if(kFi < a) {
     integral2 = beta * (1./a - 1./b);

  } else if(kFi > b) {
     double a3 = TMath::Power(a,3);
     double b3 = TMath::Power(b,3);
     integral2 = alpha/3. * (b3 - a3);

  } else {
     double a3   = TMath::Power(a,3);
     double kFi3 = TMath::Power(kFi,3);

     integral2 = alpha/3. * (kFi3 - a3) + beta * (1./kFi - 1./b);
  }
  return integral2;
}
//___________________________________________________________________________
double genie::utils::nuclear::FmArea(
                         double alpha, double beta, double kf, double pmax)
{
// Adapted from NeuGEN's fm_area() used in r_factor()

  double kf3 = TMath::Power(kf,3.);
  double sum = 4.*kPi* (alpha * kf3/3. + beta*(1./kf - 1./pmax));
  return sum;
}
//___________________________________________________________________________


