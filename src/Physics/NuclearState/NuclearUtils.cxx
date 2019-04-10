//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

 @ Sep 19, 2007 - CA
   Copied here the nuclear density methods used by intranuke so that they
   can be used by other modules as well (eg position vertex selection)
 @ Nov 30, 2007 - CA
   Added possibility to increase the nuclear size (intranuke may increase the
   nuclear size by a const times the de Broglie wavelength of a transported
   hadron). Density() methods have a new default argument (ring) which is 0
   if not explicity set.
 @ Nov 30, 2009 - CA
   Overriding NuclQELXSecSuppression() calc for deuterium and using input 
   data from S.K.Singh, Nucl. Phys. B 36, 419 (1972) [data used by Hugh for
   the Merenyi test]
 @ May 18, 2010 - CA
   Restructure NuclQELXSecSuppression() (Add RQEFG_generic()) to aid 
   reweighting. 
 @ Jul 15, 2010 - AM
   Added BindEnergy(int nucA, int nucZ), used in Intranuke
 @ Mar 18, 2016 - JJ (SD)
   Check if a local Fermi gas model should be used when calculating the
   Fermi momentum
*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearData.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/NuclearModelI.h"

using namespace genie;
using namespace genie::constants;


//____________________________________________________________________________
double genie::utils::nuclear::BindEnergy(const Target & target)
{
// Compute the average binding energy (in GeV) using the semi-empirical
// formula from Wapstra (Handbuch der Physik, XXXVIII/1)

  if(!target.IsNucleus()) return 0;

  int Z = (int) target.Z();
  int A = (int) target.A();

  return BindEnergy(A,Z);
}
//___________________________________________________________________________
double genie::utils::nuclear::BindEnergy(int nucA, int nucZ)
{
// Compute the average binding energy (in GeV) using the semi-empirical
// formula from Wapstra (Handbuch der Physik, XXXVIII/1)

  if (nucA<=0 || nucZ<=0) return 0;

  double a =  15.835;
  double b =  18.33;
  double s =  23.20;
  double d =   0.714;

  double delta = 0;                          /*E-O*/
  if (nucZ%2==1 && nucA%2==1) delta =  11.2; /*O-O*/
  if (nucZ%2==0 && nucA%2==0) delta = -11.2; /*E-E*/

  double N = (double) (nucA-nucZ);
  double Z = (double) nucZ;
  double A = (double) nucA;

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
double genie::utils::nuclear::Radius(int A, double Ro)
{
// Compute the nuclear radius (in fm)

  if(A<1) return 0;

  double Rn = Ro*TMath::Power(1.*A, 0.3333333);
  return Rn;
}
//___________________________________________________________________________
double genie::utils::nuclear::NuclQELXSecSuppression(
                string kftable, double pmax, const Interaction * interaction)
{
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  if (!target.IsNucleus()) {
     return 1.0; // no suppression for free nucleon tgt
  }

  //
  // special case: deuterium   
  //

  if (target.A() == 2) {
    NuclearData * nucldat = NuclearData::Instance();
    return nucldat->DeuteriumSuppressionFactor(interaction->Kine().Q2());
  }

  //
  // general case
  //

  int target_pdgc         = target.Pdg();
  int struck_nucleon_pdgc = target.HitNucPdg();
  int final_nucleon_pdgc  = struck_nucleon_pdgc;

  if(proc_info.IsWeakCC()) {
     final_nucleon_pdgc = pdg::SwitchProtonNeutron(struck_nucleon_pdgc);
  }

  // Check if an LFG model should be used for Fermi momentum
  // Create a nuclear model object to check the model type
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  RgKey nuclkey = "NuclearModel";
  RgAlg nuclalg = gc->GetAlg(nuclkey);
  AlgFactory * algf = AlgFactory::Instance();
  const genie::NuclearModelI* nuclModel =
    dynamic_cast<const genie::NuclearModelI*>(
			     algf->GetAlgorithm(nuclalg.name,nuclalg.config));
  // Check if the model is a local Fermi gas
  bool lfg = (nuclModel && nuclModel->ModelType(Target()) == kNucmLocalFermiGas);

  double kFi, kFf;
  if(lfg){
    double hbarc = kLightSpeed*kPlankConstant/genie::units::fermi;
    Target* tgt = interaction->InitStatePtr()->TgtPtr();
    double radius = tgt->HitNucPosition();

    int A = tgt->A();
    // kFi
    bool is_p_i = pdg::IsProton(struck_nucleon_pdgc);
    double numNuci = (is_p_i) ? (double)tgt->Z():(double)tgt->N();
    kFi = TMath::Power(3*kPi2*numNuci*
		     genie::utils::nuclear::Density(radius,A),1.0/3.0) *hbarc;
    // kFi
    bool is_p_f = pdg::IsProton(final_nucleon_pdgc);
    double numNucf = (is_p_f) ? (double)tgt->Z():(double)tgt->N();
    kFf = TMath::Power(3*kPi2*numNucf*
		     genie::utils::nuclear::Density(radius,A),1.0/3.0) *hbarc;
  }else{
    // get the requested Fermi momentum table 
    FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
    const FermiMomentumTable * kft  = kftp->GetTable(kftable);
    
    // Fermi momentum for initial, final nucleons
    kFi = kft->FindClosestKF(target_pdgc, struck_nucleon_pdgc);
    kFf = (struck_nucleon_pdgc==final_nucleon_pdgc) ? kFi : 
          kft->FindClosestKF(target_pdgc, final_nucleon_pdgc );
  }
  
  double Mn = target.HitNucP4Ptr()->M(); // can be off m/shell

  const Kinematics & kine = interaction->Kine();
  double q2 = kine.q2();

  double R = RQEFG_generic(q2, Mn, kFi, kFf, pmax);
  return R;
}
//___________________________________________________________________________
double genie::utils::nuclear::RQEFG_generic (
            double q2, double Mn, double kFi, double kFf, double pmax)
{
// Computes the nuclear suppression of the QEL differential cross section 
//
// Inputs:
//  - q2   : momentum transfer (< 0)
//  - Mn   : hit nucleon mass (nucleon can be off the mass shell)
//  - kFi  : Fermi momentum, initial state (hit) nucleon @ nuclear target
//  - kFf  : Fermi momentum, final state (recoil) nucleon @ nuclear target
//  - pmax : A Fermi momentum cutoff
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

  double R = 1.;

  // Compute magnitude of the 3-momentum transfer to the nucleon
  double Mn2    = Mn*Mn;
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
               double a, double b,  double kFi, double /*kFf*/, double /*q*/)
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
double genie::utils::nuclear::DISNuclFactor(double x, int A)
{
// Adapted from NeuGEN's nuc_factor(). Kept original comments from Hugh.

  double xv     = TMath::Min(0.75, x);
  double xv2    = xv  * xv;
  double xv3    = xv2 * xv;
  double xv4    = xv3 * xv;
  double xv5    = xv4 * xv;
  double xvp    = TMath::Power(xv, 14.417);
  double expaxv = TMath::Exp(-21.94*xv);

  double f = 1.;

  // first factor goes from free nucleons to deuterium
  if(A >= 2) {
    f= 0.985*(1.+0.422*xv - 2.745*xv2 + 7.570*xv3 - 10.335*xv4 + 5.422*xv5);
  }
  // 2nd factor goes from deuterium to iso-scalar iron
  if(A > 2) {
    f *= (1.096 - 0.364*xv - 0.278*expaxv + 2.722*xvp);
  }
  return f;
}
//___________________________________________________________________________
double genie::utils::nuclear::Density(double r, int A, double ring)
{
// [by S.Dytman]
//
  if(A>20) {
    double c = 1., z = 1.;

    if      (A ==  27) { c = 3.07; z = 0.52; }  // aluminum
    else if (A ==  28) { c = 3.07; z = 0.54; }  // silicon
    else if (A ==  40) { c = 3.53; z = 0.54; }  // argon
    else if (A ==  56) { c = 4.10; z = 0.56; }  // iron
    else if (A == 208) { c = 6.62; z = 0.55; }  // lead
    else { 
       c = TMath::Power(A,0.35); z = 0.54; 
    } //others

    LOG("Nuclear",pINFO)
	<< "r= " << r << ", ring= " << ring;
    double rho = DensityWoodsSaxon(r,c,z,ring);
    return rho;
  }
  else if (A>4) {
    double ap = 1., alf = 1.;

    if      (A ==  7) { ap = 1.77; alf = 0.327; } // lithium
    else if (A == 12) { ap = 1.69; alf = 1.08;  } // carbon
    else if (A == 14) { ap = 1.76; alf = 1.23;  } // nitrogen
    else if (A == 16) { ap = 1.83; alf = 1.54;  } // oxygen
    else  { 
      ap=1.75; alf=-0.4+.12*A; 
    }  //others- alf=0.08 if A=4

    double rho = DensityGaus(r,ap,alf,ring);
    return rho;
  }
  else {
    // helium
    double ap = 1.9/TMath::Sqrt(2.);  
    double alf=0.;    
    double rho = DensityGaus(r,ap,alf,ring);
    return rho;
  }

  return 0;
}
//___________________________________________________________________________
double genie::utils::nuclear::DensityGaus(
                             double r, double a, double alf, double ring)
{
// [adapted from neugen3 density_gaus.F written by S.Dytman]
//
// Modified harmonic osc density distribution. 
// Norm gives normalization to 1
//
// input  : radial distance in nucleus [units: fm]
// output : nuclear density            [units: fm^-3]

  ring = TMath::Min(ring, 0.3*a);

  double aeval = a + ring;
  double norm  = 1./((5.568 + alf*8.353)*TMath::Power(a,3.));  //0.0132;
  double b     = TMath::Power(r/aeval, 2.);
  double dens  = norm * (1. + alf*b) * TMath::Exp(-b);

  LOG("Nuclear", pINFO) 
        << "r = " << r << ", norm = " << norm << ", dens = " << dens 
        << ", aeval= " << aeval;

  return dens;
}
//___________________________________________________________________________
double genie::utils::nuclear::DensityWoodsSaxon(
                              double r, double c, double z, double ring)
{
// [adapted from neugen3 density_ws.F written by S.Dytman]
//
// Woods-Saxon desity distribution. Norn gives normalization to 1
//
// input  : radial distance in nucleus [units: fm]
// output : nuclear density            [units: fm^-3]

  LOG("Nuclear",pINFO)
      << "c= " << c << ", z= " << z << ",ring= " << ring;

  ring = TMath::Min(ring, 0.75*c);

  double ceval = c + ring;
  double norm  = (3./(4.*kPi*TMath::Power(c,3)))*1./(1.+TMath::Power((kPi*z/c),2));
  double dens  = norm / (1 + TMath::Exp((r-ceval)/z));

  LOG("Nuclear", pINFO) 
     << "r = " << r << ", norm = " << norm 
     << ", dens = " << dens << " , ceval= " << ceval;

  return dens;
}
//___________________________________________________________________________
double genie::utils::nuclear::BindEnergyPerNucleonParametrization(const Target & target)
{
// Compute the average binding energy per nucleon (in GeV)
if(!target.IsNucleus()) return 0;
double x = TMath::Power(target.A(),1/3.0) / target.Z();
return (0.05042772591+x*(-0.11377355795+x*(0.15159890400-0.08825307197*x)));
}
//___________________________________________________________________________
double genie::utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(const Target & target)
{
// Compute Fermi momentum for isoscalar nucleon (in GeV)
if(!target.IsNucleus()) return 0;
double x = 1.0 / target.A();
return (0.27+x*(-1.12887857491+x*(9.72670908033-39.53095724456*x)));
}
//___________________________________________________________________________

