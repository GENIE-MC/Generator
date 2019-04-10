//____________________________________________________________________________
/*
   Copyright (c) 2003-2019, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
   or see $GENIE/LICENSE

Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
STFC, Rutherford Appleton Laboratory - March 11, 2005

For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Coherent/XSection/BergerSehgalCOHPiPXSec2015.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
BergerSehgalCOHPiPXSec2015::BergerSehgalCOHPiPXSec2015() :
  XSecAlgorithmI("genie::BergerSehgalCOHPiPXSec2015")
{

}
//____________________________________________________________________________
BergerSehgalCOHPiPXSec2015::BergerSehgalCOHPiPXSec2015(string config) :
  XSecAlgorithmI("genie::BergerSehgalCOHPiPXSec2015", config)
{

}
//____________________________________________________________________________
BergerSehgalCOHPiPXSec2015::~BergerSehgalCOHPiPXSec2015()
{

}
//____________________________________________________________________________
double BergerSehgalCOHPiPXSec2015::XSec(
    const Interaction * interaction, KinePhaseSpace_t kps) const
{
  // Here we are following PRD 79, 053003 (2009) by Berger and Sehgal
  // This method computes the differential cross section represented 
  // in Eq.'s 6 (CC) and 7 (NC) from that paper.

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  bool pionIsCharged = interaction->ProcInfo().IsWeakCC();
  double xsec = 0.0;

  double E      = init_state.ProbeE(kRfLab);        // nu E
  double Q2     = kinematics.Q2();
  double y      = kinematics.y();                   // inelasticity
  double x      = kinematics.x();          
  assert(E > 0.);
  assert(y > 0.);
  assert(y < 1.);
  double ppistar = PionCOMAbsMomentum(interaction); // |Center of Mass Momentum|
  if (ppistar <= 0.0) { 
    LOG("BergerSehgalCohPi", pDEBUG) << "Pion COM momentum negative for Q2 = " << Q2 << 
      " x = " << x << " y = " << y; 
    return 0.0; 
  }
  double front  = ExactKinematicTerm(interaction);
  if (front <= 0.0) { 
    LOG("BergerSehgalCohPi", pDEBUG) << "Exact kin. form = " << front << 
      " E = " << E << " Q2 = " << Q2 << " y = " << y << " x = " << x; 
    return 0.0; 
  }

  double A      = (double) init_state.Tgt().A();   // mass number
  double A2     = TMath::Power(A, 2.);
  double A_3    = TMath::Power(A, 1./3.);
  double M      = init_state.Tgt().Mass();
  double M_pi   = pionIsCharged ? kPionMass : kPi0Mass;
  double Epi    = y*E;                             // ~pion energy
  double ma2    = TMath::Power(fMa, 2);            // "axial mass" squared
  double Ga     = ma2 / (ma2 + Q2);
  double Ga2    = TMath::Power(Ga, 2.);            // propagator term
  double Ro2    = TMath::Power(fRo * units::fermi, 2.);

  // the xsec is d^3xsec/dQ^2dydt but the only t-dependent factor 
  // is an exp(-bt) so it can be integrated analyticaly
  double Epi2   = TMath::Power(Epi, 2.);
  double R      = fRo * A_3 * units::fermi; // nuclear radius
  double R2     = TMath::Power(R, 2.);
  double b      = 0.33333 * R2;
  double MxEpi  = M * x / Epi;
  double mEpi2  = (M_pi * M_pi) / Epi2;
  double tA     = 1. + MxEpi - 0.5 * mEpi2;
  double tB     = TMath::Sqrt(1.0 + 2 * MxEpi) * TMath::Sqrt(1.0 - mEpi2);
  double tmin   = 2 * Epi2 * (tA - tB);
  double tmax   = 2 * Epi2 * (tA + tB);
  if (tmin < 1.0e-8) {
      tmin = 1.0e-8;
  }

  /* const KPhaseSpace & kphase = interaction->PhaseSpace(); */
  /* Range1D_t tl = kphase.TLim();   // TESTING! */

  double sigtot_pin  = utils::hadxs::berger::PionNucleonXSec(Epi, /* get_total = */ true, pionIsCharged);
  double sigel_pin   = utils::hadxs::berger::PionNucleonXSec(Epi, /* get_total = */ false, pionIsCharged);
  double siginel_pin = sigtot_pin - sigel_pin;

  // fabs (F_{abs}) describes the average attenuation of a pion emerging
  // from a sphere of nuclear matter with radius = R_0 A^{1/3}. it is 
  // Eq. 13 in Berger-Sehgal PRD 79, 053003
  double fabs_input  = (9.0 * A_3) / (16.0 * kPi * Ro2);
  double fabs        = TMath::Exp( -1.0 * fabs_input * siginel_pin);

  // my old hackery to get things to work, A. Mislivec provided a better alt.
  // double factor      = 0.1; // to go from 10^-37 cm^2 -> 10^-38 cm^2
  // double RS_factor   = (units::mb*A2*fabs)/(16.0*kPi) * (sigtot_pin*sigtot_pin);

  // A_RS for BS version of RS, and/or Tpi>1.0
  double RS_factor = (A2 * fabs) / (16.0 * kPi) * (sigtot_pin * sigtot_pin);

  // get the pion-nucleus cross section on carbon, fold it into differential cross section
  double tpi         = (E * y) - M_pi - ((Q2 + M_pi * M_pi) / (2 * M)); 
  double tpilow      = 0.0;
  double siglow      = 0.0;
  double tpihigh     = 0.0;
  double sighigh     = 0.0;
  double dsigdzfit   = 0.0;
  double dsigdtfit   = 0.0;
  int    xsec_stat   = 0;
  double dsig        = 0.0;
  double tstep       = 100;
  double logt_step   = TMath::Abs(log(tmax) - log(tmin)) / tstep;
  double logt        = log(tmin) - logt_step/2.0;
  double t_itt       = TMath::Exp(logt);
  double t_width     = 0.0;

  for (double t_step = 0; t_step<tstep; t_step++) {

    logt = logt + logt_step;
    t_itt = TMath::Exp(logt);
    t_width = t_itt*logt_step;

    if (tpi <= 1.0 && fRSPionXSec == false) {  
      xsec_stat =  utils::hadxs::berger::PionNucleusXSec(tpi, ppistar, t_itt, A, tpilow, siglow, tpihigh, sighigh);
      if(xsec_stat){
        LOG("BergerSehgalCohPi", pERROR) << "Call to PionNucleusXSec code failed - return xsec of 0.0";
        return 0.0;
      }
      dsigdzfit = siglow + (sighigh - siglow) * (tpi - tpilow) / (tpihigh - tpilow);
      dsigdtfit = dsigdzfit / (2.0 * ppistar * ppistar);
      // we are handed a cross section in mb, need to convert it to GeV^{-2}
      dsig +=   1.0  * front * Ga2 * t_width * dsigdtfit * units::mb;
    }
    else {
      dsig += /*factor **/ front * Ga2 * t_width * RS_factor * exp(-1.0*b*t_itt);
    }

  }
  xsec = dsig;

  // Correction for finite final state lepton mass.
  // Lepton mass modification is part of Berger-Sehgal and is not optional.
  if (pionIsCharged) {
    double C = 1.;
    // First, we need to remove the leading G_{A}^2 which is required for NC.
    xsec /= Ga2;
    // Next, \cos^2 \theta_{Cabibbo} appears in the CC xsec, but not the NC.
    xsec *= fCos8c2;
    double ml    = interaction->FSPrimLepton()->Mass();
    double ml2   = TMath::Power(ml,2);
    double Q2min = ml2 * y/(1-y);
    if(Q2 > Q2min) {
      double C1 = TMath::Power(Ga - 0.5 * Q2min / (Q2 + kPionMass2), 2);
      double C2 = 0.25 * y * Q2min * (Q2 - Q2min) / 
        TMath::Power(Q2 + kPionMass2, 2);
      C = C1 + C2;
    } else {
      C = 0.;
    }
    xsec *= (2. * C); // *2 is for CC vs NC in BS 
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("BergerSehgalCohPi", pDEBUG)
    << "\n momentum transfer .............. Q2    = " << Q2
    << "\n mass number .................... A     = " << A
    << "\n pion energy .................... Epi   = " << Epi
    << "\n propagator term ................ Ga2   = " << Ga2
    << "\n total pi+N cross section ....... sigT  = " << sigtot_pin
    << "\n inelastic pi+N cross section ... sigI  = " << siginel_pin
    << "\n nuclear size scale ............. Ro    = " << fRo
    << "\n pion absorption factor ......... Fabs  = " << fabs
    << "\n t integration range ............ [" << tmin << "," << tmax << "]"
  LOG("BergerSehgalCohPi", pINFO)
    << "d2xsec/dQ2dy[COHPi] (Q2= " << Q2 << ", y="
    << y << ", E=" << E << ") = "<< xsec;
#endif

  //----- The algorithm computes d^2xsec/dQ2dy
  // Check whether variable tranformation is needed? May be working with logs. 
  // kPSlogQ2logyfE is possible - all others will not succeed
  if(kps != kPSQ2yfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2yfE, kps);
    xsec *= J;
  }
  return xsec;
}
//____________________________________________________________________________
double BergerSehgalCOHPiPXSec2015::ExactKinematicTerm(const Interaction * interaction) const
{
  // This function is a bit inefficient but is being encapsulated as 
  // such in order to possibly migrate into a general kinematics check.
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  bool   pionIsCharged = interaction->ProcInfo().IsWeakCC();
  double M_pi          = pionIsCharged ? kPionMass : kPi0Mass;
  double E             = init_state.ProbeE(kRfLab);        // nu E
  double Q2            = kinematics.Q2();
  double y             = kinematics.y();                   // inelasticity
  double fp2           = (0.93 * M_pi)*(0.93 * M_pi); 

  double term = ((kGF2 * fp2) / (4.0 * kPi2)) * 
    ((E * (1.0 - y)) / sqrt(y*E * y*E + Q2)) * 
    (1.0 - Q2 / (4.0 * E*E * (1.0 - y)));
  return term;   
}
//____________________________________________________________________________
double BergerSehgalCOHPiPXSec2015::PionCOMAbsMomentum(const Interaction * interaction) const
{
  // This function is a bit inefficient but is being encapsulated as 
  // such in order to possibly migrate into a general kinematics check.
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  bool   pionIsCharged = interaction->ProcInfo().IsWeakCC();
  double M_pi          = pionIsCharged ? kPionMass : kPi0Mass;
  double E             = init_state.ProbeE(kRfLab);        // nu E
  double Q2            = kinematics.Q2();
  double y             = kinematics.y();                   // inelasticity
  double MT            = init_state.Tgt().Mass(); 

  double W2      = MT*MT - Q2 + 2.0 * y * E * MT;
  double arg     = (2.0*MT*(y*E - M_pi) - Q2 - M_pi*M_pi)*(2.0*MT*(y*E + M_pi) - Q2 - M_pi*M_pi);
  if (arg < 0) return arg;
  double ppistar = TMath::Sqrt(arg) / 2.0 / TMath::Sqrt(W2);

  return ppistar;
}
//____________________________________________________________________________
double BergerSehgalCOHPiPXSec2015::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool BergerSehgalCOHPiPXSec2015::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  int nu = init_state.ProbePdg();

  if (!proc_info.IsCoherent())  return false;
  if (!proc_info.IsWeak())      return false;
  if (target.HitNucIsSet())     return false;
  if (!(target.A()>1))          return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
void BergerSehgalCOHPiPXSec2015::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BergerSehgalCOHPiPXSec2015::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void BergerSehgalCOHPiPXSec2015::LoadConfig(void)
{
  GetParam( "COH-Ma",fMa ) ;
  GetParam( "COH-Ro", fRo ) ;

  double thc ;
  GetParam( "CabibboAngle", thc ) ;
  fCos8c2     = TMath::Power(TMath::Cos(thc), 2);

  // fRSPionXSec => Do not use the pion-nucleus cross section from Table 1 in PRD 79, 053003
  // Instead, use the Rein-Sehgal "style" pion-nucleon cross section and scale by A 
  // for all pion energies.
  GetParam( "COH-UseRSPionXSec", fRSPionXSec ) ;

  //-- load the differential cross section integrator
  fXSecIntegrator =
    dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

