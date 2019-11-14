//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - March 11, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Nov 21, 2007 - CA
   Was renamed to ReinSehgalCOHPiPXSec (from ReinSehgalCOHPXSec)
 @ Mar 31, 2009 - CA
   Fixed a minor bug in the C2 term controlling the forward mu- suppression
   predicted by Adler's PCAC (mpi -> mpi^2)
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
#include "Physics/Coherent/XSection/ReinSehgalCOHPiPXSec.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
ReinSehgalCOHPiPXSec::ReinSehgalCOHPiPXSec() :
XSecAlgorithmI("genie::ReinSehgalCOHPiPXSec")
{

}
//____________________________________________________________________________
ReinSehgalCOHPiPXSec::ReinSehgalCOHPiPXSec(string config) :
XSecAlgorithmI("genie::ReinSehgalCOHPiPXSec", config)
{

}
//____________________________________________________________________________
ReinSehgalCOHPiPXSec::~ReinSehgalCOHPiPXSec()
{

}
//____________________________________________________________________________
double ReinSehgalCOHPiPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();

  //----- Compute the coherent NC pi0 production d2xsec/dxdy
  //      see page 34 in Nucl.Phys.B223:29-144 (1983)
  double E      = init_state.ProbeE(kRfLab); // neutrino energy
  double x      = kinematics.x(); // bjorken x
  double y      = kinematics.y(); // inelasticity y
  double Q2     = 2.*x*y*kNucleonMass*E;  // momentum transfer Q2>0
  double A      = (double) init_state.Tgt().A(); // mass number
  double A2     = TMath::Power(A,2.);
  double A_3    = TMath::Power(A,1./3.);
  double Gf     = kGF2 * kNucleonMass / (32 * kPi3);
  double fp     = 0.93 * kPionMass; // pion decay constant
  double fp2    = TMath::Power(fp,2.);
  double Epi    = y*E; // pion energy
  double ma2    = TMath::Power(fMa,2);
  double propg  = TMath::Power(ma2/(ma2+Q2),2.); // propagator term
  double r2     = TMath::Power(fReIm,2.);
  double sTot   = utils::hadxs::TotalPionNucleonXSec(Epi); // tot. pi+N xsec
  double sTot2  = TMath::Power(sTot,2.);
  double sInel  = utils::hadxs::InelasticPionNucleonXSec(Epi); // inel. pi+N xsec
  double Ro2    = TMath::Power(fRo*units::fermi,2.);

  // effect of pion absorption in the nucleus
  double Fabs   = TMath::Exp( -9.*A_3*sInel / (16.*kPi*Ro2) );

  // the xsec in Nucl.Phys.B223:29-144 (1983) is d^3xsec/dxdydt but the only
  // t-dependent factor is an exp(-bt) so it can be integrated analyticaly
  double Epi2   = TMath::Power(Epi,2.);
  double R      = fRo * A_3 * units::fermi; // nuclear radius
  double R2     = TMath::Power(R,2.);
  double b      = 0.33333 * R2;
  double MxEpi  = kNucleonMass*x/Epi;
  double mEpi2  = kPionMass2/Epi2;
  double tA     = 1. + MxEpi - 0.5*mEpi2;
  double tB     = TMath::Sqrt(1. + 2*MxEpi) * TMath::Sqrt(1.-mEpi2);
  double tmin   = 2*Epi2 * (tA-tB);
  double tmax   = 2*Epi2 * (tA+tB);
  double tint   = (TMath::Exp(-b*tmin) - TMath::Exp(-b*tmax))/b; // t integral

  double xsec = Gf*fp2 * A2 * E*(1-y) * sTot2 * (1+r2)*propg * Fabs*tint;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalCohPi", pDEBUG)
      << "\n momentum transfer .............. Q2    = " << Q2
      << "\n mass number .................... A     = " << A
      << "\n pion energy .................... Epi   = " << Epi
      << "\n propagator term ................ propg = " << propg
      << "\n Re/Im of fwd pion scat. ampl. .. Re/Im = " << fReIm
      << "\n total pi+N cross section ....... sigT  = " << sTot
      << "\n inelastic pi+N cross section ... sigI  = " << sInel
      << "\n nuclear size scale ............. Ro    = " << fRo
      << "\n pion absorption factor ......... Fabs  = " << Fabs
      << "\n t integration range ............ [" << tmin << "," << tmax << "]"
      << "\n t integration factor ........... tint  = " << tint;
#endif

  // compute the cross section for the CC case

  if(interaction->ProcInfo().IsWeakCC()) { 
     // Check whether a modification to Adler's PCAC theorem is applied for
     // including the effect of the muon mass. 
     // See Rein and Sehgal, PCAC and the Deficit of Forward Muons in pi+ 
     // Production by Neutrinos, hep-ph/0606185
     double C = 1.;
     if(fModPCAC) {
        double ml    = interaction->FSPrimLepton()->Mass();
        double ml2   = TMath::Power(ml,2);
        double Q2min = ml2 * y/(1-y);
        if(Q2>Q2min) {
           double C1    = TMath::Power(1-0.5*Q2min/(Q2+kPionMass2), 2);
           double C2    = 0.25*y*Q2min*(Q2-Q2min)/ TMath::Power(Q2+kPionMass2,2);
           C = C1+C2;
        } else {
           C = 0.;
        }
     }
     xsec *= (2.*C); 
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinSehgalCohPi", pINFO)
         << "d2xsec/dxdy[COHPi] (x= " << x << ", y="
                       << y << ", E=" << E << ") = "<< xsec;
#endif

  //----- The algorithm computes d^2xsec/dxdy
  //      Check whether variable tranformation is needed
  if(kps!=kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSxyfE,kps);
    xsec *= J;
  }

  return xsec;
}
//____________________________________________________________________________
double ReinSehgalCOHPiPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool ReinSehgalCOHPiPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();

  int nu = init_state.ProbePdg();

  if (!proc_info.IsCoherentProduction())  return false;
  if (!proc_info.IsWeak())      return false;
  if (target.HitNucIsSet())     return false;
  if (!(target.A()>1))          return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
void ReinSehgalCOHPiPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalCOHPiPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinSehgalCOHPiPXSec::LoadConfig(void)
{

  GetParam( "COH-Ma", fMa ) ;
  GetParam( "COH-Ro", fRo ) ;
  GetParam( "COH-ReImAmpl", fReIm ) ;
  GetParam( "COH-UseModifiedPCAC", fModPCAC ) ;

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

