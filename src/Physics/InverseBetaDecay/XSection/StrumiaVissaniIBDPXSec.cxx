//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Corey Reed <cjreed \at nikhef.nl>
         Nikhef - June 22, 2009

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/InverseBetaDecay/XSection/Constants.h"
#include "Physics/InverseBetaDecay/XSection/StrumiaVissaniIBDPXSec.h"

using namespace genie;
using namespace genie::constants;

ClassImp(StrumiaVissaniIBDPXSec)

//____________________________________________________________________________
StrumiaVissaniIBDPXSec::StrumiaVissaniIBDPXSec() :
   XSecAlgorithmI("genie::StrumiaVissaniIBDPXSec"),
   fXSecIntegrator(0)
{

}
//____________________________________________________________________________
StrumiaVissaniIBDPXSec::StrumiaVissaniIBDPXSec(string config) :
   XSecAlgorithmI("genie::StrumiaVissaniIBDPXSec", config),
   fXSecIntegrator(0)
{

}
//____________________________________________________________________________
StrumiaVissaniIBDPXSec::~StrumiaVissaniIBDPXSec()
{

}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::XSec(
       const Interaction * interaction, KinePhaseSpace_t kps) const
{
  // compute the differential cross section ds/dt

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const double         Ev     = init_state.ProbeE(kRfHitNucRest);
  const Target &       target = init_state.Tgt();
  const bool           isProt = target.IsProton();
  const Kinematics &   kine   = interaction->Kine();
  const double         q2     = kine.q2();
  const double         mp     = (isProt) ? kProtonMass   : kNeutronMass;
  const double         mp2    = (isProt) ? kProtonMass2  : kNeutronMass2;
  const double         mn2    = (isProt) ? kNeutronMass2 : kProtonMass2;
  
  // calculate s-u and s-m_nucleon^2 in nucleon rest frame
  // need E_e
  const double Ee  = Ev + ( (q2 - mn2 + mp2) / (2.000*mp) );
  const double sMinusU   = ((2.000*mp*(Ev+Ee)) - kElectronMass2)
                           * (isProt ? 1.000 : -1.000);
  const double sMinusMp2 = 2.000*mp*Ev;
  
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("StrumiaVissani", pDEBUG) << "*** Ev = " << Ev 
				<< ", q2 = " << q2
				<< ", Ee = " << Ee
				<< ", s-u = " << sMinusU
				<< ", s-Mp2 = " << sMinusMp2;
#endif

        double xsec = dSigDt(sMinusU, sMinusMp2, q2);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("StrumiaVissani", pDEBUG) << "*** dSdt = " << xsec;
#endif

  // apply correction factors
  const double rdcf = RadiativeCorr(Ee);
  const double fscf = (isProt) ? 1.00 : FinalStateCorr(Ee);
  xsec *= rdcf * fscf;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("StrumiaVissani", pDEBUG) << "*** rad.corr. = " << rdcf
				<< ", fin.st.cor. = " << fscf
				<< ", xsec = " << xsec;
#endif
   
  //----- The algorithm computes dxsec/dt, t=q^2
  //      Check whether variable tranformation is needed
  if(kps!=kPSq2fE && kps!=kPSQ2fE) {
     const double J = utils::kinematics::Jacobian(interaction,kPSq2fE,kps);
     xsec *= J;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("StrumiaVissani", pDEBUG) << "*** Jacobian = " << J;
#endif
  }

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("StrumiaVissani", pDEBUG) << "*** xsec = " << xsec;
#endif
  return xsec;
}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::Integral(const Interaction * interaction) const
{
   // Compute the total cross section for a free nucleon target

   assert(interaction!=0);
   if(! this -> ValidProcess    (interaction) ) return 0.;
   if(! this -> ValidKinematics (interaction) ) return 0.;
   
   assert(fXSecIntegrator!=0);
   double xsec = fXSecIntegrator->Integrate(this, interaction);

   return xsec;
}
//____________________________________________________________________________
bool StrumiaVissaniIBDPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;
  
  // should be IBD and either nu_e + n or anu_e + p
  if (interaction->ProcInfo().IsInverseBetaDecay()) {
  
     const InitialState & init_state = interaction -> InitState();
     const Target & target = init_state.Tgt();
     if ( (target.IsProton() && pdg::IsAntiNuE(init_state.ProbePdg())) || (target.IsNeutron() && pdg::IsNuE(init_state.ProbePdg()) )) 
		return true;
  }	
  LOG("StrumiaVissani", pERROR) << "*** Should be IBD processes either nu_e + n or anu_e + p!";
  return false;
}
//____________________________________________________________________________
bool StrumiaVissaniIBDPXSec::ValidKinematics(const Interaction* interaction) const
{
   if(interaction->TestBit(kISkipKinematicChk)) return true;

   const InitialState & init_state = interaction -> InitState();
   const double Ev = init_state.ProbeE(kRfHitNucRest);
   static const double Ecutoff = kNucleonMass / 2;
   if (Ev > Ecutoff) {
      LOG("StrumiaVissani", pERROR) << "*** Ev=" << Ev
				    << " is too large for VLE QEL!";
   } else if (init_state.IsNuBarP() || init_state.IsNuN()) {
      const KPhaseSpace & kps = interaction->PhaseSpace();
      return kps.IsAboveThreshold();
   }

   return false;
}
//____________________________________________________________________________
void StrumiaVissaniIBDPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void StrumiaVissaniIBDPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void StrumiaVissaniIBDPXSec::LoadConfig(void)
{

   // cabibbo angle
   double cab ;
   GetParam( "CabibboAngle", cab ) ;
   const double cosCab = TMath::Cos(cab);
   fCosCabibbo2 = cosCab*cosCab;
   
   // form factor params
   GetParam( "QEL-FA0", fg1of0 ) ;
   double ma, mv ;
   GetParam( "QEL-Ma", ma ) ;
   GetParam( "QEL-Mv", mv ) ;
   fMa2 = ma*ma;
   fMv2 = mv*mv;

   // magnetic moments

   double mup, mun ;
   GetParam( "AnomMagnMoment-P", mup );
   GetParam( "AnomMagnMoment-N", mun );
   fNucleonMMDiff    = (mup - 1.000) - mun; // *anamolous* mag. mom. diff.

   // numeric
   int epmag ;
   GetParam("EpsilonMag", epmag ) ;
   fEpsilon = TMath::Power(10.000, -1.000 * static_cast<double>(epmag) );
   
   LOG("StrumiaVissani", pINFO) << "*** USING: cos2(Cabibbo)=" 
				<< fCosCabibbo2
				<< ", g1(0)=" << fg1of0
				<< ", Ma^2=" << fMa2
				<< ", Mv^2=" << fMv2
				<< ", xi=" << fNucleonMMDiff
				<< ", epsilon=" << fEpsilon;
   
   // load XSec Integrator
   fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
   assert(fXSecIntegrator!=0);
   
}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::dSigDt(const double sMinusU,
				      const double sMinusMnuc,
				      const double t) const
{
   // return the differential cross section dS/dt. eqn 3 from reference
   // t = q^2
   //
   // for anue + p -> e+ + n , sMinusU = s - u  and  sMinusMnuc = s - m_p^2
   // for nue  + n -> e- + p , sMinusU = u - s  and  sMinusMnuc = s - m_n^2
   // where s = (p_nu + p_p)^2 and u = (p_nu - p_n)^2
   
   const double numer = kGF2 * fCosCabibbo2;
   const double denom = (2.000 * kPi) * sMinusMnuc*sMinusMnuc;
   assert(TMath::Abs(denom) > fEpsilon);  // avoid divide by zero
   
   return ( (numer / denom) * MtxElm(sMinusU, t) );
}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::MtxElm(const double sMinusU,
				      const double t) const
{
   // return the square of the matrix element. eqn 5 from the reference paper
   // |M^2| = A(t) - (s-u)B(t) + (s-u)^2 C(t)
   //
   // factor 16 comes from eqn 6
   //
   // for anue + p -> e+ + n , sMinusU = s - u
   // for nue  + n -> e- + p , sMinusU = u - s
   // where s = (p_nu + p_p)^2 and u = (p_nu - p_n)^2
   
   // store multiply used variables to reduce number of calculations
   const double t4m      = t / k4NucMass2;
   const double t2       = t * t;
   
   const double fdenomB  = 1.000 - (t / fMv2);
   const double fdenom   = (1.000 - t4m)*fdenomB*fdenomB;
   const double g1denom  = 1.000 - (t / fMa2);
   const double g2denom  = kPionMass2 - t;
   assert(TMath::Abs(fdenom)  > fEpsilon);   // avoid divide by zero
   assert(TMath::Abs(g1denom) > fEpsilon);   // avoid divide by zero
   assert(TMath::Abs(g2denom) > fEpsilon);   // avoid divide by zero
   
   const double f1       = (1.000 - ( (1.000+fNucleonMMDiff) * t4m)) / fdenom;
   const double f2       = fNucleonMMDiff / fdenom;
   const double g1       = fg1of0 / (g1denom*g1denom);
   const double g2       = (2.000 * kNucleonMass2 * g1) / g2denom;
   
   const double f12      = f1 * f1;
   const double f124     = 4.000 * f12;
   const double f22      = f2 * f2;
   const double g12      = g1 * g1;
   const double g124     = 4.000 * g12;
   const double g22      = g2 * g2;
   const double g224meM2 = 4.000 * g22 * kElectronMass2 / kNucleonMass2;
   
   const double g1cFsumR   = g1 * (f1+f2);
   const double f1cf2      = f1 * f2;
   const double g1cg2      = g1 * g2;
   const double f1cf2R8    = f1cf2 * 8.000;
   const double g1cg2R16me = g1cg2 * 16.000 * kElectronMass2;
   
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
   LOG("StrumiaVissani", pDEBUG) << "*** t = " << t
				 << ", g2 = " << g2
				 << ", g1 = " << g1
				 << ", f1 = " << f1
				 << ", f2 = " << f2;
#endif
   
   // ok - now calculate the terms of the matrix element
   const double mat = MAterm(t, t2, f124, f22, g124, 
                             g224meM2, f1cf2R8, g1cg2R16me, g1cFsumR);
   const double mbt = MBterm(t, f1cf2, g1cg2, g1cFsumR, f22);
   const double mct = MCterm(t, f124, f22, g124);
   
   const double M2  = mat - (sMinusU * (mbt - (sMinusU * mct)));
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
   LOG("StrumiaVissani", pDEBUG) << "*** Matrix element = " << (M2 / 16.000)
				 << " [16A=" << mat << "] [16B=" << mbt
				 << "] [16C=" << mct << "]";
#endif
   
   return ( M2 / 16.000 );
}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::MAterm(const double t,
				      const double t2,
				      const double f124,
				      const double f22,
				      const double g124,
				      const double g224meM2,
				      const double f1cf2R8,
				      const double g1cg2R16me,
				      const double g1cFsumR)
{
   // return A(t) in the matrix element calculation
   // from eqn 6 (without the factor 16)
   //
   // for speed purposes, no check that |t|^2 = t2 is performed
   
   // store multiply used terms
   const double tmme2 = t - kElectronMass2;
   const double tpme2 = t + kElectronMass2;
   
   double r1  = f124 * ( k4NucMass2 + tpme2);
          r1 += g124 * (-k4NucMass2 + tpme2);
          r1 += f22  * ( (t2/kNucleonMass2) + 4.000*tpme2 );
          r1 += g224meM2 * t;
          r1 += f1cf2R8 * (2.000*t + kElectronMass2);
          r1 += g1cg2R16me;
          r1 *= tmme2;
   
   double r2  = (f124 + (t * (f22 / kNucleonMass2))) * (k4NucMass2 + tmme2);
          r2 +=  g124                                * (k4NucMass2 - tmme2);
          r2 += g224meM2 * tmme2;
          r2 += f1cf2R8 * (2.000*t - kElectronMass2);
          r2 += g1cg2R16me;
          r2 *= kNucMassDiff2;
   
   const double r3 = 32.000 * kElectronMass2 * kNucleonMass * kNucMassDiff
                        * g1cFsumR;
   
   return ( r1 - r2 - r3 );
}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::MBterm(const double t,
				      const double f1cf2,
				      const double g1cg2,
				      const double g1cFsumR,
				      const double f22)
{
   // return C(t) in the matrix element calculation
   // from eqn 6 (without the factor 16)
   
   double bterm  = 16.000 * t * g1cFsumR;
          bterm += ( k4EleMass2 * kNucMassDiff
                     * (f22 + (f1cf2 + 2.000*g1cg2)) ) / kNucleonMass2;
   return bterm;
}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::MCterm(const double t,
				      const double f124,
				      const double f22,
				      const double g124)
{
   // return C(t) in the matrix element calculation
   // from eqn 6 (without the factor 16)
   
   return ( f124 + g124 - (t * ( f22 / kNucleonMass2 )) );
}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::RadiativeCorr(const double Ee) const
{
   // radiative correction to the cross section. eqn 14 from the reference
   // only valid for Ev<<m_p!
   
   assert(Ee > fEpsilon); // must be non-zero and positive
   double rc  = 6.000 + (1.500 * TMath::Log(kProtonMass / (2.000*Ee)));
          rc += 1.200 * TMath::Power((kElectronMass / Ee), 1.500);
          rc *= kAem / kPi;
	  rc += 1.000;
   return rc;

}
//____________________________________________________________________________
double StrumiaVissaniIBDPXSec::FinalStateCorr(const double Ee) const
{
   // Sommerfeld factor; correction for final state interactions.
   // eqn 15 of the reference
   
   assert(Ee > fEpsilon); // must be non-zero and positive
   const double eta  = 2.000*kPi*kAem
                       / TMath::Sqrt(1.000 - (kElectronMass2 / (Ee*Ee)));
   const double expn = TMath::Exp(-1.000 * eta);
   assert(expn < 1.000);
   return ( eta / (1.000 - expn) );
}
//____________________________________________________________________________

