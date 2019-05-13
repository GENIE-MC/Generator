//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - February 15, 2008

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 15, 2009 - CA
   This class was first added in version 2.5.1. 
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/Diffractive/XSection/ReinDFRPXSec.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//____________________________________________________________________________
ReinDFRPXSec::ReinDFRPXSec() :
XSecAlgorithmI("genie::ReinDFRPXSec")
{

}
//____________________________________________________________________________
ReinDFRPXSec::ReinDFRPXSec(const std::string & config) :
XSecAlgorithmI("genie::ReinDFRPXSec", config)
{

}
//____________________________________________________________________________
ReinDFRPXSec::~ReinDFRPXSec()
{

}
//____________________________________________________________________________
double ReinDFRPXSec::XSec(
        const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target &       target     = init_state.Tgt();

  bool   isCC    = interaction->ProcInfo().IsWeakCC();
  double E       = init_state.ProbeE(kRfHitNucRest);          // neutrino energy
  double x       = kinematics.x();                            // bjorken x
  double y       = kinematics.y();                            // inelasticity y
  double t       = kinematics.t();                            // (magnitude of) square of four-momentum xferred to proton
  double M       = target.HitNucMass();                       //
  double Q2      = 2.*x*y*M*E;                                // momentum transfer Q2>0
  double Gf      = kGF2 * M/(16*kPi3);                        // GF/pi/etc factor
  double fp      = 0.93 * kPionMass;                          // pion decay constant (cc)
  double fp2     = TMath::Power(fp,2.);         
  double Epi     = y*E - t/(2*M);                             // pion energy.  note we use - instead of + like Rein's paper b/c our t is magnitude only
  double b       = fBeta;
  double ma2     = TMath::Power(fMa,2);
  double propg   = TMath::Power(ma2/(ma2+Q2),2.);             // propagator term
  double sTot    = utils::hadxs::TotalPionNucleonXSec(Epi, isCC);   // pi+N total cross section; CC process always produces a charged pion
  double sTot2   = TMath::Power(sTot,2.);
  double tFac    = TMath::Exp(-b*t);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("ReinDFR", pDEBUG)
    << "E = " << E << ", x = " << x << ", y = " << y << ", Q2 = " << Q2;
  LOG("ReinDFR", pDEBUG)
    << "Epi = " << Epi << ", s^{piN}_{tot} = " << sTot;
  LOG("ReinDFR", pDEBUG)
    << "b = " << b << ", t = [" << tmin << ", " << tmax << "]";
#endif

  // fixme: WARNING: don't leave this in here!!!
  // only needed for comparison to Rein's paper
//  double W2 = M*M + 2 * M * y * E - Q2;
//  if (W2 < 4)
//    return 0;

  //----- compute d^2sigma/dxdydt
  double xsec = Gf*E*fp2*(1-y)*propg*sTot2*tFac;

  // NC XS is half of CC
  if (!isCC)
    xsec *= 0.5;
 
  //----- Check whether variable tranformation is needed
  if(kps!=kPSxytfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSxytfE,kps);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("ReinDFR", pDEBUG)
     << "Jacobian for transformation to: " 
                  << KinePhaseSpace::AsString(kps) << ", J = " << J;
#endif
    xsec *= J;
  }

  //----- if requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- number of scattering centers in the target
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 
  xsec *= NNucl; 

  return xsec;
}
//____________________________________________________________________________
double ReinDFRPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool ReinDFRPXSec::ValidProcess(const Interaction * interaction) const
{
  //expect only free protons as target
  const InitialState & init_state = interaction -> InitState();
  const Target &       target     = init_state.Tgt();
  if(target.A() > 1 || target.Z() != 1)
    return false;

  if(interaction->TestBit(kISkipProcessChk)) return true;

  if(interaction->ProcInfo().IsDiffractive()) return true;
  return false;
}
//____________________________________________________________________________
void ReinDFRPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinDFRPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ReinDFRPXSec::LoadConfig(void)
{
	GetParam( "DFR-Ma", fMa ) ;
	GetParam( "DFR-Beta",fBeta ) ;

	fXSecIntegrator =
			dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator") );
	assert(fXSecIntegrator);
}
//____________________________________________________________________________

