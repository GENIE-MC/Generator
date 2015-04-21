//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Coherent/COHElasticPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//____________________________________________________________________________
COHElasticPXSec::COHElasticPXSec() :
XSecAlgorithmI("genie::COHElasticPXSec")
{

}
//____________________________________________________________________________
COHElasticPXSec::COHElasticPXSec(string config) :
XSecAlgorithmI("genie::COHElasticPXSec", config)
{

}
//____________________________________________________________________________
COHElasticPXSec::~COHElasticPXSec()
{

}
//____________________________________________________________________________
double COHElasticPXSec::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  if(!target.IsNucleus()) return 0;

  double E  = init_state.ProbeE(kRfLab);
  double Q2 = kinematics.Q2();
  int    Z  = target.Z(); 
  int    N  = target.N();

  // ...
  // ...
  // ...

  double xsec  = (0.25*kGF2/kPi2) * 
	         TMath::Power(N - (1-4*fSin2thw)*Z, 2) /* * ...  */;

  LOG("COHEl", pDEBUG)
    << "dXSec[vA,COHEl]/dQ2 (Ev = "<< E<< ", Q2 = "<< Q2 << ") = "<< xsec;

  //-- The algorithm computes dxsec/dQ2
  //   Check whether variable tranformation is needed
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
    xsec *= J;
  }

  return xsec;
}
//____________________________________________________________________________
double COHElasticPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool COHElasticPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info = interaction->ProcInfo();  
  if(!proc_info.IsCoherentElas()) return false;

  const InitialState & init_state = interaction->InitState();
  const Target & target = init_state.Tgt();
  if(!target.IsNucleus()) return false;

  return true;
}
//____________________________________________________________________________
void COHElasticPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHElasticPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHElasticPXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  double thw = fConfig->GetDoubleDef(
                       "WeinbergAngle", gc->GetDouble("WeinbergAngle"));
  fSin2thw = TMath::Power(TMath::Sin(thw), 2);

  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
