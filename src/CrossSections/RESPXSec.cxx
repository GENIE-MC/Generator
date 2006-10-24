//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - March 03, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Conventions/Constants.h"
#include "Conventions/KineVar.h"
#include "CrossSections/RESPXSec.h"
#include "CrossSections/GXSecFunc.h"
#include "Messenger/Messenger.h"
#include "Numerical/IntegratorI.h"
#include "PDG/PDGUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
RESPXSec::RESPXSec() :
XSecAlgorithmI("genie::RESPXSec")
{

}
//____________________________________________________________________________
RESPXSec::RESPXSec(string config) :
XSecAlgorithmI("genie::RESPXSec", config)
{

}
//____________________________________________________________________________
RESPXSec::~RESPXSec()
{

}
//____________________________________________________________________________
double RESPXSec::XSec(
                 const Interaction * interaction, KinePhaseSpace_t kps) const
{
  assert(kps==kPSWfE || kps==kPSQ2fE || kps==kPSq2fE);

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  GXSecFunc * func = 0;
  double xsec = 0;

  if (kps==kPSQ2fE || kps==kPSq2fE) {

    double Q2 = interaction->Kine().Q2();
    func = new Integrand_D2XSec_DWDQ2_EQ2(
                         fPartialXSecAlg, interaction, Q2);

    Range1D_t rW = utils::kinematics::KineRange(interaction, kKVW);
    assert(rW.min < rW.max && rW.min>0);

    func->SetParam(0,"W",rW);
    xsec = fIntegrator->Integrate(*func);

  } else if (kps==kPSWfE) {

    double W = interaction->Kine().W();
    func = new Integrand_D2XSec_DWDQ2_EW(
                         fPartialXSecAlg, interaction, W);

    Range1D_t rQ2 = utils::kinematics::KineRange(interaction, kKVQ2);
    assert(rQ2.min > 0 && rQ2.max > rQ2.min);

    func->SetParam(0,"Q2",rQ2);
    xsec = fIntegrator->Integrate(*func);

  } else {
   LOG("DISPXSec", pFATAL)
         << "Can not handle phase space = "  << KinePhaseSpace::AsString(kps);
    exit(1);
  }
  delete func;
  return xsec;
}
//____________________________________________________________________________
bool RESPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsResonant()) return false;
  if(!proc_info.IsWeak())     return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  return true;
}
//____________________________________________________________________________
bool RESPXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> InitState();
  double Ev  = init_state.ProbeE(kRfHitNucRest);

  double EvThr = interaction->EnergyThreshold();
  if(Ev <= EvThr) return false;

  return true;
}
//____________________________________________________________________________
void RESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESPXSec::LoadConfig(void)
{
  fPartialXSecAlg = 0;
  fIntegrator     = 0;

  //-- get the requested d^2xsec/dxdy xsec algorithm to use
  fPartialXSecAlg =
         dynamic_cast<const XSecAlgorithmI *> (this->SubAlg("DiffXSecAlg"));

  //-- get the specified integration algorithm
  fIntegrator = 
         dynamic_cast<const IntegratorI *> (this->SubAlg("Integrator"));

  assert(fPartialXSecAlg);
  assert(fIntegrator);
}
//____________________________________________________________________________

