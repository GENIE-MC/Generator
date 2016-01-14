//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ June 03, 2009 - CA
   Was first added in v2.5.1

*/
//____________________________________________________________________________

#include <TMath.h>
#include <TSystem.h>

#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"
#include "Base/XSecIntegratorI.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "GiBUU/GiBUURESPXSec.h"
#include "GiBUU/GiBUUData.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
GiBUURESPXSec::GiBUURESPXSec() :
XSecAlgorithmI("genie::GiBUURESPXSec")
{

}
//____________________________________________________________________________
GiBUURESPXSec::GiBUURESPXSec(string config) :
XSecAlgorithmI("genie::GiBUURESPXSec", config)
{

}
//____________________________________________________________________________
GiBUURESPXSec::~GiBUURESPXSec()
{

}
//____________________________________________________________________________
double GiBUURESPXSec::XSec(
  const Interaction * interaction, KinePhaseSpace_t /*kps*/) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  double xsec = 0;

/*
  // Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  double W  = kinematics.W();
  double Q2 = kinematics.Q2();

  // Under the DIS/RES joining scheme, xsec(RES)=0 for W>=Wcut
  if(fUsingDisResJoin) {
    if(W>=fWcut) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("GiBUURes", pDEBUG)
         << "RES/DIS Join Scheme: XSec[RES, W=" << W 
         << " >= Wcut=" << fWcut << "] = 0";
#endif
       return 0;
    }
  }

  // Get info about the initial state, and procces type
  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target &       target     = init_state.Tgt();

  double E      = init_state.ProbeE(kRfHitNucRest);
  double Mnuc   = target.HitNucMass(); 
  int  nucpdgc  = target.HitNucPdg();
  int  nupdgc   = init_state.ProbePdg();
  bool is_nu    = pdg::IsNeutrino     (nupdgc);
  bool is_nubar = pdg::IsAntiNeutrino (nupdgc);
  bool is_p     = pdg::IsProton       (nucpdgc);
  InteractionType_t itype = proc_info.InteractionTypeId();

  // Get the input baryon resonance
  Resonance_t resonance = interaction->ExclTag().Resonance();
  string      resname   = utils::res::AsString(resonance);
  bool        is_delta  = utils::res::IsDelta (resonance);

//  bool is_CC    = proc_info.IsWeakCC();
//  if(is_CC && !is_delta) {
//    if((is_nu && is_p) || (is_nubar && is_n)) return 0;
//  }


  // Get baryon resonance parameters
//  fBRP.RetrieveData(resonance);  
//  double Mres = fBRP.Mass();
//  double Gres = fBRP.Width();
//  int    Nres = fBRP.ResonanceIndex();

  // Get the GiBUU form factor data
  GiBUUData * gibuu_data = GiBUUData::Instance();

  // Calculate the double differential cross section d2sigma/dWdQ2

  const GiBUUData::FormFactors & ff = gibuu_data->FF();
  if(is_delta) {
    //
    // Delta resonances
    //
    double F1V = ff.F1V(Q2, resonance, nucpdgc, itype);
    double F2V = ff.F2V(Q2, resonance, nucpdgc, itype);
    double FA  = ff.FA (Q2, resonance, nucpdgc, itype);
    double FP  = ff.FP (Q2, resonance, nucpdgc, itype);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GiBUURes", pINFO) 
      << "\n F1V = " << F1V << ", F2V = " << F2V
      << ", FA = " << FA << ", FP = " << FP;
#endif

  } 
  else {
    //
    // N resonances
    //
    double C3V = ff.C3V(Q2, resonance, nucpdgc, itype);
    double C4V = ff.C4V(Q2, resonance, nucpdgc, itype);
    double C5V = ff.C5V(Q2, resonance, nucpdgc, itype);
    double C6V = ff.C6V(Q2, resonance, nucpdgc, itype);
    double C3A = ff.C3A(Q2, resonance, nucpdgc, itype);
    double C4A = ff.C4A(Q2, resonance, nucpdgc, itype);
    double C5A = ff.C5A(Q2, resonance, nucpdgc, itype);
    double C6A = ff.C6A(Q2, resonance, nucpdgc, itype);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GiBUURes", pINFO) 
      << "\n C3V = " << C3V << ", C4V = " << C4V 
      <<  ", C5V = " << C5V << ", C6V = " << C6V
      << "\n C3A = " << C3A << ", C4A = " << C4A 
      <<  ", C5A = " << C5A << ", C6A = " << C6A;
#endif

  } // Delta or N

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("GiBUURes", pINFO) 
    << "\n d2xsec/dQ2dW"  << "[" << interaction->AsString()
          << "](W=" << W << ", Q2=" << Q2 << ", E=" << E << ") = " << xsec;
#endif

  // The algorithm computes d^2xsec/dWdQ2
  // Check whether variable tranformation is needed
  if(kps!=kPSWQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSWQ2fE,kps);
    xsec *= J;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Number of scattering centers in the target
  int NNucl = (is_p) ? target.Z() : target.N();
  xsec*=NNucl; // nuclear xsec (no nuclear suppression factor)
*/

  return xsec;
}
//____________________________________________________________________________
double GiBUURESPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool GiBUURESPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const XclsTag &      xcls       = interaction->ExclTag();

  if(!proc_info.IsResonant()) return false;
  if(!proc_info.IsWeak())     return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  if (!pdg::IsProton(nuc)  && !pdg::IsNeutron(nuc))     return false;
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;

  if(!xcls.KnownResonance()) return false;

  return true;
}
//____________________________________________________________________________
void GiBUURESPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GiBUURESPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GiBUURESPXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  // Load all configuration data or set defaults

  double ma  = fConfig->GetDoubleDef( "Ma", gc->GetDouble("RES-Ma") );
  double mv  = fConfig->GetDoubleDef( "Mv", gc->GetDouble("RES-Mv") );

  fMa2 = TMath::Power(ma,2);
  fMv2 = TMath::Power(mv,2);

  //-- Use algorithm within a DIS/RES join scheme. If yes get Wcut
  fUsingDisResJoin = fConfig->GetBoolDef(
                           "UseDRJoinScheme", gc->GetBool("UseDRJoinScheme"));
  fWcut = 999999;
  if(fUsingDisResJoin) {
    fWcut = fConfig->GetDoubleDef("Wcut",gc->GetDouble("Wcut"));
  }

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________

