//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Hugh Gallagher
         Tufts University

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Physics/Strange/XSection/PaisQELLambdaPXSec.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
PaisQELLambdaPXSec::PaisQELLambdaPXSec() :
XSecAlgorithmI("genie::PaisQELLambdaPXSec")
{

}
//____________________________________________________________________________
PaisQELLambdaPXSec::PaisQELLambdaPXSec(string config) :
XSecAlgorithmI("genie::PaisQELLambdaPXSec", config)
{

}
//____________________________________________________________________________
PaisQELLambdaPXSec::~PaisQELLambdaPXSec()
{

}
//____________________________________________________________________________
double PaisQELLambdaPXSec::XSec(
                  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get kinematics & init state - compute auxiliary vars
  const Kinematics &   kinematics  = interaction->Kine();
  const InitialState & init_state  = interaction->InitState();
  const Target &       target      = init_state.Tgt();

  //neutrino energy & momentum transfer
  double E   = init_state.ProbeE(kRfHitNucRest);
  double E2  = E * E;
  double q2  = kinematics.q2();


  //resonance mass & nucleon mass
  double Mnuc  = target.HitNucMass();
  double Mnuc2 = TMath::Power(Mnuc,2);

  //----- Calculate the differential cross section dxsec/dQ^2
  double Gf        = kGF2 / (2*kPi);
  double ml        = interaction->FSPrimLepton()->Mass();
  double ml2       = TMath::Power(ml,2);
  double M1        = Mnuc;
  double M2        = (this)->MHyperon(interaction);
  double v         = (TMath::Power(M2,2) - Mnuc2 - q2) / (2*Mnuc); 
  double v2        = TMath::Power(v,2);
  double s         = Mnuc2 + 2*Mnuc*E;
  double u         = Mnuc2 + ml2 + 2*v*Mnuc - 2*Mnuc*E;

// xsec term changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = (is_neutrino) ? -1 : 1;

// Calculate the QEL form factors
  fFormFactors.Calculate(interaction);

  double F1V   = fFormFactors.F1V();
  double xiF2V = fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
//  double Fp    = fFormFactors.Fp();

// calculate w coefficients
   //start with Mass terms
  double Mp    = M2 + M1;
  double Mm    = M2 - M1;
  double Mm2   = TMath::Power(Mm, 2);
  double Mp2   = TMath::Power(Mp, 2);

   //Powers of Form Factors
  double FA2   = TMath::Power(FA, 2);
//  double FA3   = 0;  

   //Calculate W terms
 
  double w1 = (Mm2 - q2)/(4*Mnuc2)*TMath::Power((F1V + xiF2V), 2) + (Mp2 - q2)/(4*Mnuc2) * FA2;
  double w2 = FA2 + TMath::Power((F1V + xiF2V - Mp * xiF2V / (2 * Mnuc)), 2) - q2 / Mnuc2 * TMath::Power((xiF2V / 2), 2);
  double w3 = 2 * FA * (F1V + xiF2V); 

  double xsec = Gf*fSin8c2 / (16*Mnuc2*E2) * (-8*Mnuc2*q2*w1 - 4*(Mnuc2*v2 - q2)*w2 - sign*2*(s - u)*q2*w3 + (s-u)*(s-u)*w2);
  xsec = TMath::Max(xsec,0.);

  //----- The algorithm computes dxsec/dQ2
  //      Check whether variable tranformation is needed
  if(kps!=kPSQ2fE) {
    double J = utils::kinematics::Jacobian(interaction,kPSQ2fE,kps);
    xsec *= J;
  }

  //----- If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- Nuclear cross section (simple scaling here)
  int nuc   = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nuc)) ? target.Z() : target.N();
  xsec *= NNucl;

  return xsec;
}
//____________________________________________________________________________
double PaisQELLambdaPXSec::MHyperon(const Interaction * interaction) const
{
  const XclsTag & xcls = interaction->ExclTag();

  int pdgc  = xcls.StrangeHadronPdg();
  double MR = PDGLibrary::Instance()->Find(pdgc)->Mass();
  return MR;
}
//____________________________________________________________________________
double PaisQELLambdaPXSec::Integral(const Interaction * interaction) const
{

  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool PaisQELLambdaPXSec::ValidProcess(
                                        const Interaction * interaction) const
{
  // Make sure we are dealing with one of the following channels:
  // v + n --> mu+ + Sigma^{-} 
  // v + p --> mu+ + Lambda^{0}
  // v + p --> mu+ + Sigma^{0}

  if(interaction->TestBit(kISkipProcessChk)) return true;

  const XclsTag &      xcls       = interaction->ExclTag();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  bool is_exclusive_strange = (xcls.IsStrangeEvent() && !xcls.IsInclusiveStrange());
  if(!is_exclusive_strange) return false;

  if(!proc_info.IsQuasiElastic()) return false;
  if(!proc_info.IsWeak())         return false;

  bool isP = pdg::IsProton ( init_state.Tgt().HitNucPdg() );
  bool isN = pdg::IsNeutron( init_state.Tgt().HitNucPdg() );

  int pdgc = xcls.StrangeHadronPdg();

  bool can_handle = (
     (pdgc == kPdgSigmaM && isN) ||   /* v + n -> l + #Sigma^{-} */
     (pdgc == kPdgLambda && isP) ||  /* v + p -> l + #Lambda^{0} */
     (pdgc == kPdgSigma0  && isP)   /* v + p -> l + #Sigma^{0}  */
  );

  return can_handle;
}
//____________________________________________________________________________
bool PaisQELLambdaPXSec::ValidKinematics(
                                        const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state  = interaction->InitState();
  double E = init_state.ProbeE(kRfHitNucRest);

  //resonance, final state primary lepton & nucleon mass
  double MR    = this -> MHyperon  (interaction);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Mnuc  = init_state.Tgt().HitNucP4Ptr()->M();
  double Mnuc2 = TMath::Power(Mnuc,2);

  //resonance threshold
  double ER = ( TMath::Power(MR+ml,2) - Mnuc2 ) / (2*Mnuc);

  if(E <= ER) return false;

  return true;
}
//____________________________________________________________________________
void PaisQELLambdaPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PaisQELLambdaPXSec::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void PaisQELLambdaPXSec::LoadConfig(void)
{

  double thc ;
  GetParam( "CabibboAngle", thc ) ;
  fSin8c2 = TMath::Power(TMath::Sin(thc), 2);

  // load QEL form factors model
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
                                             this->SubAlg("FormFactorsAlg"));
  assert(fFormFactorsModel);
  fFormFactors.SetModel(fFormFactorsModel); // <-- attach algorithm

  // load XSec Integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);
}
//____________________________________________________________________________
