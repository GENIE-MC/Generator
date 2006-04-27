//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 05, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Base/QELFormFactors.h"
#include "Base/QELFormFactorsModelI.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/KineVar.h"
#include "Conventions/Units.h"
#include "LlewellynSmith/QELPXSec.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"
#include "Utils/NuclearUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
QELPXSec::QELPXSec() :
XSecAlgorithmI("genie::QELPXSec")
{

}
//____________________________________________________________________________
QELPXSec::QELPXSec(string config) :
XSecAlgorithmI("genie::QELPXSec", config)
{

}
//____________________________________________________________________________
QELPXSec::~QELPXSec()
{

}
//____________________________________________________________________________
double QELPXSec::XSec(const Interaction * interaction) const
{
  LOG("QELPXSec", pDEBUG) << *fConfig;

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  //----- get kinematics & init-state parameters
  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();
  const Target & target = init_state.GetTarget();

  double E  = init_state.GetProbeE(kRfStruckNucAtRest);
  double E2 = TMath::Power(E,2);
  double ml = interaction->GetFSPrimaryLepton()->Mass();
  double M  = target.StruckNucleonMass();
  double q2 = kinematics.q2();

  //----- one of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.GetProbePDGCode());
  int sign = (is_neutrino) ? -1 : 1;

  //----- calculate the QEL form factors
  fFormFactors.Calculate(interaction);    

  double F1V   = fFormFactors.F1V();
  double xiF2V = fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  double Fp    = fFormFactors.Fp();

  LOG("QELPXSec", pDEBUG) << ENDL << fFormFactors;

  //-- calculate auxiliary parameters
  double ml2     = TMath::Power(ml,    2);
  double M2      = TMath::Power(M,     2);
  double M4      = TMath::Power(M2,    2);
  double FA2     = TMath::Power(FA,    2);
  double Fp2     = TMath::Power(Fp,    2);
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double Gfactor = M2*kGF2*fCos8c2 / (8*kPi*E2);
  double s_u     = 4*E*M + q2 - ml2;
  double q2_M2   = q2/M2;

  //----- compute free nucleon differential cross section
  double A = (0.25*(ml2-q2)/M2) * (
	      (4-q2_M2)*FA2 - (4+q2_M2)*F1V2 - q2_M2*xiF2V2*(1+0.25*q2_M2)
              -4*q2_M2*F1V*xiF2V - (ml2/M2)*( 
               (F1V2+xiF2V2+2*F1V*xiF2V)+(FA2+4*Fp2+4*FA*Fp)+(q2_M2-4)*Fp2));
  double B = -1 * q2_M2 * FA*(F1V+xiF2V);
  double C = 0.25*(FA2 + F1V2 - 0.25*q2_M2*xiF2V2);

  double xsec = Gfactor * (A + sign*B*s_u/M2 + C*s_u*s_u/M4);

  LOG("QELPXSec", pDEBUG)
     << "dXSec[QEL]/dQ2 [FreeN](E = "<< E << ", Q2 = "<< -q2 << ") = "<< xsec;
  LOG("QELPXSec", pDEBUG) 
                 << "A(Q2) = " << A << ", B(Q2) = " << B << ", C(Q2) = " << C;

  //----- if requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- compute nuclear suppression factor
  //      (R(Q2) is adapted from NeuGEN - see comments therein)
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);

  //----- number of scattering centers in the target
  int nucpdgc = target.StruckNucleonPDGCode();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N(); 

  LOG("QELPXSec", pDEBUG) 
       << "Nuclear suppression factor R(Q2) = " << R << ", NNucl = " << NNucl;

  xsec *= (R*NNucl); // nuclear xsec

  LOG("QELPXSec", pDEBUG)
   << "dXSec[QEL]/dQ2 [Nuclear](E = "<< E << ", Q2 = "<< -q2 << ") = "<< xsec;

  return xsec;
}
//____________________________________________________________________________
bool QELPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->GetInitialState();
  const ProcessInfo &  proc_info  = interaction->GetProcessInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.GetTarget().StruckNucleonPDGCode();
  int  nu  = init_state.GetProbePDGCode();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool prcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
bool QELPXSec::ValidKinematics(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipKinematicChk)) return true;

  const InitialState & init_state = interaction -> GetInitialState();
  const Kinematics &   kinematics = interaction -> GetKinematics();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double Ethr = utils::kinematics::EnergyThreshold(interaction);

  LOG("QELPXSec", pDEBUG)
       << "Computing QEL dXSec/dQ2 for Ev = " << E
                            << " / Neutrino Energy Threshold = " << Ethr;
  if(E <= Ethr) {
     LOG("QELPXSec", pINFO) << "Ev = " << E << " <= Ethreshold = "<< Ethr;
     return false;
  }

  double     Q2  = kinematics.Q2();
  Range1D_t  rQ2 = utils::kinematics::KineRange(interaction, kKVQ2);

  LOG("QELPXSec", pDEBUG) << "Q2 integration range = ("
                                    << rQ2.min << ", " << rQ2.max << ")";
  bool in_range = utils::math::IsWithinLimits(Q2, rQ2);
  if(!in_range) {
     LOG("QELPXSec", pDEBUG) << "Q2 = " << Q2 << ", not in allowed range";
     return false;
  }
  return true;
}
//____________________________________________________________________________
void QELPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELPXSec::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();
  
  double thc = fConfig->GetDoubleDef(
                              "cabbibo-angle", gc->GetDouble("CabbiboAngle"));
  fCos8c2 = TMath::Power(TMath::Cos(thc), 2);

  fFormFactorsModel = 0;
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
             this->SubAlg("form-factors-alg-name", "form-factors-param-set"));
  assert(fFormFactorsModel);

  fFormFactors.SetModel(fFormFactorsModel); // <-- attach algorithm
}
//____________________________________________________________________________
