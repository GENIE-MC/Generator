//____________________________________________________________________________
/*!

\class    genie::QELPXSec

\brief    Computes the differential Quasi Elastic cross section dxsec/dq^2.\n
          Is a concrete implementation of the XSecAlgorithmI interface.\n

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Base/QELFormFactors.h"
#include "Base/QELFormFactorsModelI.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Units.h"
#include "CrossSections/QELPXSec.h"
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

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double Mnuc = init_state.GetTarget().StruckNucleonP4()->M();
  double q2   = kinematics.q2();

  //----- one of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.GetProbePDGCode());
  int sign = (is_neutrino) ? 1 : -1;

  //----- calculate the QEL form factors
  QELFormFactors form_factors;
  form_factors.SetModel(fFormFactorsModel); // <-- attach algorithm
  form_factors.Calculate(interaction);      // <-- calculate

  double F1V   = form_factors.F1V();
  double xiF2V = form_factors.xiF2V();
  double FA    = form_factors.FA();
  double Fp    = form_factors.Fp();

  LOG("QELPXSec", pDEBUG) << ENDL << form_factors;

  //-- calculate auxiliary parameters
  double ml2     = ml    * ml;
  double ml4     = ml2   * ml2;
  double Mnuc2   = Mnuc  * Mnuc;
  double Mnuc4   = Mnuc2 * Mnuc2;
  double q4      = q2   * q2;
  double s_u     = 4*E*Mnuc + q2 - ml2;
  double FA_2    = FA    * FA;
  double Fp_2    = Fp    * Fp;
  double F1V_2   = F1V   * F1V;
  double xiF2V_2 = xiF2V * xiF2V;
  double Gfactor = TMath::Power( (kGF*kCos8c)/E, 2.) / (8.*kPi);
  double ml2_q2  = ml2 - q2;

  //----- start building all dsigmaQE / dQ2 terms
  double term1 = F1V_2     * (q4 - 4*Mnuc2*ml2_q2 - ml4)     / (4*Mnuc2);
  double term2 = xiF2V_2   * (4*Mnuc2*(q4-ml4) - q4*ml2_q2 ) / (16*Mnuc4);
  double term3 = FA_2      * (q4 + 4*Mnuc2*ml2_q2 - ml4)     / (4*Mnuc2);
  double term4 = Fp_2      * ml2*q2 * ml2_q2                 / (4*Mnuc4);
  double term5 = F1V*xiF2V * (2*q4 + q2*ml2 + ml4)           / (2*Mnuc2);
  double term6 = FA*Fp     * ml2 * ml2_q2                    / (2*Mnuc2);
  double term7 = sign * FA*(F1V+xiF2V) * q2 * s_u            / Mnuc2;
  double term8 = (F1V_2 - xiF2V*xiF2V*q2/(4*Mnuc2) +FA_2 ) * s_u*s_u / (4*Mnuc2);

  //----- compute free nucleon differential cross section
  double xsec = Gfactor*(term1+term2+term3-term4+term5-term6+term7+term8);

  LOG("QELPXSec", pDEBUG)
     << "dXSec[QEL]/dQ2 [FreeN](E = "<< E << ", Q2 = "<< -q2 << ") = "<< xsec;

  //----- if requested return the free nucleon xsec even for input nuclear tgt 
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  //----- compute nuclear suppression factor
  double R = nuclear::NuclQELXSecSuppression("Default", 0.5, interaction);

  //----- number of scattering centers in the target
  const Target & tgt = init_state.GetTarget();
  int nucpdgc = tgt.StruckNucleonPDGCode();
  int NNucl = (pdg::IsProton(nucpdgc)) ? tgt.Z() : tgt.N(); 

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

  bool ccprcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  bool ncprcok = proc_info.IsWeakNC() && (isP||isN) && (isnu||isnub);
  bool prcok   = ccprcok || ncprcok;
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
  Range1D_t  rQ2 = utils::kinematics::Q2Range_M(interaction);

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
  this->LoadSubAlg();
}
//____________________________________________________________________________
void QELPXSec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadSubAlg();
}
//____________________________________________________________________________
void QELPXSec::LoadSubAlg(void)
{
  fFormFactorsModel = 0;

  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
             this->SubAlg("form-factors-alg-name", "form-factors-param-set"));
  assert(fFormFactorsModel);
}
//____________________________________________________________________________
