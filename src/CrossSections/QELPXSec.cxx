//____________________________________________________________________________
/*!

\class    genie::QELPXSec

\brief    Computes the differential Quasi Elastic cross section dxsec/dq^2.

          Is a concrete implementation of the XSecAlgorithmI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

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

using namespace genie;
using namespace genie::constants;

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

  //----- get scattering & init-state parameters

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();

  TLorentzVector * nu_p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double E    = nu_p4->Energy();
  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double Mnuc = init_state.GetTarget().StruckNucleonMass();
  double q2   = sc_params.q2();

  delete nu_p4;

  //----- "phase space" cuts

  if (E < ml) return 0.;
  if (utils::math::AreEqual( TMath::Abs(q2), 0.)) return 0.;

  //----- one of the xsec terms changes sign for antineutrinos

  int sign = 1;
  if( pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) sign = -1;

  //----- calculate the QEL form factors

  //-- instantiate a QELFormFactors object, attach a model & calculate

  const Algorithm * algbase = this->SubAlg(
                          "form-factors-alg-name", "form-factors-param-set");
  const QELFormFactorsModelI * form_factors_model =
                        dynamic_cast<const QELFormFactorsModelI *> (algbase);

  QELFormFactors form_factors;

  form_factors.SetModel(form_factors_model); // <-- attach algorithm
  form_factors.Calculate(interaction);       // <-- calculate

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
  double Gfactor = pow( (kGF*kCos8c)/E, 2.) / (8.*kPi);
  double ml2_q2  = ml2 - q2;

  //----- start building all dsigmaQE / dQ2 terms

  double term1 = F1V_2     * (q4 - 4*Mnuc2*ml2_q2 - ml4)     / (4*Mnuc2);
  double term2 = xiF2V_2   * (4*Mnuc2*(q4-ml4) - q4*ml2_q2 ) / (16*Mnuc4);
  double term3 = FA_2      * (q4 + 4*Mnuc2*ml2_q2 - ml4)     / (4*Mnuc2);
  double term4 = Fp_2      * ml2*q2 * ml2_q2                 / (4*Mnuc4);
  double term5 = F1V*xiF2V * (2*q4 + q2*ml2 + ml4)           / (2*Mnuc2);
  double term6 = FA*Fp     * ml2 * ml2_q2                    / (2*Mnuc2);
  double term7 = sign * FA*(F1V+xiF2V) * q2 * s_u            / Mnuc2;
  double term8 = ( F1V_2 - xiF2V*xiF2V*q2/(4*Mnuc2) + FA_2 ) * s_u*s_u / (4*Mnuc2);

  //----- compute differential cross section

  double CrossSection = Gfactor*(term1+term2+term3-term4+term5-term6+term7+term8);

  return CrossSection;
}
//____________________________________________________________________________

