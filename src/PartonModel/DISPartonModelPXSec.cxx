//____________________________________________________________________________
/*!

\class    genie::DISPartonModelPXSec

\brief    Massless Parton Model DIS Partial (d^2xsec/dxdy) Cross Section

\ref

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

*/
//____________________________________________________________________________

#include "AlgFactory/AlgFactory.h"
#include "Base/DISFormFactors.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PartonModel/DISPartonModelPXSec.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/KineLimits.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
DISPartonModelPXSec::DISPartonModelPXSec() :
XSecAlgorithmI()
{
  fName = "genie::DISPartonModelPXSec";
}
//____________________________________________________________________________
DISPartonModelPXSec::DISPartonModelPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::DISPartonModelPXSec";

  FindConfig();
}
//____________________________________________________________________________
DISPartonModelPXSec::~DISPartonModelPXSec()
{

}
//____________________________________________________________________________
double DISPartonModelPXSec::XSec(const Interaction * interaction) const
{
  //----- Get scattering & init-state parameters 

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();
  
  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double ml   = interaction->GetFSPrimaryLepton()->Mass();
  double Mnuc = init_state.GetTarget().StruckNucleonMass();
  double E    = p4->Energy();
  double x    = sc_params.x();
  double y    = sc_params.y();

  delete p4;

  //----- Make sure everything makes sense
  
  assert(x>0 && x<1);
  assert(y>0 && y<1);
  assert(E>0);

  //----- Compute final state invariant mass (W) and momentum transfer (Q^2)
  
  double Mnuc2 = pow(Mnuc, 2);
  double W2    = Mnuc2 + 2*Mnuc*E*y*(1-x);
  double W     = TMath::Sqrt(W2);
  double Q2    = 2*Mnuc*E*x*y;

  //----- Get the physical W and Q2 range and check whether the current W,Q2
  //      pair is allowed

  Range1D_t rW  = kine_limits::WRange     (interaction);
  Range1D_t rQ2 = kine_limits::Q2Range_xy (interaction);

  bool in_range = math_utils::IsWithinLimits(Q2, rQ2)
                                        && math_utils::IsWithinLimits(W, rW);

  if(!in_range) {
       LOG("DISXSec", pDEBUG)
             << "\n *** point (W = " << W
                           << ", Q2 = " << Q2 << " is not in physical range";
       LOG("DISXSec", pDEBUG)
             << "\n Physical W range: "
                               << "[" << rW.min << ", " << rW.max << "] GeV";
       LOG("DISXSec", pDEBUG)
             << "\n Physical Q2 range: "
                           << "[" << rQ2.min << ", " << rQ2.max << "] GeV^2";
       return 0;
   }                           

  //----- One of the xsec terms changes sign for antineutrinos
  
  int sign = 1;
  if( pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) sign = -1;

  //----- Calculate the DIS form factors
    
  //-- get the specified DISFormFactorsModelI algorithm

  const DISFormFactorsModelI *
                      form_factors_model = this->FormFactorsAlgorithm();

  //-- instantiate a DISFormFactors object

  DISFormFactors form_factors;

  form_factors.SetModel(form_factors_model); // <-- attach algorithm
  form_factors.Calculate(interaction);       // <-- calculate

  double F1   = form_factors.xF1() / x;
  double F2   = form_factors.F2();
  double F3   = form_factors.xF3() / x;
  double F4   = form_factors.F4();
  double F5   = form_factors.xF5() / x;

  //-- calculate auxiliary parameters
  
  double ml2     = ml    * ml;
  double ml4     = ml2   * ml2;
  double E2      = E     * E;
  double Gfactor = (kGF*kGF*Mnuc*E) / kPi;

  //----- Build all dsigmaQE / dQ2 terms

  double term1 = y * ( x*y + ml2/(2*E*Mnuc) );
  
  double term2 = 1 - y - Mnuc*x*y/(2*E) - ml2/(4*E2);
  
  double term3 = x*y*(1-y/2) - y*ml2/(4*Mnuc*E);
  
  double term4 = x*y*ml2/(2*Mnuc*E) + ml4/(4*Mnuc2*E2);
  
  double term5 = ml2/(2*Mnuc*E);

  //----- Compute the differential cross section 

  double CrossSection = Gfactor*( term1*F1 + term2*F2 +
                                  sign*term3*F3 + term4*F4 - term5*F5 );

  LOG("PartonModel", pDEBUG)  << form_factors;
  LOG("PartonModel", pDEBUG) 
      << "d^2xsec/dxdy (E = " << E << ", x = " << x << ", y = " << y << ") = " 
      << CrossSection;

  return CrossSection;
}
//____________________________________________________________________________
const DISFormFactorsModelI *
                        DISPartonModelPXSec::FormFactorsAlgorithm(void) const
{
  assert (
        fConfig->Exists("form-factors-alg-name")  &&
                                   fConfig->Exists("form-factors-param-set")
  );

  string alg_name, param_set;

  fConfig->Get("form-factors-alg-name",  alg_name  );
  fConfig->Get("form-factors-param-set", param_set );

  //-- get the specified DISFormFactorsModelI algorithm from the AlgFactory

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const DISFormFactorsModelI * form_factors_model =
                      dynamic_cast<const DISFormFactorsModelI *> (algbase);

  assert(form_factors_model);

  return form_factors_model;
}
//____________________________________________________________________________
