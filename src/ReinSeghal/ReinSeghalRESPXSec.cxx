//____________________________________________________________________________
/*!

\class    genie::ReinSeghalRESPXSec

\brief    Computes the double differential cross section for production of a
          single baryon resonance according to the \b Rein-Seghal model.


          The computed cross section is the d^2 xsec/ dQ^2 dW \n

          where \n
            \li \c Q^2 : momentum transfer ^ 2
            \li \c W   : invariant mass of the final state hadronic system

          If it is specified (at the external XML configuration) the cross
          section can weighted with the value of the resonance's Breit-Wigner
          distribution at the given W. The Breit-Wigner distribution type can
          be externally specified. \n

          Is a concrete implementation of the XSecAlgorithmI interface.
                              
\ref      D.Rein and L.M.Seghal, Neutrino Excitation of Baryon Resonances
          and Single Pion Production, Ann.Phys.133, 79 (1981)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 05, 2004

____________________________________________________________________________*/

#include <TMath.h>

#include "AlgFactory/AlgFactory.h"
#include "Base/SPPHelicityAmplModelI.h"
#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResParams.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Base/SPPHelicityAmpl.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/ReinSeghalRESPXSec.h"
#include "Utils/KineLimits.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSeghalRESPXSec::ReinSeghalRESPXSec() :
XSecAlgorithmI()
{
  fName = "genie::ReinSeghalRESPXSec";
}
//____________________________________________________________________________
ReinSeghalRESPXSec::ReinSeghalRESPXSec(const char * param_set) :
XSecAlgorithmI(param_set)
{
  fName = "genie::ReinSeghalRESPXSec";

  FindConfig();
}
//____________________________________________________________________________
ReinSeghalRESPXSec::~ReinSeghalRESPXSec()
{

}
//____________________________________________________________________________
double ReinSeghalRESPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalRes", pDEBUG) << *fConfig;
  
  //----- Get scattering & init-state parameters

  LOG("ReinSeghalRes", pDEBUG)
                        << "Getting scattering and initial-state parameters";

  const ScatteringParams & sc_params  = interaction -> GetScatteringParams();
  const InitialState &     init_state = interaction -> GetInitialState();

  TLorentzVector * p4 = init_state.GetProbeP4(kRfStruckNucAtRest);

  double E    = p4->Energy();
  double W    = sc_params.W();
  double q2   = sc_params.q2();
  double Mnuc = kNucleonMass; // or init_state.TargetMass(); ?
  //double ml    = interaction->GetFSPrimaryLepton()->Mass();

  delete p4;

  //-- Check energy threshold & kinematical limits in q2, W

  double EvThr = kine_limits::EnergyThreshold(interaction);

  if(E <= EvThr) {
    LOG("ReinSeghalRes", pINFO) << "E  = " << E << " < Ethr = " << EvThr;
    return 0;
  }

  //-- Check against physical range in W and Q2
  
  Range1D_t rW  = kine_limits::WRange(interaction); 
  Range1D_t rQ2 = kine_limits::Q2Range_W(interaction); 

  bool in_physical_range = math_utils::IsWithinLimits(W, rW)
                                   && math_utils::IsWithinLimits(-q2, rQ2);

  if(!in_physical_range) return 0;

  //-- Get the input baryon resonance

  Resonance_t resonance = res_utils::FromInteraction(interaction);

  //-- Instantiate a Baryon Resonance Params object & attach data-set

  const Algorithm * bds_alg_base = this->SubAlg(
             "baryonres-dataset-alg-name", "baryonres-dataset-param-set");

  const BaryonResDataSetI * bres_dataset = 
                    dynamic_cast<const BaryonResDataSetI *> (bds_alg_base);

  BaryonResParams bres_params;

  bres_params.SetDataSet(bres_dataset); // <-- attach baryon res. data set
  bres_params.RetrieveData(resonance);  // <-- get parameters for input res.

  //-- Get the resonance mass

  double Mres = bres_params.Mass();

  LOG("ReinSeghalRes", pDEBUG)
        << "Resonance = " << res_utils::AsString(resonance)
                                      << " with mass = " << Mres << " GeV";

  //-- Compute auxiliary & kinematical factors for the Rein-Seghal model
    
  double W2     = TMath::Power(W, 2);
  double Mnuc2  = TMath::Power(Mnuc, 2);
  double k      = 0.5 * (W2 - Mnuc2)/Mnuc;
  double v      = k - 0.5 * q2/Mnuc;
  double v2     = TMath::Power(v, 2);
  double Q2     = v2 - q2;
  double Q      = TMath::Sqrt(Q2);
  double Gf     = TMath::Power(kGF, 2) / ( 4 * TMath::Power(kPi, 2) );
  double Wf     = (-q2/Q2) * (W/Mnuc) * k;
  double Eprime = E - v;
  double U      = 0.5 * (E + Eprime + Q) / E;
  double V      = 0.5 * (E + Eprime - Q) / E;
  double U2     = U*U;
  double V2     = V*V;
  double UV     = U*V;

  //-- Calculate the Helicity Amplitudes

  LOG("ReinSeghalRes", pDEBUG) << "Computing SPP Helicity Amplitudes";

  const Algorithm * hamp_alg_base = this->SubAlg(
             "helicity-amplitudes-alg-name", "helicity-amplitudes-param-set");

  const SPPHelicityAmplModelI * helicity_ampl_model =
                  dynamic_cast<const SPPHelicityAmplModelI *> (hamp_alg_base);

  SPPHelicityAmpl helicity_ampl;

  helicity_ampl.SetModel(helicity_ampl_model); // <-- attach algorithm
  helicity_ampl.Calculate(interaction);        // <-- calculate

  LOG("ReinSeghalRes", pDEBUG)
         << "\n Helicity Amplitudes for ["
                 << res_utils::AsString(resonance) << "]: " << helicity_ampl;

  double amp_minus_1  = helicity_ampl.AmpMinus1 ();
  double amp_plus_1   = helicity_ampl.AmpPlus1  ();
  double amp_minus_3  = helicity_ampl.AmpMinus3 ();
  double amp_plus_3   = helicity_ampl.AmpPlus3  ();
  double amp_0_minus  = helicity_ampl.Amp0Minus ();
  double amp_0_plus   = helicity_ampl.Amp0Plus  ();
  
  //-- Calculate Helicity Cross Sections

  double xsec_left   = TMath::Power( amp_plus_3,  2. ) +
                       TMath::Power( amp_plus_1,  2. );
  double xsec_right  = TMath::Power( amp_minus_3, 2. ) +
                       TMath::Power( amp_minus_1, 2. );
  double xsec_scalar = TMath::Power( amp_0_plus,  2. ) +
                       TMath::Power( amp_0_minus, 2. );

  double scale_lr = 0.5*(kPi/k)*(Mres/Mnuc);
  double scale_sc = 0.5*(kPi/k)*(Mnuc/Mres);

  xsec_left   *=  scale_lr;
  xsec_right  *=  scale_lr;
  xsec_scalar *= (scale_sc*(-Q2/q2));

  LOG("ReinSeghalRes", pDEBUG)
      << "\n Helicity XSecs for ["<< res_utils::AsString(resonance) << "]: "
      << "\n   Sigma-Left   = " << xsec_left
      << "\n   Sigma-Right  = " << xsec_right
      << "\n   Sigma-Scalar = " << xsec_scalar;
  
  double xsec = 0;

  if ( pdg::IsNeutrino(init_state.GetProbePDGCode()) ) {

     xsec = Gf * Wf * ( U2 * xsec_left + V2 * xsec_right + 2*UV*xsec_scalar );

  } else if (pdg::IsAntiNeutrino(init_state.GetProbePDGCode()) ) {

     xsec = Gf * Wf * ( V2 * xsec_left + U2 * xsec_right + 2*UV*xsec_scalar );

  } else {
     LOG("ReinSeghalRes", pERROR) << "Probe is not (anti-)neutrino!";
     return 0;
  }

  //-- Check whether the cross section is to be weighted with a
  //   Breit-Wigner distribution (default: true)

  double bw = 1.0;

  bool weight_bw = fConfig->Exists("weight-with-breit-wigner") ?
                         fConfig->GetBool("weight-with-breit-wigner") : true;
  
  if(weight_bw) {
    
     //-- Get the Breit-Wigner weight for the current resonance
     //   (the Breit-Wigner distribution must be normalized)

     const BreitWignerI * bw_model = this->BreitWignerAlgorithm(resonance);

     bw = bw_model->Eval(W);
  } else {
      LOG("ReinSeghalRes", pINFO) << "Breit-Wigner weighting is turned-off";
  }

  double wxsec = bw * xsec; // weighted-xsec

  SLOG("ReinSeghalRes", pDEBUG)
      << "Res[" << res_utils::AsString(resonance) << "]: "
        << "<Breit-Wigner(=" << bw << ")> * <d^2 xsec/dQ^2 dW [W=" << W
          << ", q2=" << q2 << ", E=" << E << "](="<< xsec << ")> = " << wxsec;

  return wxsec;
}
//____________________________________________________________________________
const BreitWignerI * ReinSeghalRESPXSec::BreitWignerAlgorithm(
                                               Resonance_t resonance_id) const
{
// Retrieves the specifed breit-wigner algorithm

  assert(
       fConfig->Exists("breit-wigner-alg-name")         &&
       fConfig->Exists("breit-wigner-param-set-suffix")
  );

  string alg_name       = fConfig->GetString("breit-wigner-alg-name");
  string param_set_sufx = fConfig->GetString("breit-wigner-param-set-suffix");

  string resonance_name = res_utils::AsString(resonance_id);
  string param_set      = resonance_name + "-" + param_set_sufx;

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const BreitWignerI * bw = dynamic_cast<const BreitWignerI *> (algbase);

  assert(bw);

  return bw;
}
//____________________________________________________________________________
