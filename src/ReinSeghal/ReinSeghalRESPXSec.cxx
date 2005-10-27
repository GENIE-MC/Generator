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

#include "Algorithm/AlgFactory.h"
#include "BaryonResonance/BaryonResDataSetI.h"
#include "BaryonResonance/BaryonResParams.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/RefFrame.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "ReinSeghal/ReinSeghalRESPXSec.h"
#include "ReinSeghal/RSHelicityAmplModelI.h"
#include "ReinSeghal/RSHelicityAmpl.h"
#include "Utils/KineUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/Range1.h"

using namespace genie;
using namespace genie::constants;

//____________________________________________________________________________
ReinSeghalRESPXSec::ReinSeghalRESPXSec() :
XSecAlgorithmI("genie::ReinSeghalRESPXSec")
{

}
//____________________________________________________________________________
ReinSeghalRESPXSec::ReinSeghalRESPXSec(string config) :
XSecAlgorithmI("genie::ReinSeghalRESPXSec", config)
{

}
//____________________________________________________________________________
ReinSeghalRESPXSec::~ReinSeghalRESPXSec()
{

}
//____________________________________________________________________________
double ReinSeghalRESPXSec::XSec(const Interaction * interaction) const
{
  LOG("ReinSeghalRes", pDEBUG) << *fConfig;

  //----- Get kinematical & init-state parameters

  LOG("ReinSeghalRes", pDEBUG)
                        << "Getting scattering and initial-state parameters";

  const Kinematics &   kinematics = interaction -> GetKinematics();
  const InitialState & init_state = interaction -> GetInitialState();

  double E    = init_state.GetProbeE(kRfStruckNucAtRest);
  double W    = kinematics.W();
  double q2   = kinematics.q2();
  double Mnuc = kNucleonMass; // or init_state.TargetMass(); ?

  //-- Check energy threshold & kinematical limits in q2, W

  double EvThr = utils::kinematics::EnergyThreshold(interaction);

  if(E <= EvThr) {
    LOG("ReinSeghalRes", pINFO) << "E  = " << E << " < Ethr = " << EvThr;
    return 0;
  }

  //-- Check against physical range in W and Q2

  Range1D_t rW  = utils::kinematics::WRange(interaction);
  Range1D_t rQ2 = utils::kinematics::Q2Range_W(interaction);

  bool in_physical_range = utils::math::IsWithinLimits(W, rW)
                                   && utils::math::IsWithinLimits(-q2, rQ2);

  if(!in_physical_range) return 0;

  //-- Get the input baryon resonance

  assert(interaction->GetExclusiveTag().KnownResonance());
  Resonance_t resonance = interaction->GetExclusiveTag().Resonance();

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
        << "Resonance = " << utils::res::AsString(resonance)
                                      << " with mass = " << Mres << " GeV";

  //-- Compute auxiliary & kinematical factors for the Rein-Seghal model

  double W2     = TMath::Power(W, 2);
  double Mnuc2  = TMath::Power(Mnuc, 2);
  double k      = 0.5 * (W2 - Mnuc2)/Mnuc;
  double v      = k - 0.5 * q2/Mnuc;
  double v2     = TMath::Power(v, 2);
  double Q2     = v2 - q2;
  double Q      = TMath::Sqrt(Q2);
  double Gf     = kGF_2 / ( 4 * kPi_2 );
  double Wf     = (-q2/Q2) * (W/Mnuc) * k;
  double Eprime = E - v;
  double U      = 0.5 * (E + Eprime + Q) / E;
  double V      = 0.5 * (E + Eprime - Q) / E;
  double U2     = TMath::Power(U, 2);
  double V2     = TMath::Power(V, 2);
  double UV     = U*V;

  //-- Read in config options [or sets defaults]

  double zeta  = fConfig->Exists("Zeta") ? fConfig->GetDouble("Zeta") : kZeta;
  double omega = fConfig->Exists("Zeta") ? fConfig->GetDouble("Omega"): kOmega;
  double ma2   = fConfig->Exists("Ma2")  ? fConfig->GetDouble("Ma2")  : kResMa2;
  double mv2   = fConfig->Exists("Mv2")  ? fConfig->GetDouble("Mv2")  : kResMv2;

  //-- Calculate the Feynman-Kislinger-Ravndall parameters

  LOG("ReinSeghalRes", pDEBUG) << "Computing the FKR parameters";

  FKR fkr;
  fkr.SetZeta(zeta);
  fkr.SetOmega(omega);
  fkr.SetMa2(ma2);
  fkr.SetMv2(mv2);
  fkr.SetResParams(bres_params);
  fkr.Calculate(interaction);

  LOG("ReinSeghalRes", pDEBUG) << "\n FKR params for ["
                          << utils::res::AsString(resonance) << "]: " << fkr;

  //-- Calculate the Rein-Seghal Helicity Amplitudes

  LOG("ReinSeghalRes", pDEBUG) << "Computing SPP Helicity Amplitudes";

  const RSHelicityAmplModelI * helicity_ampl_model =
         dynamic_cast<const RSHelicityAmplModelI *> (this->SubAlg(
          "helicity-amplitudes-alg-name", "helicity-amplitudes-param-set"));

  RSHelicityAmpl helicity_ampl;

  helicity_ampl.SetModel(helicity_ampl_model); // <-- attach algorithm
  helicity_ampl.Calculate(interaction,fkr);    // <-- calculate

  LOG("ReinSeghalRes", pDEBUG)
         << "\n Helicity Amplitudes for ["
               << utils::res::AsString(resonance) << "]: " << helicity_ampl;

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
      << "\n Helicity XSecs for ["<< utils::res::AsString(resonance) << "]: "
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
      << "Res[" << utils::res::AsString(resonance) << "]: "
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
       fConfig->Exists("breit-wigner-alg-name") &&
       fConfig->Exists("breit-wigner-param-set-suffix")
  );

  string alg_name       = fConfig->GetString("breit-wigner-alg-name");
  string param_set_sufx = fConfig->GetString("breit-wigner-param-set-suffix");

  string resonance_name = utils::res::AsString(resonance_id);
  string param_set      = resonance_name + "-" + param_set_sufx;

  AlgFactory * algf = AlgFactory::Instance();

  const Algorithm * algbase = algf->GetAlgorithm(alg_name, param_set);

  const BreitWignerI * bw = dynamic_cast<const BreitWignerI *> (algbase);

  assert(bw);

  return bw;
}
//____________________________________________________________________________
