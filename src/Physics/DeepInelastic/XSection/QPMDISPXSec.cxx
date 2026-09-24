//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <sstream>

#include <TMath.h>
#include <TH1D.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Units.h"
#include "Physics/DeepInelastic/XSection/DISStructureFuncModelI.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/DeepInelastic/XSection/QPMDISPXSec.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"
#include "Physics/Common/PrimaryLeptonUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::constants;
//using namespace genie::units;

//____________________________________________________________________________
QPMDISPXSec::QPMDISPXSec() :
XSecAlgorithmI("genie::QPMDISPXSec")
{
  fInInitPhase = true;
}
//____________________________________________________________________________
QPMDISPXSec::QPMDISPXSec(string config) :
XSecAlgorithmI("genie::QPMDISPXSec", config)
{
  fInInitPhase = true;
}
//____________________________________________________________________________
QPMDISPXSec::~QPMDISPXSec()
{

}
//____________________________________________________________________________
double QPMDISPXSec::XSec(
     const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  // Get kinematical & init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();

  double E     = init_state.ProbeE(kRfHitNucRest);
  double ml    = interaction->FSPrimLepton()->Mass();
  double Mnuc  = init_state.Tgt().HitNucMass();
  double x     = kinematics.x();
  double y     = kinematics.y();

  double E2    = E    * E;
  double ml2   = ml   * ml;
  double ml4   = ml2  * ml2;
  double Mnuc2 = Mnuc * Mnuc;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pDEBUG)
   << "Computing d2xsec/dxdy @ E = " << E << ", x = " << x << ", y = " << y;
#endif

  // One of the xsec terms changes sign for antineutrinos @ DIS/CC

  bool is_nubar_cc = pdg::IsAntiNeutrino(init_state.ProbePdg()) &&
                     proc_info.IsWeakCC();
  int sign = (is_nubar_cc) ? -1 : 1;

  // Calculate the DIS structure functions
  fDISSF.Calculate(interaction);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pDEBUG) << fDISSF;
#endif

  //
  // Compute the differential cross section
  //

  double g2 = kGF2;
  // For EM interaction replace  G_{Fermi} with :
  // a_{em} * pi / ( sqrt(2) * sin^2(theta_weinberg) * Mass_{W}^2 }
  // See C.Quigg, Gauge Theories of the Strong, Weak and E/M Interactions,
  // ISBN 0-8053-6021-2, p.112 (6.3.57)
  // Also, take int account that the photon propagator is 1/p^2 but the
  // W propagator is 1/(p^2-Mass_{W}^2), so weight the EM case with
  // Mass_{W}^4 / q^4
  // So, overall:
  // G_{Fermi}^2 --> a_{em}^2 * pi^2 / (2 * sin^4(theta_weinberg) * q^{4})
  //
  double Q2 = utils::kinematics::XYtoQ2(E,Mnuc,x,y);
  double Q4 = Q2*Q2;
  if(proc_info.IsEM()) {
    g2 = kAem2 * kPi2 / (2.0 * fSin48w * Q4);
  }
  if (proc_info.IsWeakCC()) {
    g2 = kGF2 * kMw2 * kMw2 / TMath::Power((Q2 + kMw2), 2);
  } else if (proc_info.IsWeakNC()) {
    g2 = kGF2 * kMz2 * kMz2 / TMath::Power((Q2 + kMz2), 2);
  }
  double front_factor = (g2*Mnuc*E) / kPi;

  // Build all dxsec/dxdy terms
  double term1 = y * ( x*y + ml2/(2*E*Mnuc) );
  double term2 = 1 - y - Mnuc*x*y/(2*E) - ml2/(4*E2);
  double term3 = sign * (x*y*(1-y/2) - y*ml2/(4*Mnuc*E));
  double term4 = x*y*ml2/(2*Mnuc*E) + ml4/(4*Mnuc2*E2);
  double term5 = -1.*ml2/(2*Mnuc*E);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pDEBUG)
    << "\nd2xsec/dxdy ~ (" << term1 << ")*F1+(" << term2 << ")*F2+("
                  << term3 << ")*F3+(" << term4 << ")*F4+(" << term5 << ")*F5";
#endif

  term1 *= fDISSF.F1();
  term2 *= fDISSF.F2();
  term3 *= fDISSF.F3();
  term4 *= fDISSF.F4();
  term5 *= fDISSF.F5();

  double xsec = front_factor * (term1 + term2 + term3 + term4 + term5);
  xsec = TMath::Max(xsec,0.);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pINFO)
        << "d2xsec/dxdy[FreeN] (E= " << E
                      << ", x= " << x << ", y= " << y << ") = " << xsec;
#endif

   // The algorithm computes d^2xsec/dxdy
  // Check whether variable tranformation is needed
  if(kps!=kPSxyfE) {
    double J = utils::kinematics::Jacobian(interaction,kPSxyfE,kps);
    xsec *= J;
  }

  // If requested return the free nucleon xsec even for input nuclear tgt
  if( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  // Compute nuclear cross section (simple scaling here, corrections must
  // have been included in the structure functions)
  const Target & target = init_state.Tgt();
  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();
  xsec *= NNucl;

  // Apply scaling / if required to reach well known asymmptotic value
  if( proc_info.IsWeakCC() )  xsec *= fCCScale;
  else if( proc_info.IsWeakNC() )  xsec *= fEMScale;
  else if( proc_info.IsEM() )  xsec *= fEMScale;

  // Subtract the inclusive charm production cross section
  interaction->ExclTagPtr()->SetCharm();
  double xsec_charm = fCharmProdModel->XSec(interaction,kps);
  interaction->ExclTagPtr()->UnsetCharm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISPXSec", pINFO)
       << "Subtracting charm piece: " << xsec_charm << " / out of " << xsec;
#endif
  xsec = TMath::Max(0., xsec-xsec_charm);

  // Calculate the DIS structure functions again, but for the whole nucleon (rather than quark)
  // Do this by unsetting the hit quark (reset it afterwards)
  Target * targetPtr = init_state.TgtPtr(); //TODO merge with const ref target above?
  int qpdg = targetPtr->HitQrkPdg();
  targetPtr->UnsetHitQrkPdg();
  fDISSFNucleon.Calculate(interaction);
  targetPtr->SetHitQrkPdg(qpdg);

  return xsec;
}
//____________________________________________________________________________
double QPMDISPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool QPMDISPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const ProcessInfo & proc_info  = interaction->ProcInfo();
  if(!proc_info.IsDeepInelastic()) return false;

  const InitialState & init_state = interaction -> InitState();
  int probe_pdg = init_state.ProbePdg();
  if(!pdg::IsLepton(probe_pdg)) return false;

  if(! init_state.Tgt().HitNucIsSet()) return false;

  int hitnuc_pdg = init_state.Tgt().HitNucPdg();
  if(!pdg::IsNeutronOrProton(hitnuc_pdg)) return false;

  return true;
}
//____________________________________________________________________________
void QPMDISPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDISPXSec::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r( "QPMDISPXSec_specific", false ) ;

  RgKey xdefkey = "XSecModel@genie::EventGenerator/DIS-CC-CHARM";
  RgKey local_key = "CharmXSec" ;
  r.Set( local_key, AlgConfigPool::Instance() -> GlobalParameterList() -> GetAlg(xdefkey) ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void QPMDISPXSec::LoadConfig(void)
{
  // Access global defaults to use in case of missing parameters

  fDISSFModel = 0;
  fDISSFModel =
     dynamic_cast<const DISStructureFuncModelI *> (this->SubAlg("SFAlg"));
  assert(fDISSFModel);

  fDISSF.SetModel(fDISSFModel); // <-- attach algorithm

 // Also init the "nucleon-level" structure function calculation
  fDISSFNucleon.SetModel(fDISSFModel); // <-- attach algorithm

  // Cross section scaling factor
  GetParam( "DIS-CC-XSecScale", fCCScale ) ;
  GetParam( "DIS-NC-XSecScale", fNCScale ) ;
  GetParam( "DIS-EM-XSecScale", fEMScale ) ;

  // sin^4(theta_weinberg)
  double thw  ;
  GetParam( "WeinbergAngle", thw ) ;
  fSin48w = TMath::Power( TMath::Sin(thw), 4 );


  // Since this method would be called every time the current algorithm is
  // reconfigured at run-time, remove all the data cached by this algorithm
  // since they depend on the previous configuration

  if(!fInInitPhase) {
     Cache * cache = Cache::Instance();
     string keysubstr = this->Id().Key() + "/DIS-RES-Join";
     cache->RmMatchedCacheBranches(keysubstr);
  }
  fInInitPhase = false;

  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  // Load the charm production cross section model
  RgKey local_key = "CharmXSec" ;
  RgAlg xalg;
  GetParam( local_key, xalg ) ;
  LOG("DISXSec", pDEBUG)
     << "Loading the cross section model: " << xalg;

  fCharmProdModel = dynamic_cast<const XSecAlgorithmI *> ( this -> SubAlg(local_key) ) ;
  assert(fCharmProdModel);
}
//____________________________________________________________________________
TVector3 QPMDISPXSec::FinalLeptonPolarization(const Interaction* interaction) const
{
  /*
    Compute the final state lepton polarization for this interaction.

    References:
      [1] https://arxiv.org/pdf/hep-ph/0305324
  */

  //
  // Get event information
  //

  const Kinematics & kinematics = interaction->Kine();
  const InitialState & init_state = interaction->InitState();
  const ProcessInfo & proc_info = interaction->ProcInfo();
  const XclsTag & xcls = interaction->ExclTag();

  // Bail for NC
  if (!proc_info.IsWeakCC()) {
    return TVector3(0., 0., 0.);
  }  

  // Check no charm events end up here (would need to apply a correction to W1 if so)
  bool charm = xcls.IsCharmEvent();
  assert(("QPMDISPXSec::FinalLeptonPolarization does not support charm", !charm));

  // Get target nucleon (lab frame)
  const Target & target = init_state.Tgt(); // This is the nucelus
  const TLorentzVector nucleon_p4_lab = target.HitNucP4(); // This is the nucleon

  // Get neutrino (lab frame)
  TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
  TLorentzVector nu_4p_lab = *tempNeutrino; //TODO why this temp object?
  delete tempNeutrino;
  int nu_pdg = init_state.ProbePdg();

  // Get final state lepton (lab frame)
  const TLorentzVector lepton_4p_lab = kinematics.FSLeptonP4();


  //
  // Boost to target nucleon rest frame
  //

  // Polarization calculation is performed in target rest frame
  // Target nucleon has small momentum (Fermi motion) so is not precisely at rest, so transform 
  // to the nucelon's rest frame to perform the polarization calculation correctly.

  // Get beta corresponding to nucleon target
  TVector3 beta = nucleon_p4_lab.BoostVector();

  // Now transform the relevent 4-momenta
  TLorentzVector nucleon_p4_rest(nucleon_p4_lab);
  TLorentzVector nu_4p_rest(nu_4p_lab);
  TLorentzVector lepton_4p_rest(lepton_4p_lab);
  nucleon_p4_rest.Boost(-beta);
  nu_4p_rest.Boost(-beta);
  lepton_4p_rest.Boost(-beta);


  //
  // Get kinematic variables
  //
  
  // Note that symbols used here match [1]

  // Get Ferynman diagram definition, in the target rest frame
  TLorentzVector p = nucleon_p4_rest;
  TLorentzVector k = nu_4p_rest;
  TLorentzVector kprime = lepton_4p_rest;
  TLorentzVector q = k - kprime; //[1] eqn 5

  // Get other kinematic variables
  double Q2 = -q.Mag2();  //[1] eqn 5 //-q**2;
  double p_dot_q = p.Dot(q); // Used in multiple places, so calculating once now
  double x = Q2 / (2. * p_dot_q); // [1] eqn 10
  double M = nucleon_p4_lab.M();

  // Cross-check Q2 and x against the kinematics object         //TODO REMOVE THIS?
  double tol = 1e-3;
  assert(("Q2 mismatch", (Q2 - kinematics.Q2(true)) < tol));
  assert(("x mismatch", (x - kinematics.x(true)) < tol));


  //
  // Calculate W1-5
  //

  // Get F1-5 (nucleon level, as used by [1])
  double F1 = fDISSFNucleon.F1();
  double F2 = fDISSFNucleon.F2();
  double F3 = fDISSFNucleon.F3();
  double F4 = fDISSFNucleon.F4();
  double F5 = fDISSFNucleon.F5();

  // Get W2-5, [1] eqn 53.
  double W_common_term = pow(M, 2) / p_dot_q;
  double W2 = W_common_term * F2;
  double W3 = W_common_term * F3;
  double W4 = W_common_term * F4;
  double W5 = W_common_term * F5;

  // Get W1, which is a special case, see [1] eqn 55.
  // Includes a correction that is applied to the Björken x variable when a charm quark
  // is produced, see the last paragraph of p. 11 in [1].
  double xi = x;
  // if(charm) {
  //   xi = x / (Q2 / (Q2 + pow(m_charm, 2)));   //TODO should I handle charm in here?
  // }
  double W1 = ( 1 + (xi * W_common_term) ) * F1;

  // W6 = 0 in the Standard Model
  double W6 = 0.;


  //
  // Calculate lepton polarization
  //

  TVector3 polarization;
  genie::utils::CalculatePolarizationVectorInTargetRestFrame(
    polarization,
    nu_4p_rest,
    lepton_4p_rest, 
    pdg::IsNeutrino(nu_pdg),
    M,
    W1,
    W2,
    W3,
    W4,
    W5,
    W6
  );

  return polarization;

}
// ____________________________________________________________________________
