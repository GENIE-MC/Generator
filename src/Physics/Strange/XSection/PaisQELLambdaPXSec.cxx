//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Hugh Gallagher
 Tufts University
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
using namespace std::complex_literals;

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
  
  // Do precise calculation of lepton polarization
  GetParamDef( "PreciseLeptonPol", fIsPreciseLeptonPolarization, false ) ;

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
const TVector3 & PaisQELLambdaPXSec::FinalLeptonPolarization (const Interaction* interaction) const
{
  if (!fIsPreciseLeptonPolarization) return XSecAlgorithmI::FinalLeptonPolarization(interaction);
  double ml = interaction->FSPrimLepton()->Mass();
  // First we need access to all of the particles in the interaction
  // The particles were stored in the lab frame
  //----- get kinematics & init state - compute auxiliary vars
  const Kinematics &   kinematics  = interaction->Kine();
  const InitialState & init_state  = interaction->InitState();
  const Target &       target      = init_state.Tgt();
  
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  
  // HitNucMass() looks up the PDGLibrary (on-shell) value for the initial
  // struck nucleon
  double Mnuc  = target.HitNucMass();
  double Mnuc2 = Mnuc*Mnuc;
  
  // Note that GetProbeP4 defaults to returning the probe 4-momentum in the
  // struck nucleon rest frame, so we have to explicitly ask for the lab frame
  // here
  TLorentzVector * tempNeutrino = init_state.GetProbeP4(kRfLab);
  TLorentzVector neutrinoMom = *tempNeutrino;
  delete tempNeutrino;
  TLorentzVector inNucleonMom(*init_state.TgtPtr()->HitNucP4Ptr());
  const TLorentzVector leptonMom = kinematics.FSLeptonP4();
  
  // Ordinary 4-momentum transfer
  TLorentzVector qP4 = neutrinoMom - leptonMom;
  double Q2 = -qP4.Mag2();
  interaction->KinePtr()->SetQ2(Q2);
  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);
  double F1V   = fFormFactors.F1V();
  double xiF2V = fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
    
  //neutrino momentum transfer
  double q2  = -Q2;
  
  //----- Calculate the differential cross section dxsec/dQ^2
  //resonance mass & nucleon mass
  double Mi        = Mnuc;
  double Mf        = (this)->MHyperon(interaction);
  // calculate w coefficients
  //start with Mass terms
  double Mp    = Mf + Mi;
  double Mm    = Mf - Mi;
  double Mm2   = Mm*Mm;
  double Mp2   = Mp*Mp;
  double M     = 0.5*Mp;
  double M2    = M*M;
  
  //Powers of Form Factors
  double FA2   = FA*FA;
  
  //Calculate W terms
  double W1 = (Mm2 - q2)*TMath::Sq(F1V + xiF2V)/4/Mnuc2 + (Mp2 - q2)*FA2/4/Mnuc2;
  double W2 = FA2 + TMath::Sq(F1V + xiF2V - Mp*xiF2V/2/Mnuc) - q2*xiF2V*xiF2V/4/Mnuc2;
  double W3 = 2*FA*(F1V + xiF2V);
  double W4(0), W5(0);
  
  
  double p[4], q[4], epq[4][4], k[4], l[4], s[4], eskl[4];
  std::complex<double> jp[4], jm[4];
    
  p[0] = inNucleonMom.E();
  p[1] = inNucleonMom.Px();
  p[2] = inNucleonMom.Py();
  p[3] = inNucleonMom.Pz();
  
  q[0] = qP4.E();
  q[1] = qP4.Px();
  q[2] = qP4.Py();
  q[3] = qP4.Pz();
  
  k[0] = neutrinoMom.E();
  k[1] = -neutrinoMom.Px();
  k[2] = -neutrinoMom.Py();
  k[3] = -neutrinoMom.Pz();
  
  l[0] = leptonMom.E();
  l[1] = -leptonMom.Px();
  l[2] = -leptonMom.Py();
  l[3] = -leptonMom.Pz();
  
  s[0] = leptonMom.P()/ml;
  s[1] = -leptonMom.Vect().Unit().X()*leptonMom.E()/ml;
  s[2] = -leptonMom.Vect().Unit().Y()*leptonMom.E()/ml;
  s[3] = -leptonMom.Vect().Unit().Z()*leptonMom.E()/ml;
  
  // epsilon^\alpha\beta\gamma\delta p_\gamma q_\delta
  for (int a = 0; a < 4; a++)
  {
    for (int b = 0; b < 4; b++)
    {
        epq[a][b] = 0;
        if (b == a) continue;
        for (int g = 0; g < 4; g++)
        {
            if (g == b || g == a) continue;
            for (int d = 0; d < 4; d++)
            {
                if (d == g || d == b || d == a) continue;
                epq[a][b] += e(a,b,g,d)*(a == 0?1:-1)*(b == 0?1:-1)*p[g]*q[d];
            }
        }
    }
  }
  
  // epsilon_\alpha\beta\gamma\delta s^\beta k^\gamma l^\delta
  for (int a = 0; a < 4; a++)
  {
    eskl[a] = 0;
    for (int b = 0; b < 4; b++)
    {
        if (b == a) continue;
        for (int g = 0; g < 4; g++)
        {
            if (g == b || g == a) continue;
            for (int d = 0; d < 4; d++)
            {
                if (d == g || d == b || d == a) continue;
                double sb = s[b]*(b == 0?1:-1);
                double kg = k[g]*(g == 0?1:-1);
                double ld = l[d]*(d == 0?1:-1);
                eskl[a] += e(a,b,g,d)*sb*kg*ld;
            }
        }
    }
  }
    
  double kl = k[0]*l[0] - k[1]*l[1] - k[2]*l[2] - k[3]*l[3];
  double ks = k[0]*s[0] - k[1]*s[1] - k[2]*s[2] - k[3]*s[3];
        
  for (int a = 0; a < 4; a++)
  {
     if (is_neutrino)
     {
        jp[a] =  (l[a]*ks - s[a]*kl - 1i*eskl[a] + ml*k[a])/sqrt(kl + ml*ks);   //jp_\alpha
        jm[a] = (-l[a]*ks + s[a]*kl + 1i*eskl[a] + ml*k[a])/sqrt(kl - ml*ks);   //jm_\alpha
     }
     else
     {
        jp[a] =  (l[a]*ks - s[a]*kl + 1i*eskl[a] - ml*k[a])/sqrt(kl - ml*ks);   //jp_\alpha
        jm[a] =  (l[a]*ks - s[a]*kl + 1i*eskl[a] + ml*k[a])/sqrt(kl + ml*ks);   //jm_\alpha
     }
  }

  //Additional constants and variables
  std::complex<double> Wmunu, Wnumu, LWpp(0, 0), LWpm(0, 0), LWmp(0, 0), LWmm(0, 0);
  for(int mu = 0; mu < 4; mu++)
  {
     for(int nu = mu;nu < 4; nu++)
     {
        double Wreal = -g(mu,nu)*W1 + p[mu]*p[nu]*W2/M2 + q[mu]*q[nu]*W4/M2 + (p[mu]*q[nu] + q[mu]*p[nu])*W5/2/M2;
        double Wimag = epq[mu][nu]*W3/2/M2;
        Wmunu = Wreal - 1i*Wimag;  // W^\mu\nu
        LWpp += jp[mu]*std::conj(jp[nu])*Wmunu; // Lpp_\mu\nu*W^\mu\nu
        LWpm += jp[mu]*std::conj(jm[nu])*Wmunu; // Lpm_\mu\nu*W^\mu\nu
        LWmp += jm[mu]*std::conj(jp[nu])*Wmunu; // Lmp_\mu\nu*W^\mu\nu
        LWmm += jm[mu]*std::conj(jm[nu])*Wmunu; // Lmm_\mu\nu*W^\mu\nu
        if (mu != nu)
        {
            Wnumu = Wreal + 1i*Wimag;
            LWpp += jp[nu]*std::conj(jp[mu])*Wnumu; // Lpp_\mu\nu*W^\mu\nu
            LWpm += jp[nu]*std::conj(jm[mu])*Wnumu; // Lpm_\mu\nu*W^\mu\nu
            LWmp += jm[nu]*std::conj(jp[mu])*Wnumu; // Lmp_\mu\nu*W^\mu\nu
            LWmm += jm[nu]*std::conj(jm[mu])*Wnumu; // Lmm_\mu\nu*W^\mu\nu
        }
    }
  }
  std::complex<double> LWppmm = LWpp + LWmm;
  std::complex<double> rhopp = LWpp/LWppmm;
  std::complex<double> rhopm = LWpm/LWppmm;
  std::complex<double> rhomp = LWmp/LWppmm;
  std::complex<double> rhomm = LWmm/LWppmm;
  double PL = std::real(rhopp - rhomm);
  double PP = std::real(rhopm + rhomp);
  double PT = std::imag(rhomp - rhopm);
  
  TVector3 neutrinoMom3 = neutrinoMom.Vect();                                          
  TVector3 leptonMom3 = leptonMom.Vect();
  TVector3 Pz = leptonMom3.Unit();
  TVector3 Px = neutrinoMom3.Cross(leptonMom3).Unit();
  TVector3 Py = Pz.Cross(Px);
  TVector3 pol = PT*Px + PP*Py + PL*Pz;
  fFinalLeptonPolarization = pol;
  
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  std::cout << "PL = " << PL << ", PT = " << PT << ", PP = " << PP << "\n";
  std::cout << fFinalLeptonPolarization.Mag() << "\n";
  std::cout << "PL@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << std::endl;
  
  return fFinalLeptonPolarization;
}
//____________________________________________________________________________
inline int PaisQELLambdaPXSec::g(int a, int b) const
{
    return (a==b)*(2*(a==0) - 1);
}
//____________________________________________________________________________
inline int PaisQELLambdaPXSec::e(int a, int b, int c, int d) const
{
    return (b - a)*(c - a)*(d - a)*(c - b)*(d - b)*(d - c)/12;
}
//____________________________________________________________________________
