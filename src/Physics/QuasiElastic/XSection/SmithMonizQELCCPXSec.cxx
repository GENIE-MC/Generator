//____________________________________________________________________________
/*
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org

  Igor Kakorin <kakorin@jinr.ru>
  Joint Institute for Nuclear Research

  adapted from  fortran code provided by:

  Konstantin Kuzmin <kkuzmin@theor.jinr.ru>
  Joint Institute for Nuclear Research

  Vladimir Lyubushkin
  Joint Institute for Nuclear Research

  Vadim Naumov <vnaumov@theor.jinr.ru>
  Joint Institute for Nuclear Research

  based on code of:
  Costas Andreopoulos <c.andreopoulos \at cern.ch>
  University of Liverpool
*/
//____________________________________________________________________________

#include <sstream>
#include <string>
#include <algorithm>

#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/QuasiElastic/XSection/SmithMonizQELCCPXSec.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using std::ostringstream;

//____________________________________________________________________________
SmithMonizQELCCPXSec::SmithMonizQELCCPXSec() :
XSecAlgorithmI("genie::SmithMonizQELCCPXSec")
{

}
//____________________________________________________________________________
SmithMonizQELCCPXSec::SmithMonizQELCCPXSec(string config) :
XSecAlgorithmI("genie::SmithMonizQELCCPXSec", config)
{

}
//____________________________________________________________________________
SmithMonizQELCCPXSec::~SmithMonizQELCCPXSec()
{

}
//____________________________________________________________________________
double SmithMonizQELCCPXSec::XSec(
   const Interaction * interaction, KinePhaseSpace_t kps) const
{
  double xsec = 0. ;
  // dimension of kine phase space
  std::string s = KinePhaseSpace::AsString(kps);
  int kpsdim = s!="<|E>"?1 + std::count(s.begin(), s.begin()+s.find('}'), ','):0;
  
  if(!this -> ValidProcess (interaction) )
  {
    LOG("SmithMoniz",pWARN) << "not a valid process";
    return 0.;
  }

  if(kpsdim == 1)
  {
     if(! this -> ValidKinematics (interaction) )
     {
        LOG("SmithMoniz",pWARN) << "not valid kinematics";
        return 0.;
     }
     xsec = this->dsQES_dQ2_SM(interaction);
  }

  if(kpsdim == 2) 
  {
    xsec = this->d2sQES_dQ2dv_SM(interaction);
  }
  
  if(kpsdim == 3) 
  {
    xsec = this->d3sQES_dQ2dvdkF_SM(interaction);
  }

  
  // The algorithm computes d^1xsec/dQ2, d^2xsec/dQ2dv or d^3xsec/dQ2dvdp
  // Check whether variable tranformation is needed
  if ( kps != kPSQ2fE && kps != kPSQ2vfE ) 
  {
     double J = 1.;
     if (kpsdim == 1)
       J = utils::kinematics::Jacobian(interaction, kPSQ2fE, kps);
     else if (kpsdim == 2)
       J = utils::kinematics::Jacobian(interaction, kPSQ2vfE, kps);
     else if (kpsdim == 3)
       J = utils::kinematics::Jacobian(interaction, kPSQ2vpfE, kps);
     xsec *= J;
  }

  return xsec;

}
//____________________________________________________________________________
double SmithMonizQELCCPXSec::Integral(const Interaction * in) const
{
  return fXSecIntegrator->Integrate(this,in);

}
//____________________________________________________________________________
bool SmithMonizQELCCPXSec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();

  if(!proc_info.IsQuasiElastic()) return false;

  int  nuc = init_state.Tgt().HitNucPdg();
  int  nu  = init_state.ProbePdg();

  bool isP   = pdg::IsProton(nuc);
  bool isN   = pdg::IsNeutron(nuc);
  bool isnu  = pdg::IsNeutrino(nu);
  bool isnub = pdg::IsAntiNeutrino(nu);

  bool prcok = proc_info.IsWeakCC() && ((isP&&isnub) || (isN&&isnu));
  if(!prcok) return false;

  return true;
}
//____________________________________________________________________________
void SmithMonizQELCCPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SmithMonizQELCCPXSec::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r( "SmithMonizQELCCPXSec_specific", false ) ;
  r.Set("sm_utils_algo", RgAlg("genie::SmithMonizUtils","Default") ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void SmithMonizQELCCPXSec::LoadConfig(void)
{

  // Cross section scaling factor
  GetParamDef( "QEL-CC-XSecScale", fXSecScale, 1. ) ;

  double Vud;
  GetParam( "CKM-Vud", Vud ) ;
  fVud2 = TMath::Power( Vud, 2 );

   // load QEL form factors model
  fFormFactorsModel = dynamic_cast<const QELFormFactorsModelI *> (
                                             this->SubAlg("FormFactorsAlg"));
  assert(fFormFactorsModel);
  fFormFactors.SetModel(fFormFactorsModel); // <-- attach algorithm

   // load XSec Integrators
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

  sm_utils = const_cast<genie::SmithMonizUtils *>(
               dynamic_cast<const genie::SmithMonizUtils *>(
                 this -> SubAlg( "sm_utils_algo" ) ) ) ;

}
//____________________________________________________________________________
double SmithMonizQELCCPXSec::d3sQES_dQ2dvdkF_SM(const Interaction * interaction) const
{
  // Assuming that variables E_nu, Q2, \nu and kF are within allowable kinematic region
  // which are specified in methods: genie::utils::gsl::d2Xsec_dQ2dv::DoEval and QELEventGeneratorSM::ProcessEventRecord
  // Get kinematics & init-state parameters
  const Kinematics &  kinematics = interaction -> Kine();
  sm_utils->SetInteraction(interaction);
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  double E_nu    = init_state.ProbeE(kRfLab);
  double Q2      = kinematics.GetKV(kKVQ2);
  double v       = kinematics.GetKV(kKVv);
  double kF      = kinematics.GetKV(kKVPn);
  double kkF     = kF*kF;
  int nucl_pdg_ini = target.HitNucPdg();
  int nucl_pdg_fin = genie::pdg::SwitchProtonNeutron(nucl_pdg_ini);
  
  PDGLibrary * pdglib = PDGLibrary::Instance();
  TParticlePDG * nucl_fin = pdglib->Find( nucl_pdg_fin );

  double E_BIN   = sm_utils->GetBindingEnergy();
  double m_ini   = target.HitNucMass();
  double mm_ini  = m_ini*m_ini;
  double m_fin   = nucl_fin -> Mass();                         //  Mass of final hadron or hadron system (GeV)
  double mm_fin  = m_fin*m_fin;
  double m_tar   = target.Mass();                              //  Mass of target nucleus (GeV)
  double mm_tar  = m_tar*m_tar;
  
  // One of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int n_NT = (is_neutrino) ? +1 : -1;
  
  double E_p     = TMath::Sqrt(mm_ini+kkF)-E_BIN;
  //|\vec{q}|
  double qqv     = v*v+Q2;
  double qv      = TMath::Sqrt(qqv);
  double cosT_p  = ((v-E_BIN)*(2*E_p+v+E_BIN)-qqv+mm_ini-mm_fin)/(2*kF*qv);           //\cos\theta_p
  if (cosT_p < -1.0 || cosT_p > 1.0 ) 
  {
     return 0.0;
  }
  
  double pF      = TMath::Sqrt(kkF+(2*kF*qv)*cosT_p+qqv);
  
  double E_lep   = E_nu-v;
  double m_lep = interaction->FSPrimLepton()->Mass();
  double mm_lep = m_lep*m_lep;
  if (E_lep < m_lep) 
  {
    return 0.0;
  }
  double P_lep   = TMath::Sqrt(E_lep*E_lep-mm_lep);
  double k6 = (Q2+mm_lep)/(2*E_nu);
  double cosT_lep= (E_lep-k6)/P_lep;
  if (cosT_lep < -1.0 || cosT_lep > 1.0 ) return 0.0;
  
  double cosT_k  = (v+k6)/qv;
  if (cosT_k < -1.0 || cosT_k > 1.0 ) return 0.0;

  double b2_flux = (E_p-kF*cosT_k*cosT_p)*(E_p-kF*cosT_k*cosT_p);
  double c2_flux = kkF*(1-cosT_p*cosT_p)*(1-cosT_k*cosT_k);
  
  double k1 = fVud2*kNucleonMass2*kPi;
  double k2 = mm_lep/(2*mm_tar);
  double k7 = P_lep*cosT_lep;

  double P_Fermi = sm_utils->GetFermiMomentum();
  double FV_SM   = 4.0*TMath::Pi()/3*TMath::Power(P_Fermi, 3);
  double factor  = k1*(m_tar*kF/(FV_SM*qv*TMath::Sqrt(b2_flux-c2_flux)))*SmithMonizUtils::rho(P_Fermi, 0.0, kF)*(1-SmithMonizUtils::rho(P_Fermi, 0.01, pF));

  double a2      = kkF/kNucleonMass2;
  double a3      = a2*cosT_p*cosT_p;
  double a6      = kF*cosT_p/kNucleonMass;
  double a7      = E_p/kNucleonMass;
  double a4      = a7*a7;
  double a5      = 2*a7*a6;

  double k3      = v/qv;
  double k4      = (3*a3-a2)/qqv;
  double k5      = (a7-a6*k3)*m_tar/kNucleonMass;
  
  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);
  double F_V   = fFormFactors.F1V();
  double F_M   = fFormFactors.xiF2V();
  double F_A   = fFormFactors.FA();
  double F_P   = fFormFactors.Fp();
  double FF_V  = F_V*F_V;
  double FF_M  = F_M*F_M;
  double FF_A  = F_A*F_A;

  double t       = Q2/(4*kNucleonMass2);
  double W_1     = FF_A*(1+t)+t*(F_V+F_M)*(F_V+F_M);                     //Ref.[1], \tilde{T}_1
  double W_2     = FF_A+FF_V+t*FF_M;                                     //Ref.[1], \tilde{T}_2
  double W_3     =-2*F_A*(F_V+F_M);                                      //Ref.[1], \tilde{T}_8
  double W_4     =-0.5*F_V*F_M-F_A*F_P+t*F_P*F_P-0.25*(1-t)*FF_M;        //Ref.[1], \tilde{T}_\alpha
  double W_5     = FF_V+t*FF_M+FF_A;

  double T_1     = 1.0*W_1+(a2-a3)*0.5*W_2;                              //Ref.[1], W_1
  double T_2     = ((a2-a3)*Q2/(2*qqv)+a4-k3*(a5-k3*a3))*W_2;            //Ref.[1], W_2
  double T_3     = k5*W_3;                                               //Ref.[1], W_8
  double T_4     = mm_tar*(0.5*W_2*k4+1.0*W_4/kNucleonMass2+a6*W_5/(kNucleonMass*qv));    //Ref.[1], W_\alpha
  double T_5     = k5*W_5+m_tar*(a5/qv-v*k4)*W_2;

  double xsec    = kGF2*factor*((E_lep-k7)*(T_1+k2*T_4)/m_tar+(E_lep+k7)*T_2/(2*m_tar)
                   +n_NT*T_3*((E_nu+E_lep)*(E_lep-k7)/(2*mm_tar)-k2)-k2*T_5)
                   *(kMw2/(kMw2+Q2))*(kMw2/(kMw2+Q2))/E_nu/kPi;
  return xsec;


}
//____________________________________________________________________________
double SmithMonizQELCCPXSec::d2sQES_dQ2dv_SM(const Interaction * interaction) const
{
  Kinematics *  kinematics = interaction -> KinePtr();
  sm_utils->SetInteraction(interaction);
  const InitialState & init_state = interaction -> InitState();
  //  Assuming that the energy is greater of threshold. 
  //  See condition in method SmithMonizQELCCXSec::Integrate
  //  interaction->InitState().ProbeE(kRfLab)<sm_utils->E_nu_thr_SM()
  //  of SmithMonizQELCCXSec.cxx 
  //  if (E_nu < sm_utils->E_nu_thr_SM()) return 0;
  //  Assuming that variables Q2 and \nu are within allowable kinematic region
  //  which are specified in method: genie::utils::gsl::d2Xsec_dQ2dv::DoEval
  double Q2      = kinematics->GetKV(kKVQ2);
  double v       = kinematics->GetKV(kKVv);
  Range1D_t rkF  = sm_utils->kFQES_SM_lim(Q2,v);

  const Target & target = init_state.Tgt();
  


//  Gaussian quadratures integrate over Fermi momentum
  double R[48]= { 0.16276744849602969579e-1,0.48812985136049731112e-1,
                  0.81297495464425558994e-1,1.13695850110665920911e-1,
                  1.45973714654896941989e-1,1.78096882367618602759e-1,
                  2.10031310460567203603e-1,2.41743156163840012328e-1,
                  2.73198812591049141487e-1,3.04364944354496353024e-1,
                  3.35208522892625422616e-1,3.65696861472313635031e-1,
                  3.95797649828908603285e-1,4.25478988407300545365e-1,
                  4.54709422167743008636e-1,4.83457973920596359768e-1,
                  5.11694177154667673586e-1,5.39388108324357436227e-1,
                  5.66510418561397168404e-1,5.93032364777572080684e-1,
                  6.18925840125468570386e-1,6.44163403784967106798e-1,
                  6.68718310043916153953e-1,6.92564536642171561344e-1,
                  7.15676812348967626225e-1,7.38030643744400132851e-1,
                  7.59602341176647498703e-1,7.80369043867433217604e-1,
                  8.00308744139140817229e-1,8.19400310737931675539e-1,
                  8.37623511228187121494e-1,8.54959033434601455463e-1,
                  8.71388505909296502874e-1,8.86894517402420416057e-1,
                  9.01460635315852341319e-1,9.15071423120898074206e-1,
                  9.27712456722308690965e-1,9.39370339752755216932e-1,
                  9.50032717784437635756e-1,9.59688291448742539300e-1,
                  9.68326828463264212174e-1,9.75939174585136466453e-1,
                  9.82517263563014677447e-1,9.88054126329623799481e-1,
                  9.92543900323762624572e-1,9.95981842987209290650e-1,
                  9.98364375863181677724e-1,9.99689503883230766828e-1};

  double W[48]= { 0.00796792065552012429e-1,0.01853960788946921732e-1,
                  0.02910731817934946408e-1,0.03964554338444686674e-1,
                  0.05014202742927517693e-1,0.06058545504235961683e-1,
                  0.07096470791153865269e-1,0.08126876925698759217e-1,
                  0.09148671230783386633e-1,0.10160770535008415758e-1,
                  0.11162102099838498591e-1,0.12151604671088319635e-1,
                  0.13128229566961572637e-1,0.14090941772314860916e-1,
                  0.15038721026994938006e-1,0.15970562902562291381e-1,
                  0.16885479864245172450e-1,0.17782502316045260838e-1,
                  0.18660679627411467395e-1,0.19519081140145022410e-1,
                  0.20356797154333324595e-1,0.21172939892191298988e-1,
                  0.21966644438744349195e-1,0.22737069658329374001e-1,
                  0.23483399085926219842e-1,0.24204841792364691282e-1,
                  0.24900633222483610288e-1,0.25570036005349361499e-1,
                  0.26212340735672413913e-1,0.26826866725591762198e-1,
                  0.27412962726029242823e-1,0.27970007616848334440e-1,
                  0.28497411065085385646e-1,0.28994614150555236543e-1,
                  0.29461089958167905970e-1,0.29896344136328385984e-1,
                  0.30299915420827593794e-1,0.30671376123669149014e-1,
                  0.31010332586313837423e-1,0.31316425596861355813e-1,
                  0.31589330770727168558e-1,0.31828758894411006535e-1,
                  0.32034456231992663218e-1,0.32206204794030250669e-1,
                  0.32343822568575928429e-1,0.32447163714064269364e-1,
                  0.32516118713868835987e-1,0.32550614492363166242e-1};

  double Sum = 0;
  for(int i = 0;i<48;i++)
  {
    double kF = 0.5*(-R[i]*(rkF.max-rkF.min)+rkF.min+rkF.max);
    kinematics->SetKV(kKVPn, kF);
    Sum+=d3sQES_dQ2dvdkF_SM(interaction)*W[47-i];
    kF = 0.5*(R[i]*(rkF.max-rkF.min)+rkF.min+rkF.max);
    kinematics->SetKV(kKVPn, kF);
    Sum+=d3sQES_dQ2dvdkF_SM(interaction)*W[47-i];
  }

  double xsec = 0.5*Sum*(rkF.max-rkF.min);

  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();

  xsec *= NNucl; // nuclear xsec

  // Apply given scaling factor
  xsec *= fXSecScale;

  return xsec;

}
//____________________________________________________________________________
double SmithMonizQELCCPXSec::dsQES_dQ2_SM(const Interaction * interaction) const
{
  // Get kinematics & init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();

  double E  = init_state.ProbeE(kRfHitNucRest);
  double E2 = TMath::Power(E,2);
  double ml = interaction->FSPrimLepton()->Mass();
  double M  = target.HitNucMass();
  double q2 = kinematics.q2();

  // One of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  int sign = (is_neutrino) ? -1 : 1;

  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);

  double F1V   = fFormFactors.F1V();
  double xiF2V = fFormFactors.xiF2V();
  double FA    = fFormFactors.FA();
  double Fp    = fFormFactors.Fp();


  // Calculate auxiliary parameters
  double ml2     = TMath::Power(ml,    2);
  double M2      = TMath::Power(M,     2);
  double M4      = TMath::Power(M2,    2);
  double FA2     = TMath::Power(FA,    2);
  double Fp2     = TMath::Power(Fp,    2);
  double F1V2    = TMath::Power(F1V,   2);
  double xiF2V2  = TMath::Power(xiF2V, 2);
  double Gfactor = M2*kGF2*fVud2*(kMw2/(kMw2-q2))*(kMw2/(kMw2-q2)) / (8*kPi*E2);
  double s_u     = 4*E*M + q2 - ml2;
  double q2_M2   = q2/M2;

  // Compute free nucleon differential cross section
  double A = (0.25*(ml2-q2)/M2) * (
         (4-q2_M2)*FA2 - (4+q2_M2)*F1V2 - q2_M2*xiF2V2*(1+0.25*q2_M2)
              -4*q2_M2*F1V*xiF2V - (ml2/M2)*(
               (F1V2+xiF2V2+2*F1V*xiF2V)+(FA2+4*Fp2+4*FA*Fp)+(q2_M2-4)*Fp2));
  double B = -1 * q2_M2 * FA*(F1V+xiF2V);
  double C = 0.25*(FA2 + F1V2 - 0.25*q2_M2*xiF2V2);

  double xsec = Gfactor * (A + sign*B*s_u/M2 + C*s_u*s_u/M4);

  // Apply given scaling factor
  xsec *= fXSecScale;

  // Pauli-correction factor for deuterium, we formally apply this factor for He-3 and tritium, 
  // because RFG model is not applicable for them. 
  if (1<target.A() && target.A()<4)
  {
    double Q2 = -q2;
    double fQES_Pauli = 1.0-0.529*TMath::Exp((Q2*(228.0-531.0*Q2)-48.0)*Q2);
    xsec *= fQES_Pauli;
  }

  int nucpdgc = target.HitNucPdg();
  int NNucl = (pdg::IsProton(nucpdgc)) ? target.Z() : target.N();

  xsec *= NNucl; // nuclear xsec

  // Apply radiative correction to the cross section for IBD processes
  // Refs:
  // 1) I.S. Towner, Phys. Rev. C 58 (1998) 1288;
  // 2) J.F. Beacom, S.J. Parke, Phys. Rev. D 64 (2001) 091302;
  // 3) A. Kurylov, M.J. Ramsey-Musolf, P. Vogel, Phys. Rev. C 65 (2002) 055501;
  // 4) A. Kurylov, M.J. Ramsey-Musolf, P. Vogel, Phys. Rev. C 67 (2003) 035502.
  double rc = 1.0;
  if ( (target.IsProton() && pdg::IsAntiNuE(init_state.ProbePdg())) || (target.IsNeutron() && pdg::IsNuE(init_state.ProbePdg()) ))
  {
    const double mp  = kProtonMass;
    const double mp2 = kProtonMass2;
    const double mn2 = kNeutronMass2;
    const double Ee  = E + ( (q2 - mn2 + mp2) / 2.0 / mp );
    assert(Ee > 0.0); // must be non-zero and positive
    rc  = 6.0 + (1.5 * TMath::Log(kProtonMass / 2.0 / Ee));
    rc += 1.2 * TMath::Power((kElectronMass / Ee), 1.5);
    rc *= kAem / kPi;
    rc += 1.0;
  }

  xsec *= rc;
  return xsec;
}
