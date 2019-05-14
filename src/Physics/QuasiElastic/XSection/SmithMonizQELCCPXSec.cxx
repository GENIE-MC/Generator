//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          adapted from  fortran code provided by 
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>, Joint Institute for Nuclear Research
          Vladimir Lyubushkin, Joint Institute for Nuclear Research
          Vadim Naumov <vnaumov@theor.jinr.ru>, Joint Institute for Nuclear Research
          based on code of Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab
           
 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________
#include <sstream>

#include <TMath.h>


#include "Framework/Algorithm/AlgFactory.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/QuasiElastic/XSection/QELFormFactors.h"
#include "Physics/QuasiElastic/XSection/QELFormFactorsModelI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Conventions/KineVar.h"
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
  
  if(! this -> ValidProcess (interaction) ) 
  {
    LOG("SmithMoniz",pWARN) << "not a valid process"; 
    return 0.;
  }
  
  if(kps == kPSQ2fE) 
  {
     if(! this -> ValidKinematics (interaction) )
     {
        LOG("SmithMoniz",pWARN) << "not valid kinematics"; 
        return 0.;
     }
     return this->dsQES_dQ2_SM(interaction);
  }
  
  if(kps == kPSQ2vfE) {
    return this->d2sQES_dQ2dv_SM(interaction);
  }
  
  return 0;
  
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
  // Get kinematics & init-state parameters
  const Kinematics &   kinematics = interaction -> Kine();
  
  double kF      = kinematics.GetKV(kKVPn);
  double kkF     = kF*kF;
  double P_Fermi, E_nuBIN;
  
  E_nuBIN = sm_utils->GetBindingEnergy();
  
  double E_p     = TMath::Sqrt(fmm_ini+kkF)-E_nuBIN;
  double cosT_p  = ((fv-E_nuBIN)*(2*E_p+fv+E_nuBIN)-fqqv+fmm_ini-fmm_fin)/(2*kF*fqv);           //\cos\theta_p
  if (cosT_p < -1.0 || cosT_p > 1.0 ) return 0.0;
  double pF      = TMath::Sqrt(kkF+(2*kF*fqv)*cosT_p+fqqv);
  double b2_flux = (E_p-kF*fcosT_k*cosT_p)*(E_p-kF*fcosT_k*cosT_p);
  double c2_flux = kkF*(1-cosT_p*cosT_p)*(1-fcosT_k*fcosT_k);
  
  P_Fermi        = sm_utils->GetFermiMomentum();
  double FV_SM   = 4.0*TMath::Pi()/3*TMath::Power(P_Fermi, 3);
  double factor  = fk1*(fm_tar*kF/(FV_SM*fqv*TMath::Sqrt(b2_flux-c2_flux)))*SmithMonizUtils::rho(P_Fermi, 0.0, kF)*(1-SmithMonizUtils::rho(P_Fermi, 0.01, pF));    
  
  double a2      = kkF/kNucleonMass2;
  double a3      = a2*cosT_p*cosT_p;
  double a6      = kF*cosT_p/kNucleonMass;
  double a7      = E_p/kNucleonMass;
  double a4      = a7*a7;
  double a5      = 2*a7*a6;
  
  double k3      = fv/fqv;
  double k4      = (3*a3-a2)/fqqv;
  double k5      = (a7-a6*k3)*fm_tar/kNucleonMass;
  
  double T_1     = 1.0*fW_1+(a2-a3)*0.5*fW_2;                              //Ref.[1], W_1
  double T_2     = ((a2-a3)*fQ2/(2*fqqv)+a4-k3*(a5-k3*a3))*fW_2;           //Ref.[1], W_2
  double T_3     = k5*fW_3;                                                //Ref.[1], W_8
  double T_4     = fmm_tar*(0.5*fW_2*k4+1.0*fW_4/kNucleonMass2+a6*fW_5/(kNucleonMass*fqv));    //Ref.[1], W_\alpha
  double T_5     = k5*fW_5+fm_tar*(a5/fqv-fv*k4)*fW_2;
  
  double xsec    = kGF2*factor*((fE_lep-fk7)*(T_1+fk2*T_4)/fm_tar+(fE_lep+fk7)*T_2/(2*fm_tar)
                   +fn_NT*T_3*((fE_nu+fE_lep)*(fE_lep-fk7)/(2*fmm_tar)-fk2)-fk2*T_5)
                   *(kMw2/(kMw2+fQ2))*(kMw2/(kMw2+fQ2))/fE_nu/kPi;
  return xsec;
  
  
}
//____________________________________________________________________________
double SmithMonizQELCCPXSec::d2sQES_dQ2dv_SM(const Interaction * interaction) const
{
  Kinematics *  kinematics = interaction -> KinePtr();
  sm_utils->SetInteraction(interaction);
  fQ2      = kinematics->GetKV(kKVQ2);
  fv       = kinematics->GetKV(kKVv);
  Range1D_t rkF = sm_utils->kFQES_SM_lim(fQ2,fv);
  
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  PDGLibrary * pdglib = PDGLibrary::Instance();
  
  // One of the xsec terms changes sign for antineutrinos
  bool is_neutrino = pdg::IsNeutrino(init_state.ProbePdg());
  fn_NT = (is_neutrino) ? +1 : -1;
  
  int nucl_pdg_ini = target.HitNucPdg();
  double m_ini  = target.HitNucMass();                         
  fmm_ini = TMath::Power(m_ini,    2);
  int nucl_pdg_fin = genie::pdg::SwitchProtonNeutron(nucl_pdg_ini);
  TParticlePDG * nucl_fin = pdglib->Find( nucl_pdg_fin );
  double m_fin  = nucl_fin -> Mass();                         //  Mass of final hadron or hadron system (GeV)
  fmm_fin = TMath::Power(m_fin,    2);
  fm_tar = target.Mass();                                     //  Mass of target nucleus (GeV)
  fmm_tar = TMath::Power(fm_tar,    2);
  
  fE_nu  = init_state.ProbeE(kRfLab);
  fE_lep   = fE_nu-fv;
  double m_lep = interaction->FSPrimLepton()->Mass();
  double mm_lep = m_lep*m_lep;
  if (fE_lep < m_lep) return 0.0;
  double P_lep   = TMath::Sqrt(fE_lep*fE_lep-mm_lep);
  double k6 = (fQ2+mm_lep)/(2*fE_nu);
  double cosT_lep= (fE_lep-k6)/P_lep;
  if (cosT_lep < -1.0 || cosT_lep > 1.0 ) return 0.0;
  //|\vec{q}|
  fqqv     = fv*fv+fQ2;
  fqv      = TMath::Sqrt(fqqv);
  fcosT_k  = (fv+k6)/fqv;
  if (fcosT_k < -1.0 || fcosT_k > 1.0 ) return 0.0;
  
  fk1 = fVud2*kNucleonMass2*kPi;
  fk2 = mm_lep/(2*fmm_tar);
  fk7 = P_lep*cosT_lep;
  
  
  // Calculate the QEL form factors
  fFormFactors.Calculate(interaction);
  fF_V   = fFormFactors.F1V();
  fF_M   = fFormFactors.xiF2V();
  fF_A   = fFormFactors.FA();
  fF_P   = fFormFactors.Fp();
  fFF_V  = fF_V*fF_V;
  fFF_M  = fF_M*fF_M;
  fFF_A  = fF_A*fF_A;
  
  double t = fQ2/(4*kNucleonMass2);
  fW_1     = fFF_A*(1+t)+t*(fF_V+fF_M)*(fF_V+fF_M);                   //Ref.[1], \tilde{T}_1
  fW_2     = fFF_A+fFF_V+t*fFF_M;                                     //Ref.[1], \tilde{T}_2
  fW_3     =-2*fF_A*(fF_V+fF_M);                                      //Ref.[1], \tilde{T}_8
  fW_4     =-0.5*fF_V*fF_M-fF_A*fF_P+t*fF_P*fF_P-0.25*(1-t)*fFF_M;    //Ref.[1], \tilde{T}_\alpha
  fW_5     = fFF_V+t*fFF_M+fFF_A;
  
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
  
  // Deuterium and tritium is a special case
  if (target.A()>1 && target.A()<4)
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

