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
#include <TMath.h>
#include <Math/IFunction.h>
#include <Math/RootFinder.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/RefFrame.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Framework/Utils/Range1.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"

using namespace genie;
using namespace genie::constants;
using std::ostringstream;

//____________________________________________________________________________
SmithMonizUtils::SmithMonizUtils() :
Algorithm("genie::SmithMonizUtils")
{

}
//____________________________________________________________________________
SmithMonizUtils::SmithMonizUtils(string config) :
Algorithm("genie::SmithMonizUtils", config)
{

}
//____________________________________________________________________________
SmithMonizUtils::~SmithMonizUtils()
{

}
//____________________________________________________________________________
void SmithMonizUtils::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SmithMonizUtils::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SmithMonizUtils::LoadConfig(void)
{


  GetParam( "FermiMomentumTable", fKFTable);
  GetParam( "RFG-UseParametrization", fUseParametrization);


  // load removal energy for specific nuclei from either the algorithm's
  // configuration file or the UserPhysicsOptions file.
  // if none is used use Wapstra's semi-empirical formula.
  const std::string keyStart = "RFG-NucRemovalE@Pdg=";

  RgIMap entries = GetConfig().GetItemMap();

  for(RgIMap::const_iterator it = entries.begin(); it != entries.end(); ++it)
  {
    const std::string& key = it->first;
    int pdg = 0;
    int Z = 0;
    if (0 == key.compare(0, keyStart.size(), keyStart.c_str()))
    {
      pdg = atoi(key.c_str() + keyStart.size());
      Z = pdg::IonPdgCodeToZ(pdg);
    }
    if (0 != pdg && 0 != Z)
    {
      ostringstream key_ss ;
      key_ss << keyStart << pdg;
      RgKey rgkey   = key_ss.str();
      double eb;
      GetParam( rgkey, eb ) ;
      eb = TMath::Max(eb, 0.);
      fNucRmvE.insert(map<int,double>::value_type(Z,eb));
    }
  }


}
//____________________________________________________________________________
// Set the variables necessary for further calculations
void SmithMonizUtils::SetInteraction(const Interaction * interaction)
{

  fInteraction = interaction;
  // get kinematics & init-state parameters
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  PDGLibrary * pdglib = PDGLibrary::Instance();
  
  // neutrino energy (GeV)
  E_nu = interaction->InitState().ProbeE(kRfLab); 

  assert(target.HitNucIsSet());
  // get lepton&nuclear masses (init & final state nucleus)
  
  // mass of final charged lepton (GeV)
  m_lep = interaction->FSPrimLepton()->Mass();
  mm_lep     = TMath::Power(m_lep,    2);
  int nucl_pdg_ini = target.HitNucPdg();
  m_ini  = target.HitNucMass();
  mm_ini = TMath::Power(m_ini,    2);
  int nucl_pdg_fin = genie::pdg::SwitchProtonNeutron(nucl_pdg_ini);
  TParticlePDG * nucl_fin = pdglib->Find( nucl_pdg_fin );
  // mass of final hadron or hadron system (GeV)
  m_fin  = nucl_fin -> Mass();
  mm_fin = TMath::Power(m_fin,    2);
  // mass of target nucleus (GeV)
  m_tar = target.Mass();
  mm_tar = TMath::Power(m_tar,    2);

  // RFG is not applied for A<4
  if (target.A()<4)
  {
    E_BIN = P_Fermi = m_rnu = mm_rnu = 0;
    return;
  }

  bool is_p = pdg::IsProton(nucl_pdg_ini);
  int Zi = target.Z();
  int Ai = target.A();
  int Zf = (is_p) ? Zi-1 : Zi;
  int Af = Ai-1;
  TParticlePDG * nucl_f = pdglib->Find( pdg::IonPdgCode(Af, Zf) );
  if(!nucl_f)
  {
    LOG("SmithMoniz", pFATAL)
    << "Unknwown nuclear target! No target with code: "
    << pdg::IonPdgCode(Af, Zf) << " in PDGLibrary!";
    exit(1);
  }
  // mass of residual nucleus (GeV)
  m_rnu = nucl_f -> Mass();
  mm_rnu = TMath::Power(m_rnu, 2);

  int Z = target.Z();
  int A = target.A();
  int N = A-Z;


  // maximum value of Fermi momentum of target nucleon (GeV)
  if (A < 6 || !fUseParametrization)
  {
     // look up the Fermi momentum for this Target
     FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
     const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
     P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucl_pdg_ini);
  }
  else
  {
    // define the Fermi momentum for this Target
    //
    P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
    // correct the Fermi momentum for the struck nucleon
    if(is_p) P_Fermi *= TMath::Power( 2.*Z/A, 1./3);
    else
           P_Fermi *= TMath::Power( 2.*N/A, 1./3);
  }

  // neutrino binding energy (GeV)
  if (target.A() < 6 || !fUseParametrization)
  {
     map<int,double>::const_iterator it = fNucRmvE.find(Z);
     if(it != fNucRmvE.end()) E_BIN = it->second;
       else E_BIN = utils::nuclear::BindEnergyPerNucleon(target);
  }
  else
    E_BIN = utils::nuclear::BindEnergyPerNucleonParametrization(target);




}
//____________________________________________________________________________
// Get the neutrino energy threshold
double SmithMonizUtils::E_nu_thr_SM(void) const
{

  Func1D<SmithMonizUtils> QEL_EnuMin_SM_(*this, &SmithMonizUtils::QEL_EnuMin_SM);

  // maximum of function call
  const int MFC = 10000;
  const double EPSABS = 0.;
  const double EPSREL = 1.0e-08;
  
  // Energy threshold of scattering on nucleus (Eq. (1) of Ref. 3)
  double E_min = ((m_lep + m_rnu + m_fin)*(m_lep + m_rnu + m_fin) - mm_tar)/2/m_tar;
  
  // Energy threshold of scattering on bound nucleon (Eq. (2) of Ref. 3)
  double E_min2 = ((m_lep + m_fin)*(m_lep + m_fin)-mm_ini-E_BIN*E_BIN+2*E_BIN*TMath::Sqrt(mm_ini+P_Fermi*P_Fermi))/2/(TMath::Sqrt(mm_ini+P_Fermi*P_Fermi)-E_BIN+P_Fermi);
  
  E_min = TMath::Max(E_min, E_min2);
  
  // if nu_1>nu_max at E_min then we try to find energy in range (E_min, 50.) where nu_1=nu_max and put E_min equal to it.
  // nu_1   -- minimal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3), 
  // nu_max -- maximal energy transfer on nucleus (see Eq. (9) of Ref.3)
  if (QEL_EnuMin_SM(E_min) > 0)
  {
    // C++ analog of fortran function Enu_rf= DZEROX(E_min,Enu_2,EPS,MFC,QEL_EnuMin_SM,1)
    ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
    //convergence is reached using tolerance = 2 *( epsrel * abs(x) + epsabs)
    if ( rfgb.Solve(QEL_EnuMin_SM_, E_min, 50., MFC, EPSABS, EPSREL) )
    {
       E_min = rfgb.Root();
    }
  }
  E_min = TMath::Max(E_min, 0.);
  return E_min;

}
//____________________________________________________________________________
// The auxiliary function for determining energy threshold
double SmithMonizUtils::QEL_EnuMin_SM(double Enu) const
{
  // return the minimum of nu_1-nu_max as function of Q2 in range ( Q2_lim1(Enu), Q2_lim2(Enu) )
  const double EPS  = 1.0e-06;
  const double Delta= 1.0e-14;
  const double Precision = std::numeric_limits<double>::epsilon();
  double s = 2*Enu*m_tar+mm_tar;
  double W2 = (m_rnu+m_fin)*(m_rnu+m_fin);
  // neutrino energy in CMS (see Eq. (4) of Ref.3)
  double E_nu_CM = (s-mm_tar)/2/TMath::Sqrt(s);
  // final lepton energy and momentum at W2_min (see Eqs. (5) and (6) of Ref.3)
  double E_l_CM  = (s+mm_lep-W2)/2/TMath::Sqrt(s);
  double P_l_CM  = E_l_CM>m_lep?TMath::Sqrt(E_l_CM*E_l_CM-mm_lep):Precision;
  // minimal and maximal allowed Q2 for scattering on nucleus (see Eqs. (3) of Ref.3)
  double Q2_lim1 = 2*E_nu_CM*(E_l_CM-P_l_CM)-mm_lep;
  double Q2_lim2 = 2*E_nu_CM*(E_l_CM+P_l_CM)-mm_lep;
  // C++ analog of fortran function DMINFC(Q2lim1_SM,Q2_lim1,Q2_lim2,EPS,Delta,Q2_0,F_MIN,LLM)
  Func1Dpar<SmithMonizUtils> Q2lim1_SM_(*this, &SmithMonizUtils::Q2lim1_SM, Enu);
  double Q2_0,F_MIN;
  bool LLM;
  // find minimum of nu_1-nu_max as function of Q2 in range (Q2_lim1,Q2_lim2)
  DMINFC(Q2lim1_SM_,Q2_lim1,Q2_lim2,EPS,Delta,Q2_0,F_MIN,LLM);
  return F_MIN;
}
//____________________________________________________________________________
// The auxiliary function for determining Q2-range
double SmithMonizUtils::Q2lim1_SM(double Q2, double Enu) const
{
  // maximal energy transfer (see Eq. (9) of Ref.3) 
  double nu_max = Enu*Q2/(Q2+mm_lep)-(Q2+mm_lep)/4/Enu;
  
  double E = sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = 0.5*(Q2+mm_fin-b);
  // minimal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3), 
  double nu_1 = (a*(E-E_BIN)-P_Fermi*TMath::Sqrt(a*a+Q2*b))/b;
  return nu_1-nu_max;

}
//____________________________________________________________________________
// The auxiliary function for determining Q2-range
double SmithMonizUtils::Q2lim2_SM(double Q2) const
{
  // minimal energy transfer for scattering on nucleus (see Eq. (7) of Ref.3)
  double nu_min = ((m_rnu+m_fin)*(m_rnu+m_fin)+Q2-mm_tar)/(2*m_tar);
  
  double E = sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = (Q2+mm_fin-b)*0.5;
  // maximal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3)
  double nu_2  = (a*(E-E_BIN)+P_Fermi*TMath::Sqrt(a*a+Q2*b))/b;
  return nu_min-nu_2;

}
//____________________________________________________________________________
// Return allowed Q2-range
Range1D_t SmithMonizUtils::Q2QES_SM_lim(void) const
{

  // here we find Q2_min and Q2_max such that
  // 0. Q2_min>=0
  // 1. nu_1(Q2_min)<=nu_max(Q2_min)
  // 2. nu_1(Q2_max)<=nu_max(Q2_max)
  // 3. nu_min(Q2_min)<=nu_2(Q2_min)
  // 4. Q2_min>=the value of minimal Q2 defined for scattering on nucleus 
  // 5. Q2_max<=the value of maximal Q2 defined for scattering on nucleus (see Eqs. (3) of Ref.3)
  // nu_1 --   minimal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3), 
  // nu_max -- maximal energy transfer on nucleus (see Eq. (9) of Ref.3),
  // nu_2 --   maximal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3), 
  // nu_min -- minimal energy transfer on nucleus (see Eq. (7) of Ref.3)
  
  // maximum of function call
  const int MFC = 1000;
  const double EPS  = 1.0e-08;
  const double Delta= 1.0e-14;
  const double EPSABS = 0;
  const double EPSREL = 1.0e-08;
  const double Precision = std::numeric_limits<double>::epsilon();
  // if the nucleus mass is less than 4 then this is a special case
  if (E_BIN == 0 && P_Fermi == 0)
  {
    double s = 2*E_nu*m_ini+mm_ini;
    // minimal W2 for scattering on nucleus (see Eq. (6) of Ref.3)
    double W2 = mm_fin;
    // neutrino energy in CMS (see Eq. (4) of Ref.3)
    double E_nu_CM = (s-mm_ini)/2/TMath::Sqrt(s);
    // final lepton energy and momentum at W2_min (see Eqs. (5) of Ref.3)
    double E_l_CM  = (s+mm_lep-W2)/2/TMath::Sqrt(s);
    double P_l_CM  = E_l_CM>m_lep?TMath::Sqrt(E_l_CM*E_l_CM-mm_lep):Precision;
    // minimal and maximal allowed Q2 for scattering on nucleus (see Eqs. (3) of Ref.3)
    double Q2_min = 2*E_nu_CM*(E_l_CM-P_l_CM)-mm_lep;
    double Q2_max = 2*E_nu_CM*(E_l_CM+P_l_CM)-mm_lep;
    Q2_min= TMath::Max(Q2_min,0.0);
    Range1D_t R(Q2_min,Q2_max);
    return R;
  }
  double s = 2*E_nu*m_tar+mm_tar;
  // minimal W2 for scattering on nucleus (see Eq. (6) of Ref.3)
  double W2 = (m_rnu+m_fin)*(m_rnu+m_fin);
  // neutrino energy in CMS (see Eq. (4) of Ref.3)
  double E_nu_CM = (s-mm_tar)/2/TMath::Sqrt(s);
  // final lepton energy and momentum at W2_min (see Eqs. (5) of Ref.3)
  double E_l_CM  = (s+mm_lep-W2)/2/TMath::Sqrt(s);
  double P_l_CM  = E_l_CM>m_lep?TMath::Sqrt(E_l_CM*E_l_CM-mm_lep):Precision;
  // minimal and maximal allowed Q2 for scattering on nucleus (see Eqs. (3) of Ref.3)
  double Q2_min = 2*E_nu_CM*(E_l_CM-P_l_CM)-mm_lep;
  double Q2_max = 2*E_nu_CM*(E_l_CM+P_l_CM)-mm_lep;
  double F_MIN, Q2_0;
  bool LLM;
  // C++ analog of fortran function DMINFC(Q2lim1_SM,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM)
  Func1Dpar<SmithMonizUtils> Q2lim1_SM_(*this, &SmithMonizUtils::Q2lim1_SM, E_nu);
  // if minimum of nu_1-nu_max>0 then exit with error, because it's impossible
  DMINFC(Q2lim1_SM_,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM);
  if (F_MIN>0)
  {
      LOG("SmithMoniz", pFATAL)
      << "No overlapped area for energy " << E_nu << "\n" <<
      "Q2_min=" << Q2_min << " Q2_max=" << Q2_max << "\n" <<
      "Q2_0=" << Q2_0 << " F_MIN=" << F_MIN;
      exit(1);
  }
  // at Q2_0 here we have: nu_1(Q2_0)-nu_max(Q2_0)<0 
  // if nu_1(Q2_min)-nu_max(Q2_min)>0 we find corrected Q2_min_cor>Q2_min where nu_1(Q2_min_cor)-nu_max(Q2_min_cor)=0
  // (it is always possible because of above conditions) 
  if (Q2lim1_SM(Q2_min, E_nu)>0)
  {
    //C++ analog of fortran function Q2_RF=DZEROX(Q2_min,Q2_0,EPS,MFC,Q2lim1_SM,1)
    ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
    // convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
    if (rfgb.Solve(Q2lim1_SM_, Q2_min, Q2_0, MFC, EPSABS, EPSREL))
    {
      Q2_min= rfgb.Root();
    }
  }
  // if nu_1(Q2_max)-nu_max(Q2_max)>0 we find Q2_max_cor<Q2_max where nu_1(Q2_max_cor)-nu_max(Q2_max_cor)=0
  // (it is always possible because of above conditions) 
  if(Q2lim1_SM(Q2_max, E_nu)>0)
  {
     // C++ analog of fortran function Q2_RF=DZEROX(Q2_0,Q2_max,Eps,MFC,Q2lim1_SM,1)
     ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
     //convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
     if (rfgb.Solve(Q2lim1_SM_, Q2_0, Q2_max, MFC, EPSABS, EPSREL))
     {
       Q2_max= rfgb.Root();
     }
  }
  Func1D<SmithMonizUtils> Q2lim2_SM_(*this, &SmithMonizUtils::Q2lim2_SM);
  // if nu_min(Q2_min)-nu_2(Q2_min)>0 and nu_min(Q2_max)-nu_2(Q2_max)>0 then set Q2_min=Q2_max (it makes xsec equal to zero).
  if (Q2lim2_SM(Q2_min)>0)
  {
     if(Q2lim2_SM(Q2_max)>0)
     {
           LOG("SmithMoniz", pWARN) << "The RFG model is not applicable! The cross section is set zero!";
           Q2_min = Q2_max;
     }
     // here we have nu_min(Q2_min)-nu_2(Q2_min)>0 and nu_min(Q2_max)-nu_2(Q2_max)<0 or vice versa
     // so we always can find Q2_min_cor where nu_min(Q2_min_cor)-nu_2(Q2_min_cor)=0
     else
     {
       // C++ analog of fortran function Q2_RF = DZEROX(Q2_min,Q2_max,Eps,MFC,Q2lim2_SM,1)
       ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
       // convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
       if (rfgb.Solve(Q2lim2_SM_, Q2_min,Q2_max, MFC, EPSABS, EPSREL))
       {
          Q2_min= rfgb.Root();
       }
     }
  }
  Q2_min = TMath::Max(Q2_min,0.0);

  Range1D_t R(Q2_min,Q2_max);
  return R;

}
//____________________________________________________________________________
// Return allowed v-range for given Q2
Range1D_t SmithMonizUtils::vQES_SM_lim(double Q2) const
{

  // minimal energy transfer for scattering on nucleus (see Eq. (7) of Ref.3)
  double nu_min= ((m_rnu+m_fin)*(m_rnu+m_fin)+Q2-mm_tar)/2/m_tar;
  
  // if the target is nucleon then nu_min=nu_max=(m_fin^2+Q^2-m_ini^2)/(2*m_ini)
  if (E_BIN == 0 && P_Fermi == 0)
    return Range1D_t(nu_min, nu_min);
  
  // maximal energy transfer (see Eq. (9) of Ref.3)
  double nu_max = E_nu*Q2/(Q2+mm_lep)-(Q2+mm_lep)/4/E_nu;
  
  // now we find limits for bound nucleon
  double E = TMath::Sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = (Q2+mm_fin-b)*0.5;
  double tmp1 = a*(E-E_BIN);
  double tmp2  = P_Fermi*TMath::Sqrt(a*a+Q2*b);
  // minimal and maximal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3)
  double nu_1 = (tmp1-tmp2)/b;
  double nu_2 = (tmp1+tmp2)/b;
  // for minimal energy transfer we take maximum of corresponding values on nucleus and bound nucleon
  nu_min= TMath::Max(nu_min,nu_1);
  // for maximal energy transfer we take minimum of corresponding values on nucleus and bound nucleon
  nu_max= TMath::Min(nu_max,nu_2);
  
  if (nu_min<=nu_max)
    return Range1D_t(nu_min,nu_max);
  else 
  // to avoid machine precision errors
    return Range1D_t(0.5*(nu_min+nu_max),0.5*(nu_min+nu_max));

}
//____________________________________________________________________________
// Return minimum of low edge of v-range for given neutrino energy
double SmithMonizUtils::vQES_SM_min(double Q2_min, double Q2_max) const
{
  // maximum of function call
  const double EPS  = 1.0e-08;
  const double Delta= 1.0e-14;
   
  if (E_BIN == 0 && P_Fermi == 0)
    return vQES_SM_lim(Q2_min).min;
  
  double F_MIN, Q2_0;
  bool LLM;
  // C++ analog of fortran function DMINFC(Q2lim1_SM,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM)
  Func1D<SmithMonizUtils> vlim1_SM_(*this, &SmithMonizUtils::vlim1_SM);
  DMINFC(vlim1_SM_,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM);
  
  return F_MIN;
  
}
//____________________________________________________________________________
// Return maximum of low edge of v-range for given neutrino energy
double SmithMonizUtils::vQES_SM_max(double Q2_min, double Q2_max) const
{
  // maximum of function call
  const double EPS  = 1.0e-08;
  const double Delta= 1.0e-14;
     
  if (E_BIN == 0 && P_Fermi == 0)
    return vQES_SM_lim(Q2_max).min;
  
  double F_MIN, Q2_0;
  bool LLM;
  // C++ analog of fortran function DMINFC(Q2lim1_SM,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM)
  Func1D<SmithMonizUtils> vlim2_SM_(*this, &SmithMonizUtils::vlim2_SM);
  DMINFC(vlim2_SM_,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM);
  
  return -F_MIN;
  
}
//____________________________________________________________________________
// The auxiliary function for determining minimum of low edge of v-range
double SmithMonizUtils::vlim1_SM(double Q2) const
{
  // minimal energy transfer for scattering on nucleus (see Eq. (7) of Ref.3)
  double nu_min = ((m_rnu+m_fin)*(m_rnu+m_fin)+Q2-mm_tar)/(2*m_tar);
  
  double E = sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = (Q2+mm_fin-b)*0.5;
  // minimal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3)
  double nu_1  = (a*(E-E_BIN)-P_Fermi*TMath::Sqrt(a*a+Q2*b))/b;
  nu_min= TMath::Max(nu_min,nu_1);
  return nu_min;
}
//____________________________________________________________________________
// The auxiliary function for determining maximum of up edge of v-range
double SmithMonizUtils::vlim2_SM(double Q2) const
{
  // maximal energy transfer (see Eq. (9) of Ref.3) 
  double nu_max = E_nu*Q2/(Q2+mm_lep)-(Q2+mm_lep)/4/E_nu;
  
  double E = sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = 0.5*(Q2+mm_fin-b);
  // maximal energy transfer for bound nucleon (see Eqs. (11) of Ref. 3)
  double nu_2 = (a*(E-E_BIN)+P_Fermi*TMath::Sqrt(a*a+Q2*b))/b;
  nu_max= TMath::Min(nu_max,nu_2);
  return -nu_max;
}
//____________________________________________________________________________
// Return allowed Fermi momentum range for given Q2 and v
Range1D_t SmithMonizUtils::kFQES_SM_lim(double Q2, double nu) const
{
  double qv = TMath::Sqrt(nu*nu+Q2);
  double c_f = (nu-E_BIN)/qv;
  double d_f = (E_BIN*E_BIN-2*nu*E_BIN-Q2+mm_ini-mm_fin)/(2*qv*m_ini);
  // minimal energy of initial nucleon (see Eq. (13) of Ref.3)
  double Ef_min= TMath::Max(m_ini*(c_f*d_f+TMath::Sqrt(1.0-c_f*c_f+d_f*d_f))/(1.0-c_f*c_f), TMath::Sqrt(P_Fermi*P_Fermi+mm_ini)-nu);
  double kF_min= P_Fermi!=0?TMath::Sqrt(TMath::Max(Ef_min*Ef_min-mm_ini,0.0)):0.;
  double kF_max= P_Fermi;
  Range1D_t R;
  if (kF_min<=kF_max)
    R = Range1D_t(kF_min,kF_max);
  else
    R = Range1D_t(0.5*(kF_min+kF_max),0.5*(kF_min+kF_max));
  return R;

}
//____________________________________________________________________________
// C++ implementation of DMINC function from CERN library
void SmithMonizUtils::DMINFC(Functor1D &F, double A,double B, double EPS, double DELTA, double &X, double &Y, bool &LLM) const
{
  const double W5 = TMath::Sqrt(5);
  const double HV = (3-W5)/2;
  const double HW = (W5-1)/2;
  const double R1 = 1.0;
  const double HF = R1/2;

  int N = -1;
  if(A!=B) N = TMath::Nint(2.08*TMath::Log(TMath::Abs((A-B)/EPS)));
  double C = A;
  double D = B;
  if(A>B)
  {
     C = B;
     D = A;
  }
  bool LLT = true;
  bool LGE = true;
  double V, FV, W, FW, H;
  while(true)
  {
     H = D-C;
     if(N<0)
     {
       X = HF*(C+D);
       Y = F(X);
       LLM = TMath::Abs(X-A)>DELTA && TMath::Abs(X-B)>DELTA;
       return;
     }
     if(LLT)
     {
       V = C+HV*H;
       FV = F(V);
     }
     if(LGE)
     {
        W = C+HW*H;
        FW = F(W);
     }
     if(FV<FW)
     {
        LLT = true;
        LGE = false;
        D = W;
        W = V;
        FW = FV;
     }
     else
     {
        LLT = false;
        LGE = true;
        C = V;
        V = W;
        FV = FW;
     }
     N = N-1;
  }
}
//____________________________________________________________________________
// Density of Fermi gas, for case T_Fermi=0 is a step functio, which is blurred at T_Fermi!=0
double SmithMonizUtils::rho(double P_Fermi, double T_Fermi, double p)
{

  if (T_Fermi==0)
  {
    //Pure Fermi gaz with T_Fermi=0
    if(p<=P_Fermi)
      return 1.0;
    else
      return 0.0;
  }
  else
  {
      //Fermi-Dirac distribution
      return 1.0/(1.0 + TMath::Exp(-(P_Fermi-p)/T_Fermi));
  }


}
//____________________________________________________________________________
double SmithMonizUtils::GetBindingEnergy(void) const
{
  return E_BIN;
}
//____________________________________________________________________________
double SmithMonizUtils::GetFermiMomentum(void) const
{
  return P_Fermi;
}
//____________________________________________________________________________
double SmithMonizUtils::GetTheta_k(double v, double qv) const
{
  return TMath::ACos((v + (mm_lep-v*v+qv*qv)/2/E_nu)/qv);
}
//____________________________________________________________________________
double SmithMonizUtils::GetTheta_p(double pv, double v, double qv, double &E_p) const
{
  E_p = TMath::Sqrt(mm_ini+pv*pv)-E_BIN;
  if (pv != 0)
    return TMath::ACos(((v-E_BIN)*(2*E_p+v+E_BIN)-qv*qv+mm_ini-mm_fin)/(2*pv*qv));
  else
    return 0;
}
//____________________________________________________________________________
double SmithMonizUtils::PhaseSpaceVolume(KinePhaseSpace_t ps) const
{
   double vol = 0.0;
   if (ps == kPSQ2fE)
   {
     Range1D_t rQ2 = Q2QES_SM_lim();
     vol = rQ2.max - rQ2.min;
     return vol;
   }
   else if (ps == kPSQ2vpfE)
   {
     const int kNQ2 = 101;
     const int kNv  = 101;
     Range1D_t rQ2 = Q2QES_SM_lim();
     double dQ2 = (rQ2.max-rQ2.min)/(kNQ2-1);
     for(int iq2=0; iq2<kNQ2-1; iq2++)
     {
       double Q2 = rQ2.min + iq2*dQ2;
       Range1D_t rv  = vQES_SM_lim(Q2);
       double dv = (rv.max-rv.min)/(kNv-1);
       for(int iv=0; iv<kNv-1; iv++)
       {
         double v = rv.min + iv*dv;
         Range1D_t rkF  = kFQES_SM_lim(Q2,v);
         double dkF = (rkF.max-rkF.min);
         vol += (dQ2*dv*dkF);
       }
     }
     return vol;
   }
   else
     return 0;
}
