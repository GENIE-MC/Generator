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
  
  
  // Load removal energy for specific nuclei from either the algorithm's
  // configuration file or the UserPhysicsOptions file.
  // If none is used use Wapstra's semi-empirical formula.
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
  // Get kinematics & init-state parameters
  // unused // const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  const Target & target = init_state.Tgt();
  PDGLibrary * pdglib = PDGLibrary::Instance();
  
  E_nu = interaction->InitState().ProbeE(kRfLab);         //  Neutrino energy (GeV)
  
  assert(target.HitNucIsSet());
  // get lepton&nuclear masses (init & final state nucleus)
  m_lep = interaction->FSPrimLepton()->Mass();          //  Mass of final charged lepton (GeV)
  mm_lep     = TMath::Power(m_lep,    2);
  int nucl_pdg_ini = target.HitNucPdg();
  m_ini  = target.HitNucMass();
  mm_ini = TMath::Power(m_ini,    2);
  int nucl_pdg_fin = genie::pdg::SwitchProtonNeutron(nucl_pdg_ini);
  TParticlePDG * nucl_fin = pdglib->Find( nucl_pdg_fin );
  m_fin  = nucl_fin -> Mass();                         //  Mass of final hadron or hadron system (GeV)
  mm_fin = TMath::Power(m_fin,    2);
  m_tar = target.Mass();                             //  Mass of target nucleus (GeV)
  mm_tar = TMath::Power(m_tar,    2);

  // For hydrogen and deuterium RFG is not applied
  if (target.A()<3)
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

  m_rnu = nucl_f -> Mass();                          //  Mass of residual nucleus (GeV)
  mm_rnu = TMath::Power(m_rnu, 2);
  
  int Z = target.Z();
  int A = target.A();
  int N = A-Z;
  
  
  // Maximum value of Fermi momentum of target nucleon (GeV)
  if (A < 6 || !fUseParametrization)
  {
     //-- look up the Fermi momentum for this Target
     FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
     const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
     P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucl_pdg_ini);
  }
  else
  {
    //-- define the Fermi momentum for this Target
    //
    P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
    //-- correct the Fermi momentum for the struck nucleon
    if(is_p) P_Fermi *= TMath::Power( 2.*Z/A, 1./3);
    else
           P_Fermi *= TMath::Power( 2.*N/A, 1./3);
  }
  
  // Neutrino binding energy (GeV)
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
  
  
  const int MFC = 10000;      //  Maximum of function call
  const double EPSABS = 0;
  const double EPSREL = 1.0e-08;
  double E_min = ((m_lep + m_rnu + m_fin)*(m_lep + m_rnu + m_fin) - mm_tar)/(2*m_tar);
  double Enu_2 = 5.0e+00;
  double Enu_rf;
  if (QEL_EnuMin_SM(E_min) > 0)
  {
    // C++ analog of fortran function Enu_rf= DZEROX(E_min,Enu_2,EPS,MFC,QEL_EnuMin_SM,1)
    ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
    rfgb.Solve(QEL_EnuMin_SM_, E_min, Enu_2, MFC, EPSABS, EPSREL); //convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
    Enu_rf = rfgb.Root();
  }
  else
  {
    Enu_rf = -1.0e+01;
  }
  E_min = TMath::Max(E_min,Enu_rf);
  if(E_min < 0)
  {
    E_min = 0;
    LOG("SmithMoniz", pDEBUG) << "E_min = " << E_min;
  }
  return E_min;

}
//____________________________________________________________________________
// The auxiliary function for determining energy threshold
double SmithMonizUtils::QEL_EnuMin_SM(double Enu) const
{
  const double EPS  = 1.0e-06;
  const double Delta= 1.0e-14;
  const double Precision = std::numeric_limits<double>::epsilon();
  Enu_in = Enu;
  double s = 2*Enu*m_tar+mm_tar;
  double W2 = (m_rnu+m_fin)*(m_rnu+m_fin);
  double c  = 0.5*(W2+mm_lep-mm_tar*(W2-mm_lep)/s);
  double sqrtD  = TMath::Sqrt(TMath::Max(Precision,LambdaFUNCTION(1.0, mm_lep/s, W2/s)));
  double tmp    = 0.5*(s-mm_tar);
  double Q2_lim1= tmp*(1.0-sqrtD)-c;
  double Q2_lim2= tmp*(1.0+sqrtD)-c;
  // C++ analog of fortran function DMINFC(Q2lim1_SM,Q2_lim1,Q2_lim2,EPS,Delta,Q2_0,F_MIN,LLM)
  Func1D<SmithMonizUtils> Q2lim1_SM_(*this, &SmithMonizUtils::Q2lim1_SM);
  double Q2_0,F_MIN;
  bool LLM;
  DMINFC(Q2lim1_SM_,Q2_lim1,Q2_lim2,EPS,Delta,Q2_0,F_MIN,LLM);
  return F_MIN;
}
//____________________________________________________________________________
// The auxiliary function for determining Q2-range
double SmithMonizUtils::Q2lim1_SM(double Q2) const
{
  double s = 2*Enu_in*m_tar+mm_tar;
  double nu_max = (s-mm_tar-mm_lep*(s-mm_tar)/(Q2+mm_lep)-mm_tar*(Q2+mm_lep)/(s-mm_tar))/(2*m_tar);
  double E = sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = 0.5*(Q2+mm_fin-b);
  double nu_1 = (a*(E-E_BIN)-P_Fermi*TMath::Sqrt(a*a+Q2*b))/b;
  return nu_1-nu_max;

}
//____________________________________________________________________________
// The auxiliary function for determining Q2-range
double SmithMonizUtils::Q2lim2_SM(double Q2) const
{
  double nu_min = ((m_rnu+m_fin)*(m_rnu+m_fin)+Q2-mm_tar)/(2*m_tar);
  double E = sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = (Q2+mm_fin-b)*0.5;
  double nu_2  = (a*(E-E_BIN)+P_Fermi*TMath::Sqrt(a*a+Q2*b))/b;
  return nu_min-nu_2;

}
//____________________________________________________________________________
// Return allowed Q2-range
Range1D_t SmithMonizUtils::Q2QES_SM_lim(void) const
{


  const int MFC = 1000;      //  Maximum of function call
  const double EPS  = 1.0e-08;
  const double Delta= 1.0e-14;
  const double EPSABS = 0;
  const double EPSREL = 1.0e-08;

  Enu_in = E_nu;
  double s = 2*E_nu*m_tar+mm_tar;
  double W2 = (m_rnu+m_fin)*(m_rnu+m_fin);
  double c = 0.5*(W2+mm_lep-mm_tar*(W2-mm_lep)/s);
  double sqrtD = TMath::Sqrt(LambdaFUNCTION(1.0, mm_lep/s, W2/s));
  double tmp = 0.5*(s-mm_tar);
  double Q2_min = tmp*(1.0-sqrtD)-c;
  double Q2_max = tmp*(1.0+sqrtD)-c;
  // if the nucleus is hydrogen or deuterium then there is no need in further calculation
  if (E_BIN == 0 && P_Fermi == 0)
  {
    Q2_min= TMath::Max(Q2_min,0.0);
    Range1D_t R(Q2_min,Q2_max);
    return R;
  }
  double F_MIN, Q2_0;
  bool LLM;
  // C++ analog of fortran function DMINFC(Q2lim1_SM,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM)
  Func1D<SmithMonizUtils> Q2lim1_SM_(*this, &SmithMonizUtils::Q2lim1_SM);
  DMINFC(Q2lim1_SM_,Q2_min,Q2_max,EPS,Delta,Q2_0,F_MIN,LLM);
  if (F_MIN>0)
  {
     LOG("SmithMoniz", pFATAL)
     << "No overlapped area for energy " << E_nu << "\n" <<
     "Q2_min=" << Q2_min << " Q2_max=" << Q2_max << "\n" <<
     "Q2_0=" << Q2_0 << " F_MIN=" << F_MIN;
     exit(1);
  }
  if (Q2lim1_SM(Q2_min)>0)
  {
    //C++ analog of fortran function Q2_RF=DZEROX(Q2_min,Q2_0,EPS,MFC,Q2lim1_SM,1)
    ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
    rfgb.Solve(Q2lim1_SM_, Q2_min, Q2_0, MFC, EPSABS, EPSREL); //convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
    double Q2_RF = rfgb.Root();
    Q2_min= TMath::Max(Q2_min,Q2_RF);
  }
  if(Q2lim1_SM(Q2_max)>0)
  {
     // C++ analog of fortran function Q2_RF=DZEROX(Q2_0,Q2_max,Eps,MFC,Q2lim1_SM,1)
     ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
     rfgb.Solve(Q2lim1_SM_, Q2_0, Q2_max, MFC, EPSABS, EPSREL); //convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
     double Q2_RF = rfgb.Root();
     Q2_max= TMath::Min(Q2_max,Q2_RF);
  }
  Func1D<SmithMonizUtils> Q2lim2_SM_(*this, &SmithMonizUtils::Q2lim2_SM);
  if (Q2lim2_SM(Q2_min)>0)
  {
     if(Q2lim2_SM(Q2_max)>0)
     {
             LOG("SmithMoniz", pFATAL) << "No overlapped area for energy " << E_nu;
             exit(1);
     }
     else
     {
       // C++ analog of fortran function Q2_RF = DZEROX(Q2_min,Q2_max,Eps,MFC,Q2lim2_SM,1)
       ROOT::Math::RootFinder rfgb(ROOT::Math::RootFinder::kGSL_BRENT);
       rfgb.Solve(Q2lim2_SM_, Q2_min,Q2_max, MFC, EPSABS, EPSREL); //convergence is reached using  tolerance = 2 *( epsrel * abs(x) + epsabs)
       double Q2_RF = rfgb.Root();
       Q2_min= TMath::Max(Q2_min,Q2_RF);
     }
  }

  Q2_min= TMath::Max(Q2_min,0.0);
  Range1D_t R(Q2_min,Q2_max);
  return R;

}
//____________________________________________________________________________
// Return allowed v-range for given Q2
Range1D_t SmithMonizUtils::vQES_SM_lim(double Q2) const
{
  double s = 2*E_nu*m_tar+mm_tar;
  double nu_min= ((m_rnu+m_fin)*(m_rnu+m_fin)+Q2-mm_tar)/(2*m_tar);
  double nu_max= (s-mm_tar-mm_lep*(s-mm_tar)/(Q2+mm_lep)-mm_tar*(Q2+mm_lep)/(s-mm_tar))/(2*m_tar);
  double E = TMath::Sqrt(P_Fermi*P_Fermi+mm_ini);
  double b = (E-E_BIN)*(E-E_BIN)-P_Fermi*P_Fermi;
  double a = (Q2+mm_fin-b)*0.5;
  double tmp1 = a*(E-E_BIN);
  double tmp2  = P_Fermi*TMath::Sqrt(a*a+Q2*b);
  double nu_1 = (tmp1-tmp2)/b;
  double nu_2 = (tmp1+tmp2)/b;
  nu_min= TMath::Max(nu_min,nu_1);
  nu_max= TMath::Min(nu_max,nu_2);
  Range1D_t R;
  if (E_BIN == 0 && P_Fermi == 0)
    nu_max=nu_min;
  if (nu_min<=nu_max)
    R = Range1D_t(nu_min,nu_max);
  else
    R = Range1D_t(0.5*(nu_min+nu_max),0.5*(nu_min+nu_max));
  return R;

}
//____________________________________________________________________________
// Return allowed Fermi momentum range for given Q2 and v
Range1D_t SmithMonizUtils::kFQES_SM_lim(double Q2, double nu) const
{
  double qv = TMath::Sqrt(nu*nu+Q2);
  double c_f = (nu-E_BIN)/qv;
  double d_f = (E_BIN*E_BIN-2*nu*E_BIN-Q2+mm_ini-mm_fin)/(2*qv*m_ini);
  double Ef_min= m_ini*(c_f*d_f+TMath::Sqrt(1.0-c_f*c_f+d_f*d_f))/(1.0-c_f*c_f);
  double kF_min= TMath::Sqrt(TMath::Max(Ef_min*Ef_min-mm_ini,0.0));
  double kF_max= P_Fermi;
  Range1D_t R;
  if (kF_min<=kF_max)
    R = Range1D_t(kF_min,kF_max);
  else
    R = Range1D_t(0.5*(kF_min+kF_max),0.5*(kF_min+kF_max));
  return R;

}
//____________________________________________________________________________
double SmithMonizUtils::LambdaFUNCTION(double a, double b, double c) const
{
  return a*a+b*b+c*c-2*(a*b+b*c+a*c);
}
//____________________________________________________________________________
// C++ implementation of DMINC function from CERN library
void SmithMonizUtils::DMINFC(Func1D<SmithMonizUtils> &F, double A,double B, double EPS, double DELTA, double &X, double &Y, bool &LLM) const
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

  if (T_Fermi==0)                             //Pure Fermi gaz with T_Fermi=0
    if(p<=P_Fermi) 
      return 1.0;
    else
      return 0.0;
  else                                        //Fermi-Dirac distribution
      return 1.0/(1.0 + TMath::Exp(-(P_Fermi-p)/T_Fermi));


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
double SmithMonizUtils::vQES_SM_low_bound(double Q2) const
{
        return vQES_SM_lim(Q2).min;
}
//____________________________________________________________________________
double SmithMonizUtils::vQES_SM_upper_bound(double Q2) const
{
  return -1.0*vQES_SM_lim(Q2).max;
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
   else if (ps == kPSQ2vfE)
   {
     Range1D_t rQ2 = Q2QES_SM_lim();
     const int kNQ2 = 101;
     double dQ2 = (rQ2.max-rQ2.min)/(kNQ2-1);
     for(int iq2=0; iq2<kNQ2; iq2++)
     {
       double Q2 = (rQ2.min + iq2*dQ2);
       Range1D_t rvlims  = vQES_SM_lim(Q2);
       double dv = (rvlims.max-rvlims.min);
       vol += (dQ2*dv);
     }
     return vol;
   }
   else
     return 0;
}
