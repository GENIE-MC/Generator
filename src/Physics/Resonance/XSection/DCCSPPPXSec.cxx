//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research

  For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TSystem.h>
#include <TMath.h>
#include <Math/SpecFuncMathMore.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/Units.h"
#include <Framework/Conventions/KinePhaseSpace.h>
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Physics/Resonance/XSection/DCCSPPPXSec.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/Resonance/XSection/DCCSPPPXSec.h"

#include <fstream>
#include <limits>

using namespace genie;
using namespace genie::constants;
using namespace genie::units;
using namespace std::complex_literals;

//____________________________________________________________________________
DCCSPPPXSec::DCCSPPPXSec() :
XSecAlgorithmI("genie::DCCSPPPXSec")
{

}
//____________________________________________________________________________
DCCSPPPXSec::DCCSPPPXSec(std::string config) :
XSecAlgorithmI("genie::DCCSPPPXSec", config)
{

}
//____________________________________________________________________________
DCCSPPPXSec::~DCCSPPPXSec()
{

}
//____________________________________________________________________________
double DCCSPPPXSec::XSec(const Interaction * interaction, KinePhaseSpace_t kps) const
{
  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const InitialState & init_state = interaction -> InitState();
  const ProcessInfo &  proc_info  = interaction -> ProcInfo();
  const Target & target = init_state.Tgt();

  double xsec = 0;
  //-- Get 1pi exclusive channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);

  int  nucpdgc   = target.HitNucPdg();
//  int  probepdgc = init_state.ProbePdg();
  int  helicity  = init_state.ProbeHelicity();
//  bool is_nubar  = pdg::IsAntiNeutrino     (probepdgc);
    bool is_EM     = proc_info.IsEM();
//  bool is_CC     = proc_info.IsWeakCC();
//  bool is_NC     = proc_info.IsWeakNC();
  bool is_p      = pdg::IsProton  (nucpdgc);


  // Get kinematical parameters
  const Kinematics & kinematics = interaction -> Kine();
  double Eli_L = init_state.ProbeE(kRfHitNucRest);
  double ml    = interaction->FSPrimLepton()->Mass();
  double ml2   = ml*ml;
  double Q2    = kinematics.Q2();
  double W     = kinematics.W();
  double W2    = W*W;
  
  // dimension of kine phase space
  std::string s = KinePhaseSpace::AsString(kps);
  int kpsdim = s!="<|E>"?1 + std::count(s.begin(), s.begin()+s.find('}'), ','):0;
  if (kpsdim < 3 || kpsdim > 4) return 0.;
  // Pion angles should be given in pi-N rest frame
  double cos_theta_pi = kinematics.GetKV(kKVctp);
  double phi_pi = 0.;
  if (kpsdim == 4)
    phi_pi = kinematics.GetKV(kKVphip);

  double sin_theta_pi   = TMath::Sqrt(1. - cos_theta_pi*cos_theta_pi);
//  double SinHalfTheta  = TMath::Sqrt((1 - CosTheta_pi)/2);
//  double CosHalfTheta  = TMath::Sqrt((1 + CosTheta_pi)/2);

  PDGLibrary * pdglib = PDGLibrary::Instance();
  double Mi    = pdglib->Find(SppChannel::InitStateNucleon(spp_channel))->Mass();
  double Mi2   = Mi*Mi;
  double Mf    = pdglib->Find(SppChannel::FinStateNucleon(spp_channel))->Mass();
  double Mf2   = Mf*Mf;
  double m_pi  = pdglib->Find(SppChannel::FinStatePion(spp_channel))->Mass();
  double m_pi2 = m_pi*m_pi;


  // Eq. 14 of ref. 4
  double nu                = (W2 - Mi2 - Q2)/2./W;
  double nu_L              = nu*W/Mi;
  // Eq. 15 of ref. 4
  double q                 = TMath::Sqrt(Q2 + nu*nu);
  double qL                = TMath::Sqrt(Q2 + nu_L*nu_L);
  
  double k0                = (W2 - Mf2 + m_pi2)/2./W;
  // Eq. 16 of ref. 4
  double k                 = TMath::Sqrt(k0*k0 - m_pi2);
  // Eq. 17 of ref. 4
  double q_gamma           = (W2 - Mi2)/2./W;
  // Eq. 3 of ref. 4
  double qL_gamma          = (W2 - Mi2)/2./Mi;
  
  double Elf_L             = Eli_L - nu_L;
  double pli_L             = TMath::Sqrt(Eli_L*Eli_L - ml2);
  double plf_L             = TMath::Sqrt(Elf_L*Elf_L - ml2);
  double cos_theta_l       = (Eli_L*Elf_L - Q2/2 - ml2)/pli_L/plf_L;
  // Eq. 4 of ref. 4
  double epsilon = 0.;
  if (cos_theta_l > -1.)
  {
    double tan_half_theta_l2  = (1. - cos_theta_l)/(1. + cos_theta_l);
    epsilon                   = 1./(1. + 2.*qL*qL*tan_half_theta_l2/Q2);
  }
  // Eq. 2 of ref. 4
  double Gamma             = kAem*qL_gamma*Elf_L/2./kPi2/Q2/Eli_L/(1. - epsilon);
    
  
  std::complex<double> F1, F2,F3, F4, F7, F8;
  for (unsigned int L = 0; L<5; L++)
  {
     VAmpl vampl = Amplitudes(W, Q2, L, spp_channel);
     // Multipole amplitudes, in units (Fermi/1000)
     std::complex<double> & ELp =  vampl[0];
     std::complex<double> & ELm =  vampl[1];
     std::complex<double> & MLp =  vampl[2];
     std::complex<double> & MLm =  vampl[3];
     std::complex<double> & SLp =  vampl[4];
     std::complex<double> & SLm =  vampl[5];
     double & x = cos_theta_pi;

     // Eq. 26 of ref. 4
     F1 += dPdx(L+1, x)*(ELp + 1.*L*MLp) + dPdx(L-1, x)*(ELm + 1.*(L+1)*MLm);
     // Eq. 27 of ref. 4
     F2 += dPdx(L, x)*(1.*(L+1)*MLp + 1.*L*MLm);
     // Eq. 28 of ref. 4
     F3 += d2Pdx2(L+1, x)*(ELp - MLp) + d2Pdx2(L-1, x)*(ELm + MLm);
     // Eq. 29 of ref. 4
     F4 += d2Pdx2(L, x)*(MLp - MLm - ELp - ELm);
     // Eq. 32 of ref. 4
     F7 += dPdx(L, x)*(1.*L*SLm - 1.*(L+1)*SLp) ;
     // Eq. 33 of ref. 4
     F8 += (L+1)*dPdx(L+1, x)*SLp - L*dPdx(L-1, x)*SLm;      
  }
  
  // The Pauli matrices times by imaginary unit between the final and initial nucleon spins along z-axis 
  std::vector<std::vector<std::complex<double> > > ISx = { {0., 1i}, {1i, 0.} };
  std::vector<std::vector<std::complex<double> > > ISy = { {0., 1.}, {-1., 0.} };
  std::vector<std::vector<std::complex<double> > > ISz = { {1i, 0.}, {0., -1i} };
  
  std::vector<std::vector<std::complex<double> > > Fx(2, std::vector<std::complex<double> >(2));
  std::vector<std::vector<std::complex<double> > > Fy(2, std::vector<std::complex<double> >(2));
  std::vector<std::vector<std::complex<double> > > F0(2, std::vector<std::complex<double> >(2));
  
  for (int sf = 0; sf<2; sf++)
  {
    for (int si = 0; si<2; si++)
    {
      // Eq. 23 of ref. 4
      Fx[sf][si] = ISx[sf][si]*(F1 - F2*cos_theta_pi + F4*sin_theta_pi*sin_theta_pi) + 
                   ISz[sf][si]*sin_theta_pi*(F2 + F3 + F4*cos_theta_pi);
      // Eq. 24 of ref. 4
      Fy[sf][si] = ISy[sf][si]*(F1 - F2*cos_theta_pi) - F2*sin_theta_pi;
      // Eq. 25 of ref. 4
      F0[sf][si] = ISz[sf][si]*(F7*cos_theta_pi + F8) - ISy[sf][si]*F2*sin_theta_pi;
    }
  }
  
  // Eq. 21 of ref. 4
  double ReFxx, ReFyy, ReF00, ReFx0, ImFx0;
  std::complex<double> Fxx, Fyy, F00, Fx0;
  for (int sf = 0; sf<2; sf++)
  {
    for (int si = 0; si<2; si++)
    {
      Fxx += Fx[sf][si]*std::conj(Fx[sf][si]);
      Fyy += Fy[sf][si]*std::conj(Fy[sf][si]);
      F00 += F0[sf][si]*std::conj(F0[sf][si]);
      Fx0 += Fx[sf][si]*std::conj(F0[sf][si]);
    }
  }
  ReFxx = 0.5*Fxx.real();
  ReFyy = 0.5*Fyy.real();
  ReF00 = 0.5*F00.real();
  ReFx0 = 0.5*Fx0.real();
  ImFx0 = 0.5*Fx0.imag();
  
  // Eq. 9 of ref. 4
  double dsigmaTdOmega_pi    = k*(ReFxx + ReFyy)/q_gamma/2.;
  // Eq. 10 of ref. 4
  double dsigmaLdOmega_pi    = k*Q2*ReF00/q_gamma/q/q;
  // Eq. 11 of ref. 4
  double dsigmaTTdOmega_pi   = k*(ReFxx - ReFyy)/q_gamma/2.;
  // Eq. 12 of ref. 4
  double dsigmaLTdOmega_pi   = -k*TMath::Sqrt(Q2)*ReFx0/q_gamma/q;
  // Eq. 13 of ref. 4
  double dsigmaLTpdOmega_pi  = k*TMath::Sqrt(Q2)*ImFx0/q_gamma/q;

  double dsigmavdOmega_pi;
  if (kpsdim == 4)
  {
    // Eq. 8 of ref. 4
    dsigmavdOmega_pi = 
       dsigmaTdOmega_pi + epsilon*dsigmaLdOmega_pi +
       epsilon*dsigmaTTdOmega_pi*TMath::Cos(2.*phi_pi) +
       TMath::Sqrt(2*epsilon*(1. + epsilon))*dsigmaLTdOmega_pi*TMath::Cos(phi_pi) + 
       TMath::Sqrt(2*epsilon*(1. - epsilon))*dsigmaLTpdOmega_pi*TMath::Sin(phi_pi)*helicity;
  }
  else if (kpsdim == 3)
  {
    // Eq. 18 of ref. 4
    dsigmavdOmega_pi = dsigmaTdOmega_pi + epsilon*dsigmaLdOmega_pi;
    dsigmavdOmega_pi*= 2*kPi;
  }
  // convert to hbarc=1 units
  dsigmavdOmega_pi*=fermi2*1e-6;
  // Eq. 1 of ref. 4
  xsec = Gamma*dsigmavdOmega_pi;
  // convert to d^4sig/dQ2dWdcos_theta_pidphi_pi
  xsec*=kPi*W/Mi/pli_L/plf_L;
 
  // The algorithm computes d^4xsec/dWdQ2dCostThetadPhi or d^3xsec/dWdQ2dCostTheta
  // Check whether variable tranformation is needed
  if ( kps != kPSWQ2ctpphipfE && kps != kPSWQ2ctpfE )
  {
     double J = 1.;
     if (kpsdim == 3)
       J = utils::kinematics::Jacobian(interaction, kPSWQ2ctpfE, kps);
     else if (kpsdim == 4)
       J = utils::kinematics::Jacobian(interaction, kPSWQ2ctpphipfE, kps);
     xsec *= J;
  }

  // Apply given scaling factor
  if      (is_EM) { xsec *= fXSecScaleEM; }

  // If requested return the free nucleon xsec even for input nuclear tgt
  if ( interaction->TestBit(kIAssumeFreeNucleon) ) return xsec;

  int Z = target.Z();
  int A = target.A();
  int N = A-Z;

  // Take into account the number of scattering centers in the target
  int NNucl = (SppChannel::InitStateNucleon(spp_channel) == kPdgProton) ? Z : N;
  xsec*=NNucl; // nuclear xsec (no nuclear suppression symmetry_factor)

  if ( fUsePauliBlocking && A!=1 && kps == kPSWQ2ctpfE )
  {
    // Calculation of Pauli blocking according to refs. 9-11
    double P_Fermi = 0.0;

    // Maximum value of Fermi momentum of target nucleon (GeV)
    if ( A<6 || ! fUseRFGParametrization )
    {
        // look up the Fermi momentum for this target
        FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
        const FermiMomentumTable * kft = kftp->GetTable(fKFTable);
        P_Fermi = kft->FindClosestKF(pdg::IonPdgCode(A, Z), nucpdgc);
     }
     else {
        // define the Fermi momentum for this target
        P_Fermi = utils::nuclear::FermiMomentumForIsoscalarNucleonParametrization(target);
        // correct the Fermi momentum for the struck nucleon
        if(is_p) { P_Fermi *= TMath::Power( 2.*Z/A, 1./3); }
        else     { P_Fermi *= TMath::Power( 2.*N/A, 1./3); }
     }

     double FactorPauli_RES = 1.0;

     if (P_Fermi > 0.)
     {
        if ( 2*P_Fermi < q-k )
           FactorPauli_RES = 1.0;
        if ( 2*P_Fermi >= q+k )
           FactorPauli_RES = ((3*q*q + k*k)/(2*P_Fermi) - (5*TMath::Power(q,4) + TMath::Power(k,4) + 10*q*q*k*k)/(40*TMath::Power(P_Fermi,3)))/(2*q);
        if ( 2*P_Fermi >= q-k && 2*P_Fermi <= q+k )
           FactorPauli_RES = ((k + q)*(k + q) - 4*P_Fermi*P_Fermi/5 - TMath::Power(q - k, 3)/(2*P_Fermi)+TMath::Power(q - k, 5)/(40*TMath::Power(P_Fermi, 3)))/(4*k*q);
     }

     xsec *= FactorPauli_RES;
  }
  return xsec;

}

//____________________________________________________________________________
double DCCSPPPXSec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool DCCSPPPXSec::ValidProcess(const Interaction * interaction) const
{

  if(interaction->TestBit(kISkipProcessChk)) return true;

  //-- Get the requested SPP channel
  SppChannel_t spp_channel = SppChannel::FromInteraction(interaction);
  if( spp_channel == kSppNull ) {
    return false;
  }

  return true;

}
//____________________________________________________________________________
bool DCCSPPPXSec::ValidKinematics(const Interaction * interaction) const
{
  // call only after ValidProcess
  if ( interaction->TestBit(kISkipKinematicChk) ) return true;

  const KPhaseSpace  & kps        = interaction->PhaseSpace();
  
  // Get kinematical parameters
  const InitialState & init_state = interaction -> InitState();
  const Kinematics & kinematics = interaction -> Kine();
  double Enu = init_state.ProbeE(kRfHitNucRest);
  double W    = kinematics.W();
  double Q2   = kinematics.Q2();

  if (Enu < kps.Threshold())
    return false;
  
  Range1D_t Wl  = kps.WLim_SPP();
  Range1D_t Q2l = kps.Q2Lim_W_SPP();
      
  // model restrictions
  Wl.min  = TMath::Max (Wl.min,  1.08);
  Wl.max  = TMath::Min (Wl.max,  2.00);
  Q2l.min = TMath::Max (Q2l.min, 0.00);
  Q2l.max = TMath::Min (Q2l.max, 3.00);
  
  if (W < Wl.min || W > Wl.max)
    return false;

  if (Q2 < Q2l.min || Q2 > Q2l.max)
    return false;

  return true;

}
//____________________________________________________________________________
void DCCSPPPXSec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPPXSec::Configure(std::string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DCCSPPPXSec::LoadConfig(void)
{

  // Cross section scaling symmetry_factors
  this->GetParam( "DCC-EM-XSecScale", fXSecScaleEM );
  
  GetParamDef( "WarnIfMissing", fWarnIfMissing, true );
 
  // Either a data path relative to the root GENIE folder
  // or an absolute path can be used. Find out which
  // option was chosen.
  std::string path_type;
  GetParamDef( "DataPathType", path_type, std::string("relative") );
 
  // Right now, there can only be a single data path
  // specified. We use a vector of paths to allow for
  // easy expansion later.
  std::string data_path;
  GetParam( "DataPath", data_path );
 
  // Convert the relative path to an absolute one if needed
  if ( path_type == "relative" ) {
    data_path = std::string( gSystem->Getenv("GENIE") ) + '/' + data_path;
  }
 
  fDataPaths.push_back( data_path );

  GetParam("FermiMomentumTable", fKFTable);
  GetParam("RFG-UseParametrization", fUseRFGParametrization);
  GetParam("UsePauliBlockingForRES", fUsePauliBlocking);

  // Load the differential cross section integrator
  fXSecIntegrator = dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
//____________________________________________________________________________
DCCSPPPXSec::CPtrDT DCCSPPPXSec::GetDataTable(SppChannel_t spp_chn) const
{

  // First check to see if the data table object already exists in memory.
  // If it does, return the existing object.
  if ( fTables.count(spp_chn) ) return fTables.find(spp_chn)->second.get();

  // If not, try to create it
  CPtrDT dt = BuildDataTable(spp_chn);

  if ( !dt && fWarnIfMissing ) 
  {
    LOG("DCCSPPPXSec", pWARN) << "Unable to create a data table"
      << " for channel " << SppChannel::AsString(spp_chn);
  }

  return dt;
}
//____________________________________________________________________________
DCCSPPPXSec::CPtrDT DCCSPPPXSec::BuildDataTable(SppChannel_t spp_chn) const
{
  bool table_ok = true;

  std::string table_file_basename = GetDataTableFileBasename( spp_chn );

  // Tables values are represented using a 2D grid that is stored in a data
  // file. Get the full path to the file, or an empty string if it could not
  // be found. Also set the table_ok flag to false if the file could not be
  // found.
  std::string full_file_name = FindDataTableFile(table_file_basename, table_ok);

  if ( table_ok ) 
  {

    // Create the new table object
    LOG("DCCSPPPXSec", pINFO) << "Loading the data table file" << full_file_name;

    fTables[spp_chn] = ParseDataTableFile( full_file_name );

    // Return a pointer to the newly-created table object
    return fTables[spp_chn].get();
  }
  else 
  {
    // If we couldn't make the table, store a nullptr to avoid
    // unsuccessful repeat attempts. These can otherwise slow things down
    // for no good reason.
    fTables[spp_chn] = NULL;

    if ( fWarnIfMissing ) 
    {
      LOG("DCCSPPPXSec", pERROR) << "The table data file \""
        << full_file_name << "\" requested for channel "
        << SppChannel::AsString(spp_chn)  << " could not be found";
    }
 
  }

  // If there was a problem, return a null pointer
  return NULL;
}
//____________________________________________________________________________
DCCSPPPXSec::UPtrDT DCCSPPPXSec::ParseDataTableFile( std::string full_file_name ) const
{
  UPtrDT dt = std::make_unique< std::vector<std::vector<double> > >();
  std::ifstream infile( full_file_name.c_str() );
  if (infile) 
  {
    std::string line;
    unsigned int line_number = 0;
    while (std::getline(infile, line))
    {
      line_number++;
      unsigned int NQ2 = (line_number - 8)/306;
      unsigned int NL = (line_number - 306*NQ2 - 8)/61;
      unsigned int NW = (line_number - 306*NQ2 - 61*NL - 8);
      if (NQ2<=25 && NL<=4 && NW<=58)
      {
        std::vector<double> row(12);
        for (int i=0; i<12; i++)
        {
           row[i] = std::stod(line.substr(14+12*i, 11)); 
        }
        dt->push_back(row);
      } 
    }
  }
  
  return dt;
}
//____________________________________________________________________________
void DCCSPPPXSec::GetTablePos(double W, double Q2, unsigned int L, TablePos & tabpos) const
{
  unsigned int NW;
  if (1.08 <= W && W <= 1.16)
  {
    NW = int(floor(50*W)) - 54;
    if (NW == 50*W - 54)
    {
      tabpos.hi_W = tabpos.lo_W = 1.08 + 0.02*NW;
    } 
    else
    {
      tabpos.lo_W = 1.08 + 0.02*NW;
      tabpos.hi_W = 1.10 + 0.02*NW;
    }
  }
  else if (1.16 < W && W <= 1.24)
  {
    NW = int(floor(200*W)) - 228;
    if (NW == 200*W - 228)
    {
      tabpos.hi_W = tabpos.lo_W = 1.14 + 0.005*NW;
    } 
    else
    {
      tabpos.lo_W = 1.14  + 0.005*NW;
      tabpos.hi_W = 1.145 + 0.005*NW;
    }
  }
  else if (1.24 < W && W <= 2.0)
  {
    NW = int(floor(50*W)) - 42;
    if (NW == 50*W - 42)
    {
      tabpos.hi_W = tabpos.lo_W = 0.84 + 0.02*NW;
    } 
    else
    {
      tabpos.lo_W = 0.84 + 0.02*NW;
      tabpos.hi_W = 0.86 + 0.02*NW;
    }
  }
  else
  {
     LOG("DCCSPPPXSec", pERROR) 
               << "Can't find row in table for given W";
     exit (1);
  }
  
  unsigned int NQ2;  
  if (0.0 <= Q2 && Q2 <= 0.2)
  {
    NQ2 = int(floor(50*Q2));
    if (NQ2 == 50*Q2)
    {
      tabpos.hi_Q2 = tabpos.lo_Q2 = 0.02*NQ2;
    } 
    else
    {
      tabpos.lo_Q2 = 0.02*NQ2;
      tabpos.hi_Q2 = 0.02 + 0.02*NQ2;
    }
  }
  else if (0.2 < Q2 && Q2 <= 1.0)
  {
    NQ2 = int(floor(10*Q2)) + 8;
    if (NQ2 == 10*Q2 + 8)
    {
      tabpos.hi_Q2 = tabpos.lo_Q2 = -0.8 + 0.1*NQ2;
    } 
    else
    {
      tabpos.lo_Q2 = -0.8 + 0.1*NQ2;
      tabpos.hi_Q2 = -0.7 + 0.1*NQ2;
    }
  }
  else if (1.0 < Q2 && Q2 <= 1.6)
  {
    NQ2 = int(floor(5*Q2)) + 13;
    if (NQ2 == 5*Q2 + 13)
    {
      tabpos.hi_Q2 = tabpos.lo_Q2 = -2.6 + 0.2*NQ2;
    } 
    else
    {
      tabpos.lo_Q2 = -2.6 + 0.2*NQ2;
      tabpos.hi_Q2 = -2.4 + 0.2*NQ2;
    }
  }
  else if (1.6 < Q2 && Q2 <= 2.4)
  {
    NQ2 = int(floor(2.5*Q2)) + 17;
    if (NQ2 == 2.5*Q2 + 17)
    {
      tabpos.hi_Q2 = tabpos.lo_Q2 = -6.8 + 0.4*NQ2;
    } 
    else
    {
      tabpos.lo_Q2 = -6.8 + 0.4*NQ2;
      tabpos.hi_Q2 = -6.4 + 0.4*NQ2;
    }
  }
  else if (2.4 < Q2 && Q2 <= 3.0)
  {
    NQ2 = int(floor(Q2/0.3)) + 15;
    if (NQ2 == Q2/0.3 + 15)
    {
      tabpos.hi_Q2 = tabpos.lo_Q2 = -4.5 + 0.3*NQ2;
    } 
    else
    {
      tabpos.lo_Q2 = -4.5 + 0.3*NQ2;
      tabpos.hi_Q2 = -4.2 + 0.3*NQ2;
    }
  }
  else
  {
     LOG("DCCSPPPXSec", pERROR) 
               << "Can't find row in table for given Q2";
     exit (1);
  }
  tabpos.lo_row_W  = 295*NQ2 + 59*L + NW;
  if (tabpos.hi_W!=tabpos.lo_W && tabpos.hi_Q2!=tabpos.lo_Q2)
  {
     tabpos.hi_row_W  = tabpos.lo_row_W + 1;
     tabpos.lo_row_Q2 = tabpos.lo_row_W + 295;
     tabpos.hi_row_Q2 = tabpos.lo_row_Q2 + 1;
  }
  else if (tabpos.hi_W==tabpos.lo_W && tabpos.hi_Q2!=tabpos.lo_Q2)
  {
     tabpos.hi_row_W  = tabpos.lo_row_W;
     tabpos.hi_row_Q2 = tabpos.lo_row_Q2 = tabpos.lo_row_W + 295;
  }
  else if (tabpos.hi_W!=tabpos.lo_W && tabpos.hi_Q2==tabpos.lo_Q2)
  {
     tabpos.hi_row_W  = tabpos.lo_row_W + 1;
     tabpos.lo_row_Q2 = tabpos.lo_row_W;
     tabpos.hi_row_Q2 = tabpos.hi_row_W;
  }
  else
    tabpos.hi_row_W  = tabpos.lo_row_Q2 = tabpos.hi_row_Q2 = tabpos.lo_row_W;
  
}
//____________________________________________________________________________
std::string DCCSPPPXSec::FindDataTableFile(const std::string &basename, bool &ok) const
{

  for (size_t p = 0; p < fDataPaths.size(); ++p) 
  {
    const std::string& path = fDataPaths.at( p );
    std::string full_name = path + '/' + basename;
    // is file exist?
    if ( std::ifstream(full_name.c_str()).good() ) 
      return full_name;
  }

  // A matching file could not be found
  ok = false;
  return std::string();
 
}
//____________________________________________________________________________
std::string DCCSPPPXSec::GetDataTableFileBasename(SppChannel_t spp_chn) const
{
  switch (spp_chn)
  {
    case kSpp_lp_em_10010:
        return "q2-mul-gp-pi0p.dat";
    case kSpp_lp_em_01100:
        return "q2-mul-gp-pipn.dat";
    case kSpp_ln_em_01010:
        return "q2-mul-gn-pi0n.dat";
    case kSpp_ln_em_10001:
        return "q2-mul-gn-pimp.dat";
    default:
        return std::string();
  }
}
//____________________________________________________________________________
DCCSPPPXSec::VAmpl DCCSPPPXSec::Amplitudes(double W, double Q2, unsigned int L, SppChannel_t spp_chn) const
{
  CPtrDT pdt = GetDataTable(spp_chn);
  VAmpl vampl(6);
  TablePos tabpos;
  GetTablePos(W, Q2, L, tabpos);
  
  // bilinear interpolation
  double x1  = tabpos.lo_W;
  double x2  = tabpos.hi_W;
  double y1  = tabpos.lo_Q2;
  double y2  = tabpos.hi_Q2;
  for (int i = 0; i < 6; i++)
  {
    std::complex<double> a = 0.;
    for (int j = 0; j < 2; j++)
    {
       double z;
       if (x1!=x2 && y1!=y2)
       {
         double z11 = (*pdt)[tabpos.lo_row_W][2*i+j];
         double z21 = (*pdt)[tabpos.hi_row_W][2*i+j];
         double z12 = (*pdt)[tabpos.lo_row_Q2][2*i+j];
         double z22 = (*pdt)[tabpos.hi_row_Q2][2*i+j];
         double z1  = z11*(x2 - W)/(x2 - x1) + z21*(W - x1)/(x2 - x1);
         double z2  = z12*(x2 - W)/(x2 - x1) + z22*(W - x1)/(x2 - x1);
         z   = z1*(y2 - Q2)/(y2 - y1) + z2*(Q2 - y1)/(y2 - y1);
       }
       else if (x1==x2 && y1!=y2)
       {
          double z11 = (*pdt)[tabpos.lo_row_W][2*i+j];
          double z12 = (*pdt)[tabpos.lo_row_Q2][2*i+j];
          double z1 = z11;
          double z2 = z12;
          z   = z1*(y2 - Q2)/(y2 - y1) + z2*(Q2 - y1)/(y2 - y1);
       }
       else if (x1!=x2 && y1==y2)
       {
          double z11 = (*pdt)[tabpos.lo_row_W][2*i+j];
          double z21 = (*pdt)[tabpos.hi_row_W][2*i+j];
          double z1 = z11;
          double z2 = z21;
          z   = z1*(x2 - W)/(x2 - x1) + z2*(W - x1)/(x2 - x1);
       }
       else
       {
          double z11 = (*pdt)[tabpos.lo_row_W][2*i+j];
          z = z11;
       }
       a += j?z*std::complex<double>(0., 1.):z;
    }
    vampl[i] = a;
  }
  return vampl;
}
//____________________________________________________________________________
double DCCSPPPXSec::dPdx (int L, double x) const
{
     if (L<=0)
       return 0.;
       
     if (x == 1.0)
       x -= std::numeric_limits<double>::epsilon();
     if (x == -1.0)
       x += std::numeric_limits<double>::epsilon();

     return ROOT::Math::assoc_legendre(L, 1, x)*TMath::Power(1 - x*x, -0.5);
}
//____________________________________________________________________________
double DCCSPPPXSec::d2Pdx2 (int L, double x) const
{
     if (L<=1)
       return 0.;
       
     if (x == 1.0)
       x -= std::numeric_limits<double>::epsilon();
     if (x == -1.0)
       x += std::numeric_limits<double>::epsilon();

     return ROOT::Math::assoc_legendre(L, 2, x)*TMath::Power(1 - x*x, -1.0);
}
//____________________________________________________________________________
