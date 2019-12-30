//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
           Konstantin Kuzmin <kkuzmin@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
           Vadim Naumov <vnaumov@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
          
  
  For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>

#include <vector>


#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Conventions/Constants.h"
#include "Physics/Resonance/EventGen/RSPPEventGenerator.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
RSPPEventGenerator::RSPPEventGenerator() :
KineGeneratorWithCache("genie::RSPPEventGenerator")
{

}
//___________________________________________________________________________
RSPPEventGenerator::RSPPEventGenerator(string config) :
KineGeneratorWithCache("genie::RSPPEventGenerator", config)
{
  
}
//___________________________________________________________________________
RSPPEventGenerator::~RSPPEventGenerator()
{
  
}
//___________________________________________________________________________
void RSPPEventGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
    
  
  LOG("RSPPEventGen", pINFO) << "Generating resonance single pion production event kinematics...";

  if(fGenerateUniformly) {
    LOG("RSPPEventGen", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }
  
  //-- Get the interaction from the GHEP record
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction -> InitState();
  const KPhaseSpace& kps = interaction->PhaseSpace();

  // Access the target from the interaction summary
  //const Target & tgt = init_state.Tgt();
  Target * tgt = interaction -> InitStatePtr()->TgtPtr();
  
  // get masses of nucleon and pion
  PDGLibrary * pdglib = PDGLibrary::Instance();
  // imply isospin symmetry  
  double mpi  = (pdglib->Find(kPdgPiP)->Mass() + pdglib->Find(kPdgPi0)->Mass() + pdglib->Find(kPdgPiM)->Mass())/3;
  double mpi2 = mpi*mpi; 
  double M = (pdglib->Find(kPdgProton)->Mass() + pdglib->Find(kPdgNeutron)->Mass())/2;
  double M2  = M*M;
  // mass of final lepton
  double ml   = interaction->FSPrimLepton()->Mass();
  double ml2  = ml*ml;
  
  // 4-momentum of neutrino in lab frame
  TLorentzVector k1(*(init_state.GetProbeP4(kRfLab)));
  // 4-momentum of hit nucleon in lab frame
  TLorentzVector p1(*(evrec->HitNucleon())->P4());
  
  TLorentzVector p1_copy(p1);
  
  // set temporarily initial nucleon on shell
  p1.SetE(TMath::Sqrt(p1.P()*p1.P() + M*M));
  
  tgt->SetHitNucP4(p1);
  
  // neutrino 4-momentun in nucleon rest frame
  TLorentzVector k1_HNRF = k1;
  k1_HNRF.Boost(-p1.BoostVector());
  // neutrino energy in nucleon rest frame
  double Ev = k1_HNRF.E();
  // initial nucleon 4-momentun in nucleon rest frame
  TLorentzVector p1_HNRF(0,0,0,M);
  
  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();
  
  // function gives differential cross section and depends on reduced variables W,Q2,cos(theta) and phi -> 1
  genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E * f   = new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fXSecModel, interaction);
  
  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  double xsec = -1;
  double xin[4];
  
  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);
 
  // generate W, Q2, cos(theta) and phi by accept-reject method
  unsigned int iter = 0;
  bool accept = false;
  while(1) 
  {
     iter++;
     if(iter > 100*kRjMaxIterations) {
         LOG("RSPPEventGen", pWARN)
              << "*** Could not select a valid kinematics variable after "
                                                    << iter << " iterations";
         evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
         genie::exceptions::EVGThreadException exception;
         exception.SetReason("Couldn't select kinematics");
         exception.SwitchOnFastForward();
         throw exception;
     }

     xin[0] = rnd->RndKine().Rndm();
     xin[1] = rnd->RndKine().Rndm();
     xin[2] = rnd->RndKine().Rndm();
     xin[3] = rnd->RndKine().Rndm();
     
     
     //-- Computing cross section for the current kinematics
     xsec = -(*f)(xin);
     
     //-- Decide whether to accept the current kinematics
     if(!fGenerateUniformly)
     {
       this->AssertXSecLimits(interaction, xsec, xsec_max);
       double t = xsec_max * rnd->RndKine().Rndm();
       accept = (t < xsec);
     }
     else 
     {
        accept = (xsec>0);
     }
     
     // If the generated kinematics are accepted, finish-up module's job
     if(accept)
     {
       // reset 'trust' bits
       interaction->ResetBit(kISkipProcessChk);
       interaction->ResetBit(kISkipKinematicChk);
       break;
     }
     iter++;
   }

  // W,Q2,cos(theta) and phi from reduced variables
  Range1D_t Wl  = kps.WLim_RSPP();
  Range1D_t Q2l = kps.Q2Lim_W_RSPP();
  double W  = Wl.min + (Wl.max - Wl.min)*xin[0];
  double W2 = W*W;
  interaction->KinePtr()->SetW(W);
  double Q2 = Q2l.min + (Q2l.max - Q2l.min)*xin[1];
  double CosTheta_isb = -1. + 2.*xin[2];
  double SinTheta_isb = (1 - CosTheta_isb*CosTheta_isb)<0?0:TMath::Sqrt(1 - CosTheta_isb*CosTheta_isb);
  double Phi_isb = 2*kPi*xin[3];
  
  // compute x,y for selected W,Q2
  double x=-1, y=-1;
  kinematics::WQ2toXY(Ev,M,W,Q2,x,y);
  
  if(fGenerateUniformly)
  {
    double vol     = (Wl.max-Wl.min)*(Q2l.max-Q2l.min)*4*kPi;
    double totxsec = evrec->XSec();
    double wght    = (vol/totxsec)*xsec;
    wght *= evrec->Weight();
    evrec->SetWeight(wght);
  }


  // set the cross section for the selected kinematics
  evrec->SetDiffXSec(xsec,kPSWQ2ctpphipfE);

  // lock selected kinematics & clear running values
  interaction->KinePtr()->SetQ2(Q2, true);
  interaction->KinePtr()->SetW (W,  true);
  interaction->KinePtr()->Setx (x,  true);
  interaction->KinePtr()->Sety (y,  true);
  interaction->KinePtr()->ClearRunningValues();
  
  
  
  // Kinematical values of all participating particles in the isobaric frame
  double Enu_isb = (Ev*M - (ml2 + Q2)/2)/W;
  double El_isb  = (Ev*M - (ml2 + W2 - M2)/2)/W;
  double v_isb   = (W2 - M2 - Q2)/2/W;
  double q_isb   = TMath::Sqrt(v_isb*v_isb + Q2);
  double kz1_isb = 0.5*(q_isb+(Enu_isb*Enu_isb - El_isb*El_isb + ml2)/q_isb);
  double kz2_isb = kz1_isb - q_isb;
  double kx1_isb = (Enu_isb*Enu_isb - kz1_isb*kz1_isb)<0?0:TMath::Sqrt(Enu_isb*Enu_isb - kz1_isb*kz1_isb);
  double Epi_isb = (W2 + mpi2 - M2)/W/2;
  double ppi_isb = (Epi_isb - mpi)<0?0:TMath::Sqrt(Epi_isb*Epi_isb - mpi2);
  double E1_isb  = (W2 + Q2 + M2)/2/W;
  double E2_isb  = W - Epi_isb;
  
  // 4-momentum of all particles in the isobaric frame
  TLorentzVector k1_isb(kx1_isb, 0, kz1_isb, Enu_isb);
  TLorentzVector k2_isb(kx1_isb, 0, kz2_isb, El_isb);
  TLorentzVector p1_isb(0, 0, -q_isb, E1_isb);
  TLorentzVector p2_isb(-ppi_isb*SinTheta_isb*TMath::Cos(Phi_isb), -ppi_isb*SinTheta_isb*TMath::Sin(Phi_isb), -ppi_isb*CosTheta_isb, E2_isb);
  TLorentzVector pi_isb( ppi_isb*SinTheta_isb*TMath::Cos(Phi_isb),  ppi_isb*SinTheta_isb*TMath::Sin(Phi_isb),  ppi_isb*CosTheta_isb, Epi_isb);
  
  
  // boost from isobaric frame to hit necleon rest frame
  TVector3 boost = -p1_isb.BoostVector();
  k1_isb.Boost(boost);
  k2_isb.Boost(boost);
  p2_isb.Boost(boost);
  pi_isb.Boost(boost);

  
  // rotation to align 3-momentum of all particles
  TVector3 rot_vect = k1_isb.Vect().Cross(k1_HNRF.Vect());
  double rot_angle =  k1_isb.Vect().Angle(k1_HNRF.Vect());
  k2_isb.Rotate(rot_angle, rot_vect);
  p2_isb.Rotate(rot_angle, rot_vect);
  pi_isb.Rotate(rot_angle, rot_vect);

  // boost to laboratory frame
  boost = p1.BoostVector();
  k2_isb.Boost(boost);
  p2_isb.Boost(boost);
  pi_isb.Boost(boost);
  
  
  tgt->SetHitNucP4(p1_copy);

  TLorentzVector x4l(*(evrec->Probe())->X4());
  // add final lepton
  evrec->AddParticle(interaction->FSPrimLeptonPdg(), kIStStableFinalState, evrec->ProbePosition(), -1, -1, -1, k2_isb, x4l);
  
  GHepStatus_t ist = (tgt->IsNucleus()) ? kIStHadronInTheNucleus : kIStStableFinalState;
  // add final nucleon
  evrec->AddParticle(this->GetRecoilNucleonPdgCode(interaction), ist, evrec->HitNucleonPosition(), -1, -1, -1, p2_isb, x4l);
  // add final pion
  evrec->AddParticle(this->GetFinalPionPdgCode(interaction), ist, evrec->HitNucleonPosition(), -1, -1, -1, pi_isb, x4l);
  
  delete f;
  

  return;

}
//___________________________________________________________________________
int RSPPEventGenerator::GetRecoilNucleonPdgCode(Interaction * interaction) const
{
   const XclsTag & xcls = interaction->ExclTag();
   if (xcls.NProtons() == 1)
     return kPdgProton;
   else
     return kPdgNeutron;

}
//___________________________________________________________________________
int RSPPEventGenerator::GetFinalPionPdgCode(Interaction * interaction) const
{
   const XclsTag & xcls = interaction->ExclTag();
   if (xcls.NPiPlus() == 1)
     return kPdgPiP;
   else if (xcls.NPiMinus() == 1)
     return kPdgPiM;
   return kPdgPi0;

}
//___________________________________________________________________________
void RSPPEventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSPPEventGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RSPPEventGenerator::LoadConfig(void)
{
  // Safety factor for the maximum differential cross section
  this->GetParamDef("MaxXSec-SafetyFactor", fSafetyFactor, 1.25);

  // Minimum energy for which max xsec would be cached, forcing explicit
  // calculation for lower eneries
  this->GetParamDef("Cache-MinEnergy", fEMin, 0.5);


  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  this->GetParamDef("MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999.);
  assert(fMaxXSecDiffTolerance>=0);

  // Generate kinematics uniformly over allowed phase space and compute
  // an event weight?
  this->GetParamDef("UniformOverPhaseSpace", fGenerateUniformly, false);
  
  GetParamDef( "NumOfKnots_W",         N_W,  40);
  GetParamDef( "NumOfKnots_Q2",        N_Q2, 40);
  GetParamDef( "NumOfKnots_CosTheta",  N_CosTheta,  40);
  GetParamDef( "NumOfKnots_Phi",       N_Phi, 40);

}
//____________________________________________________________________________
double RSPPEventGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
  const InitialState & init_state = interaction -> InitState();
  
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
  ROOT::Math::IBaseFunctionMultiDim * f = new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fXSecModel, interaction);
  min->SetFunction( *f );
  min->SetMaxFunctionCalls(10000);  // for Minuit/Minuit2
  min->SetMaxIterations(10000);     // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(0);
  double step[4] = {1e-7,1e-7,1e-7,1e-7};
  double variable[4];
  double xin[4];
  
  // coarse minimum search
  double Fmin = 0.;
  int ixWmin = -1, ixQ2min = -1, ixCosThetamin = -1, ixPhimin = -1;
  for (int ixW = 0; ixW <= N_W; ixW++)
  {
    xin[0] = 1.*ixW/N_W;
    for (int ixQ2 = 0; ixQ2 <= N_Q2; ixQ2++)
    {
      xin[1] = 1.*ixQ2/N_Q2;
      for (int ixCosTheta = 0; ixCosTheta <= N_CosTheta; ixCosTheta++)
      {
        xin[2] = 1.*ixCosTheta/N_CosTheta;
        for (int ixPhi = 0; ixPhi <= N_Phi; ixPhi++)
        {
          xin[3] = 1.*ixPhi/N_Phi;
          double F = (*f)(xin);
          if (F < Fmin)
          {
            Fmin = F;
            ixWmin = ixW;
            ixQ2min = ixQ2;
            ixCosThetamin = ixCosTheta;
            ixPhimin = ixPhi;
          }
        }
      }
    }
  }

  
  //delicate minimum search
  int ixWdl        = TMath::Max(ixWmin - 1,0);
  int ixWul        = TMath::Min(ixWmin + 1,N_W);
  int ixQ2dl       = TMath::Max(ixQ2min - 1,0);
  int ixQ2ul       = TMath::Min(ixQ2min + 1,N_Q2);
  int ixCosThetadl = TMath::Max(ixCosThetamin - 1,0);
  int ixCosThetaul = TMath::Min(ixCosThetamin + 1,N_CosTheta);
  int ixPhidl      = TMath::Max(ixPhimin - 1,0);
  int ixPhiul      = TMath::Min(ixPhimin + 1,N_Phi);
    
  variable[0] = 1.*ixWmin/N_W;
  variable[1] = 1.*ixQ2min/N_Q2;
  variable[2] = 1.*ixCosThetamin/N_CosTheta;
  variable[3] = 1.*ixPhimin/N_Phi;
  min->SetVariable(0,"W",         variable[0], step[0]);
  min->SetVariable(1,"Q2",        variable[1], step[1]);
  min->SetVariable(2,"CosThetaPi",variable[2], step[2]);
  min->SetVariable(3,"PhiPi",     variable[3], step[3]);
  min->SetVariableLimits(0, 1.*ixWdl/N_W,               1.*ixWul/N_W);
  min->SetVariableLimits(1, 1.*ixQ2dl/N_Q2,             1.*ixQ2ul/N_Q2);
  min->SetVariableLimits(2, 1.*ixCosThetadl/N_CosTheta, 1.*ixCosThetaul/N_CosTheta);
  min->SetVariableLimits(3, 1.*ixPhidl/N_Phi,           1.*ixPhiul/N_Phi);
  min->Minimize();
  double min_xsec = min->MinValue();
  
  if (min_xsec > Fmin)
    min_xsec = Fmin;
     
  delete f;


  return -min_xsec;
}
//____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::d4XSecMK_dWQ2CosThetaPhi_E(
     const XSecAlgorithmI * m, const Interaction * interaction) :
ROOT::Math::IBaseFunctionMultiDim(), fModel(m)
{

  isZero = false;
  fInteraction = const_cast<Interaction*>(interaction);
  // skip process and kinematic checks
  fInteraction->SetBit(kISkipProcessChk);
  fInteraction->SetBit(kISkipKinematicChk);
  
  kps = fInteraction->PhaseSpacePtr();
  
  // Get kinematical parameters
  const InitialState & init_state = interaction -> InitState();
  double Enu = init_state.ProbeE(kRfHitNucRest);


  if (Enu < kps->Threshold_RSPP())
  {
    isZero = true;
    return;
  }
  
  Wl  = kps->WLim_RSPP();

}
genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::~d4XSecMK_dWQ2CosThetaPhi_E()
{

}
unsigned int genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::NDim(void) const
{
  return 4;
}
double genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::DoEval(const double * xin) const
{

// outputs:
//   differential cross section [10^-38 cm^2/GeV^3] for resonance single pion production production
//
  if (isZero) return 0.;
  
  double W  = Wl.min + (Wl.max - Wl.min)*xin[0];
  fInteraction->KinePtr()->SetW(W);
   
  Range1D_t Q2l = kps->Q2Lim_W_RSPP(); 
   
  double Q2 = Q2l.min + (Q2l.max - Q2l.min)*xin[1];
  fInteraction->KinePtr()->SetQ2(Q2);
  
  fInteraction->KinePtr()->SetKV(kKVctp, -1. + 2.*xin[2]); // cosine of pion theta in resonance rest frame
  
  fInteraction->KinePtr()->SetKV(kKVphip , 2.*kPi*xin[3]); // pion phi in resonance rest frame
    
  double xsec = -fModel->XSec(fInteraction, kPSWQ2ctpphipfE)*(Wl.max-Wl.min)*(Q2l.max-Q2l.min)*4*kPi;
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::Clone() const
{
  return
    new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fModel,fInteraction);
}

