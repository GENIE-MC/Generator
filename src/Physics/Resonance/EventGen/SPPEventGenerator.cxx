//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
           Vadim Naumov <vnaumov@theor.jinr.ru >,  Joint Institute for Nuclear Research \n
          
  
  For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>

#include <vector>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/BaryonResUtils.h"
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
#include "Physics/Resonance/EventGen/SPPEventGenerator.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
SPPEventGenerator::SPPEventGenerator() :
KineGeneratorWithCache("genie::SPPEventGenerator")
{

}
//___________________________________________________________________________
SPPEventGenerator::SPPEventGenerator(string config) :
KineGeneratorWithCache("genie::SPPEventGenerator", config)
{
  
}
//___________________________________________________________________________
SPPEventGenerator::~SPPEventGenerator()
{
  
}
//___________________________________________________________________________
void SPPEventGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
    
  
  LOG("SPPEventGen", pINFO) << "Generating resonance single pion production event kinematics...";

  if(fGenerateUniformly) {
    LOG("SPPEventGen", pNOTICE)
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
  genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E * f   = new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fXSecModel, interaction, fWcut);
  
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
         LOG("SPPEventGen", pWARN)
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
  Range1D_t Wl  = kps.WLim_SPP_iso();
  if (fWcut >= Wl.min)
    Wl.max = TMath::Min(fWcut,Wl.max);
  Range1D_t Q2l = kps.Q2Lim_W_SPP_iso();
  double W  = Wl.min + (Wl.max - Wl.min)*xin[0];
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
  
  
  double W2 = W*W;
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
  
  
  // boost from isobaric frame to hit nucleon rest frame
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
int SPPEventGenerator::GetRecoilNucleonPdgCode(Interaction * interaction) const
{
   const XclsTag & xcls = interaction->ExclTag();
   if (xcls.NProtons() == 1)
     return kPdgProton;
   else
     return kPdgNeutron;

}
//___________________________________________________________________________
int SPPEventGenerator::GetFinalPionPdgCode(Interaction * interaction) const
{
   const XclsTag & xcls = interaction->ExclTag();
   if (xcls.NPiPlus() == 1)
     return kPdgPiP;
   else if (xcls.NPiMinus() == 1)
     return kPdgPiM;
   return kPdgPi0;

}
//___________________________________________________________________________
void SPPEventGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SPPEventGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void SPPEventGenerator::LoadConfig(void)
{
  
    // Safety factor for the maximum differential cross section
  this->GetParamDef("MaxXSec-SafetyFactor", fSafetyFactor, 1.03);
  this->GetParamDef("Maximum-Depth", fMaxDepth, 3);

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
  
  this->GetParamDef("Wcut", fWcut, -1.);
  

}
//____________________________________________________________________________
double SPPEventGenerator::ComputeMaxXSec(const Interaction * interaction) const
{
   KPhaseSpace * kps = interaction->PhaseSpacePtr();
   Range1D_t Wl = kps->WLim_SPP_iso();
   ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
   ROOT::Math::IBaseFunctionMultiDim * f = new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fXSecModel, interaction, fWcut);
   min->SetFunction( *f );
   min->SetMaxFunctionCalls(10000);  // for Minuit/Minuit2
   min->SetMaxIterations(10000);     // for GSL
   min->SetTolerance(1e-3);
   min->SetPrintLevel(0);
   double min_xsec = 0.;
   double xsec;
   double step = 1e-7;
   // a heuristic algorithm for maximum search
   int total_cells = (TMath::Power(16, fMaxDepth) - 1)/15;
   vector<Cell> cells(total_cells);
   
   for (int dep = 0; dep < fMaxDepth; dep++)
   {
     int aux = TMath::Power(16, dep) - 1;
     for (int cell = aux/15; cell <= 16*aux/15 ; cell++)
     {
       if (cell == 0)
       {
          cells[cell].Vertex1 = Vertex(0., 0., 0., 0.);
          cells[cell].Vertex2 = Vertex(1., 1., 1., 1.);
       }
       double x1m = (cells[cell].Vertex1.x1 + cells[cell].Vertex2.x1)/2;
       double x2m = (cells[cell].Vertex1.x2 + cells[cell].Vertex2.x2)/2;
       double x3m = (cells[cell].Vertex1.x3 + cells[cell].Vertex2.x3)/2;
       double x4m = (cells[cell].Vertex1.x4 + cells[cell].Vertex2.x4)/2;
       min->SetVariable(0, "x1", x1m, step);
       min->SetVariable(1, "x2", x2m, step);
       min->SetVariable(2, "x3", x3m, step);
       min->SetVariable(3, "x4", x4m, step);
       min->SetVariableLimits(0, cells[cell].Vertex1.x1, cells[cell].Vertex2.x1);
       min->SetVariableLimits(1, cells[cell].Vertex1.x2, cells[cell].Vertex2.x2);
       min->SetVariableLimits(2, cells[cell].Vertex1.x3, cells[cell].Vertex2.x3);
       min->SetVariableLimits(3, cells[cell].Vertex1.x4, cells[cell].Vertex2.x4);
       min->Minimize();
       xsec = min->MinValue();
       if (xsec < min_xsec)
         min_xsec = xsec;
       const double *xs = min->X();
       Vertex minv(xs[0], xs[1], xs[2], xs[3]);
       if (minv == cells[cell].Vertex1 || minv == cells[cell].Vertex2)
          minv = Vertex (x1m, x2m, x3m, x4m);
       if (dep < fMaxDepth - 1)
         for (int i = 0; i < 16; i++)
         {
            int child = 16*cell + i + 1;
            cells[child].Vertex1 = minv;
            cells[child].Vertex2 = Vertex ((i>>0)%2?cells[cell].Vertex1.x1:cells[cell].Vertex2.x1,
                                           (i>>1)%2?cells[cell].Vertex1.x2:cells[cell].Vertex2.x2,
                                           (i>>2)%2?cells[cell].Vertex1.x3:cells[cell].Vertex2.x3,
                                           (i>>3)%2?cells[cell].Vertex1.x4:cells[cell].Vertex2.x4);
         }
     }
  }
  Resonance_t res = genie::utils::res::FromString("P33(1232)");
  const InitialState & init_state = interaction -> InitState();
  double Enu = init_state.ProbeE(kRfHitNucRest);
  // other heuristic algorithm for maximum search to fix flaws of the first
  int N3 = 2;
  int N4 = 4;
  double x2max;
  if (Enu < 1.)
    x2max = 1.;
  else 
    x2max = 1./3;
  double dW = Wl.max - Wl.min;  
  double MR  = utils::res::Mass(res);
  double WR  = utils::res::Width(res);
  double x1 = (MR - Wl.min)/dW;
  double x1min = (MR - WR - Wl.min)/dW;
  if (x1min > 1)
  {
     delete f;
     return -min_xsec;
  }
  x1min = x1min<0?0:x1min;
  double x1max = (MR + WR - Wl.min)/dW;
  if (x1max < 0)
  {
     delete f;
     return -min_xsec;
  }
  x1max = x1max>1?1:x1max;
  if (x1 < x1min || x1 > x1max) x1=0.5*(x1min + x1max);
  for (int i3 = 0; i3 < N3; i3++)
  {
    double x3 = 1.*i3;
    double x3min = .5*i3;
    double x3max = .5*(i3 + 1);
    for (int i4 = 0; i4 <= N4; i4++)
    {
      double x4 = 1.*i4/N4;
      double x4min = 1.*i4/N4;
      double x4max = 1.*(i4 + 1)/N4;
      if (i4 == N4)
      {
        x4min = 3./N4;
        x4max = 1;
      }
      min->SetVariable(0, "x1", x1, step);
      min->SetVariable(1, "x2", 1./6, step);
      min->SetVariable(2, "x3", x3, step);
      min->SetVariable(3, "x4", x4, step);
      min->SetVariableLimits(0, x1min, x1max);
      min->SetVariableLimits(1, 0, x2max);
      min->SetVariableLimits(2, x3min, x3max);
      min->SetVariableLimits(3, x4min, x4max);
      min->Minimize();
      xsec = min->MinValue();
      if (xsec < min_xsec)
        min_xsec = xsec;
    }
  }  
  
  delete f;

  return -fSafetyFactor*min_xsec;
}
//____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::d4XSecMK_dWQ2CosThetaPhi_E(
     const XSecAlgorithmI * m, const Interaction * interaction, double wcut) :
ROOT::Math::IBaseFunctionMultiDim(), fModel(m), fWcut(wcut)
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


  if (Enu < kps->Threshold_SPP_iso())
  {
    isZero = true;
    return;
  }
  
  Wl  = kps->WLim_SPP_iso();
  if (fWcut >= Wl.min)
    Wl.max = TMath::Min(fWcut,Wl.max);
  
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
   
  Range1D_t Q2l = kps->Q2Lim_W_SPP_iso(); 
  double Q2 = Q2l.min + (Q2l.max - Q2l.min)*xin[1];
  fInteraction->KinePtr()->SetQ2(Q2);
  
  fInteraction->KinePtr()->SetKV(kKVctp, -1. + 2.*xin[2]); // cosine of pion theta in resonance rest frame
  
  fInteraction->KinePtr()->SetKV(kKVphip , 2.*kPi*xin[3]); // pion phi in resonance rest frame
    
  double xsec = fModel->XSec(fInteraction, kPSWQ2ctpphipfE);
  xsec *= 4*kPi*(Wl.max - Wl.min)*(Q2l.max - Q2l.min);
  return -xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::Clone() const
{
  return
    new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fModel,fInteraction, fWcut);
}
