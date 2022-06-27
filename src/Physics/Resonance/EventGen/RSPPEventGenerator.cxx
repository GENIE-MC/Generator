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
  Range1D_t Q2l;
  double W, Q2, CosTheta_isb, SinTheta_isb, Phi_isb;
  Range1D_t Wl  = kps.WLim_RSPP();
  if (fWcut >= Wl.min)
    Wl.max = TMath::Min(fWcut,Wl.max);
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

     
     W  = Wl.min + (Wl.max - Wl.min)*rnd->RndKine().Rndm();
     interaction->KinePtr()->SetW(W);
     Q2l = kps.Q2Lim_W_RSPP();
     Q2 = Q2l.min + (Q2l.max - Q2l.min)*rnd->RndKine().Rndm();
     CosTheta_isb = -1. + 2.*rnd->RndKine().Rndm();
     SinTheta_isb = (1 - CosTheta_isb*CosTheta_isb)<0?0:TMath::Sqrt(1 - CosTheta_isb*CosTheta_isb);
     Phi_isb = 2*kPi*rnd->RndKine().Rndm();
     
     xin[0] = W;
     xin[1] = Q2;
     xin[2] = CosTheta_isb;
     xin[3] = Phi_isb;
     
     
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
  
  fResList.Clear();
  string resonances ;
  GetParam( "ResonanceNameList", resonances ) ;
  fResList.DecodeFromNameList(resonances);
  
  // Safety factor for the maximum differential cross section
  this->GetParamDef("MaxXSec-SafetyFactor", fSafetyFactor, 1.25);
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
double RSPPEventGenerator::ComputeMaxXSec(
                                       const Interaction * in) const
{
   Interaction *interaction = const_cast<Interaction*>(in);
   KPhaseSpace * kps = interaction->PhaseSpacePtr();
   Range1D_t Wl = kps->WLim_RSPP();
   if (fWcut >= Wl.min)
     Wl.max = TMath::Min(fWcut,Wl.max);
   double dW = Wl.max - Wl.min;
   const InitialState & init_state = interaction -> InitState();
   ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
   ROOT::Math::IBaseFunctionMultiDim * f = new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fXSecModel, interaction);
   min->SetFunction( *f );
   min->SetMaxFunctionCalls(10000);  // for Minuit/Minuit2
   min->SetMaxIterations(10000);     // for GSL
   min->SetTolerance(0.001);
   min->SetPrintLevel(0);
   double min_xsec = 0.;
   double xsec;
   double step = 1e-7;
   double Enu = init_state.ProbeE(kRfHitNucRest);
   if (Enu <= 15)
   {
     // low-energy heuristic algorithm for maximum search
     int N3 = 2;
     int N4 = 4;
     for (auto res : fResList)
     {
       double MR  = utils::res::Mass(res);
       double WR  = utils::res::Width(res);
       double W = MR;
       double Wmin = TMath::Max(Wl.min, MR - WR);
       double Wmax = TMath::Min(Wl.max, MR + WR);
       interaction->KinePtr()->SetW(W);
       Range1D_t Q2l = kps->Q2Lim_W_RSPP();
       double dQ2 = Q2l.max-Q2l.min;
       double Q2 = Q2l.min + dQ2/6.;
       double Q2min = Q2l.min;
       double Q2max = Q2l.min + dQ2/(Enu < 1?1.:3.);
       for (int i3 = 0; i3 < N3; i3++)
       {
         double cost    = i3 - 1;
         double costmin = i3 - 1;
         double costmax = i3;
         for (int i4 = 0; i4 < N4; i4++)
         {
           double phi = 2*kPi*i4/N4;
           double phimin = 2*kPi*i4/N4;
           double phimax = 2*kPi*(i4 + 1)/N4;
           
           min->SetVariable(0, "x1", W, (Wmax - Wmin)*step);
           min->SetVariable(1, "x2", Q2, (Q2max - Q2min)*step);
           min->SetVariable(2, "x3", cost, (costmax - costmin)*step);
           min->SetVariable(3, "x4", phi, (phimax - phimin)*step);
           min->SetVariableLimits(0, Wmin, Wmax);
           min->SetVariableLimits(1, Q2min, Q2max);
           min->SetVariableLimits(2, costmin, costmax);
           min->SetVariableLimits(3, phimin, phimax);
           min->Minimize();
           xsec = min->MinValue();
           if (xsec < min_xsec)
             min_xsec = xsec;
        }
      }
    }
   }
   else
   {
     // high-energy heuristic algorithm for maximum search
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
         
         double W  = Wl.min + dW*x1m;
         interaction->KinePtr()->SetW(W);
         Range1D_t Q2l = kps->Q2Lim_W_RSPP();
         double dQ2 = Q2l.max-Q2l.min;
         
         min->SetVariable(0, "x1", W, dW*TMath::Abs(cells[cell].Vertex2.x1-cells[cell].Vertex1.x1)*step);
         min->SetVariable(1, "x2", Q2l.min + dQ2*x2m, dQ2*TMath::Abs(cells[cell].Vertex2.x2-cells[cell].Vertex1.x2)*step);
         min->SetVariable(2, "x3", -1. + 2.*x3m, 2.*TMath::Abs(cells[cell].Vertex2.x3-cells[cell].Vertex1.x3)*step);
         min->SetVariable(3, "x4", 2*kPi*x4m, 2.*kPi*TMath::Abs(cells[cell].Vertex2.x4-cells[cell].Vertex1.x4)*step);
         min->SetVariableLimits(0, Wl.min + dW*cells[cell].Vertex1.x1, Wl.min + dW*cells[cell].Vertex2.x1);
         min->SetVariableLimits(1, Q2l.min + dQ2*cells[cell].Vertex1.x2, Q2l.min + dQ2*cells[cell].Vertex2.x2);
         min->SetVariableLimits(2, -1. + 2.*cells[cell].Vertex1.x3, -1. + 2.*cells[cell].Vertex2.x3);
         min->SetVariableLimits(3, 2*kPi*cells[cell].Vertex1.x4, 2*kPi*cells[cell].Vertex2.x4);
         min->Minimize();
         xsec = min->MinValue();
         if (xsec < min_xsec)
           min_xsec = xsec;
         const double *xs = min->X();
         Vertex minv((xs[0]-Wl.min)/dW, (xs[1]-Q2l.min)/dQ2, (xs[2]+1)/2., xs[3]/2./kPi);
         if (minv == cells[cell].Vertex1 || minv == cells[cell].Vertex2)
           minv = Vertex (x1m, x2m, x3m, x4m);
     
         if (dep < fMaxDepth - 1)
         {
     
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
     }
  }
  
  delete f;

  return -fSafetyFactor*min_xsec;
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
  
  fInteraction->KinePtr()->SetW(xin[0]); 
   
  fInteraction->KinePtr()->SetQ2(xin[1]);
  
  fInteraction->KinePtr()->SetKV(kKVctp,  xin[2]); // cosine of pion theta in resonance rest frame
  
  fInteraction->KinePtr()->SetKV(kKVphip, xin[3]); // pion phi in resonance rest frame
    
  double xsec = -fModel->XSec(fInteraction, kPSWQ2ctpphipfE);
  return xsec/(1E-38 * units::cm2);
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::Clone() const
{
  return
    new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fModel,fInteraction);
}

