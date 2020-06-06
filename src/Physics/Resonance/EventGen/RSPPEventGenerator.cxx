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

#include <sstream>

#include <TMath.h>
#include <TFile.h>
#include <TKey.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <TH1.h>


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

using std::ostringstream;

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
RSPPEventGenerator::RSPPEventGenerator() :
EventRecordVisitorI("genie::RSPPEventGenerator")
{
   OpenFile();
}
//___________________________________________________________________________
RSPPEventGenerator::RSPPEventGenerator(string config) :
EventRecordVisitorI("genie::RSPPEventGenerator", config)
{
   OpenFile();
}
//___________________________________________________________________________
RSPPEventGenerator::~RSPPEventGenerator()
{
  fFoam_file->Close();
  delete fFoam_file;
}
//___________________________________________________________________________
void RSPPEventGenerator::OpenFile(void)
{
  TDirectory *& cur_dir = TDirectory::CurrentDirectory();
  fFoam_file = new TFile("foam_grid.root","RECREATE","foam");
  if (cur_dir)
    TDirectory::Cd(cur_dir->GetPath());
}
//___________________________________________________________________________
void RSPPEventGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  
  TDirectory *& cur_dir = TDirectory::CurrentDirectory();
  fFoam_file->cd();
  
  
  LOG("RSPPEventGen", pINFO) << "Generating resonance single pion production event kinematics...";
  
  //-- Get the interaction from the GHEP record
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction -> InitState();
  const KPhaseSpace& kps = interaction->PhaseSpace();

  // Access the target from the interaction summary
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
      
  bool isDestroy = false;
  
  double xin[4];
  double xsec = -1;
  double wght = 1.;
  if(fGenerateUniformly)
  {
     LOG("RSPPEventGen", pNOTICE)
           << "Generating kinematics uniformly over the allowed phase space";
     
     // generate W, Q2, cos(theta) and phi by accept-reject method
     unsigned int iter = 0;
     bool accept = false;
     while(1) 
     {
        iter++;
        if(iter > kRjMaxIterations) {
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
        xsec = (*f)(xin);
        
        accept = (xsec>0);
        
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
  }
  else
  {
       Interaction * ref_int = new Interaction (*interaction);
       double Eth = kps.Threshold_RSPP();
       int num_grid = TMath::Floor(.5 + (Ev - Eth)/fNextGriddE);
       InitialState * ref_init_state = ref_int -> InitStatePtr();
       Target * ref_tgt = ref_init_state->TgtPtr();
       ref_init_state->SetProbeP4(TLorentzVector(0., 0., 0., Eth + num_grid*fNextGriddE));
       ref_tgt->SetHitNucP4(p1_HNRF);
       genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E * rho = new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fXSecModel, ref_int);
       string key = CreateKey(num_grid, ref_int);
       fFoam = (TFoam*) fFoam_file->Get(key.c_str());
       if (!fFoam)
       {
          fFoam  = new TFoam(key.c_str());   // Create Simulator
          isDestroy = true;
          fFoam->SetkDim(       rho->NDim());      // Mandatory!!!Dimension of the integration space. Must be redefined!
          fFoam->SetnCells(     fFOAMNCells);    // optional. No of allocated number of cells. default: 1000
          fFoam->SetnSampl(     fFOAMNSampl);    // optional. No. of MC events in the cell MC exploration. default: 200
          fFoam->SetnBin(       fFOAMNBin);      // optional. No. of bins in edge-histogram in cell exploration. default: 8
          fFoam->SetOptRej(     fFOAMOptRej);    // optional. if 0 then weighted; if 1 then weight=1 for MC events. default: 1
          fFoam->SetOptDrive(   fFOAMOptDrive);  // optional. Maximum weight reduction, =1 for variance reduction. default: 2
          fFoam->SetEvPerBin(   fFOAMEvPerBin);  // optional. Maximum number of the effective weight=1 events/bin, if 0 deactivates this option. default: 25
          fFoam->SetChat(       fFOAMChat);      // optional  =0,1,2 is the `‘chat level’' in the standard output. default: 1
          fFoam->SetMaxWtRej(   fFOAMMaxWtRej);  // optional  Maximum weight used to get weight=1 MC events. default: 1.1
          fFoam->SetPseRan(&rnd->RndKine());     // Set random number generator
          fFoam->SetRho(rho);
          fFoam->Initialize();                   // Initialize simulator
       }
       else
           fFoam->SetRho(rho);
       
       if (isDestroy)
         fFoam->Write(key.c_str());     // Writing Foam on the disk, TESTING PERSISTENCY!!!
       
       fFoam->MakeEvent();                           // generate MC event
       fFoam->GetMCvect(xin);
       wght = fFoam->GetMCwt();
       xsec = (*f)(xin);
       wght *= xsec/(*rho)(xin);
       if (isDestroy)
         delete fFoam;
       fFoam = 0;
       delete rho;
       delete ref_int;
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
    wght    = (vol/totxsec)*xsec;
  }

  wght *= evrec->Weight();
  evrec->SetWeight(wght);

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

  fFoam_file->Write();
  TDirectory::Cd(cur_dir->GetPath());
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
string RSPPEventGenerator::CreateKey(int num_grid, Interaction * inter) const
{
  // build the cache branch key as: namespace::algorithm/config/interaction

  const InitialState & init_state = inter -> InitState();
    
  // Access the target from the interaction summary
  const Target & tgt = init_state.Tgt();
  const ProcessInfo & proc_info = inter -> ProcInfo();
  const XclsTag & xcls_tag = inter -> ExclTag();
  ostringstream interaction;
  interaction << num_grid << "_";
  interaction << "nu:"  << init_state.ProbePdg() << "_";
  interaction << "N:" << tgt.HitNucPdg() << "_";
  interaction << "proc:" << proc_info.InteractionTypeAsString() 
              << "," << proc_info.ScatteringTypeAsString()  << "_";
  
  string xcls = xcls_tag.AsString();
  interaction << xcls;
  
  return interaction.str();

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
  //// Safety factor for the maximum differential cross section
  //this->GetParamDef("MaxXSec-SafetyFactor", fSafetyFactor, 1.25);

  // Interval between reference energies of neutrino 
  // (in hit nucleon rest frame) at which the foam grid is built
  this->GetParamDef("EnergyWidthOfNextGrid", fNextGriddE, 0.5);
  if (fNextGriddE <= 0)
    fNextGriddE = 0.5;
  // No of allocated number of cells  
  this->GetParamDef("FOAM_NCells",   fFOAMNCells, 1000);
  // No. of MC events in the cell MC exploration  
  this->GetParamDef("FOAM_NSampl",   fFOAMNSampl,  200);
  // No. of bins in edge-histogram in cell exploration  
  this->GetParamDef("FOAM_NBin",     fFOAMNBin,      8);
  // OptRej = 0, weighted; OptRej=1, wt=1 MC events    
  this->GetParamDef("FOAM_OptRej",   fFOAMOptRej,    1);
  // Maximum weight reduction, =1 for variance reduction    
  this->GetParamDef("FOAM_OptDrive", fFOAMOptDrive,  2);
  // Maximum number of the effective wt=1 events/bin, EvPerBin=0 deactivates this option
  this->GetParamDef("FOAM_EvPerBin", fFOAMEvPerBin, 25);
  // =0,1,2 is the `‘chat level’' in the standard output   
  this->GetParamDef("FOAM_Chat",     fFOAMChat,      0);
  // Maximum weight used to get w=1 MC events   
  this->GetParamDef("FOAM_MaxWtRej", fFOAMMaxWtRej, 1.1);
    
  // Generate kinematics uniformly over allowed phase space and compute
  // an event weight?
  this->GetParamDef("UniformOverPhaseSpace", fGenerateUniformly, false);

}
//____________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::d4XSecMK_dWQ2CosThetaPhi_E(
     const XSecAlgorithmI * m, const Interaction * interaction) :
ROOT::Math::IBaseFunctionMultiDim(), TFoamIntegrand(), fModel(m)
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
    
  double xsec = fModel->XSec(fInteraction, kPSWQ2ctpphipfE)*(Wl.max-Wl.min)*(Q2l.max-Q2l.min)*4*kPi;
  return xsec/(1E-38 * units::cm2);
}
double genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::Density(int ndim, double * xin)
{
  return this->DoEval(xin); 
}
ROOT::Math::IBaseFunctionMultiDim *
   genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E::Clone() const
{
  return
    new genie::utils::gsl::d4XSecMK_dWQ2CosThetaPhi_E(fModel,fInteraction);
}

