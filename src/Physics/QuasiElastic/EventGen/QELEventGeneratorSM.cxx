//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Igor Kakorin <kakorin@jinr.ru>
 Joint Institute for Nuclear Research

 adapted from fortran code provided by:

 Konstantin Kuzmin <kkuzmin@theor.jinr.ru>,
 Joint Institute for Nuclear Research,
 Institute for Theoretical and Experimental Physics

 Vadim Naumov <vnaumov@theor.jinr.ru>,
 Joint Institute for Nuclear Research

 based on code of:
 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <TMath.h>
#include <Math/Factory.h>
#include <Math/Minimizer.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Physics/QuasiElastic/EventGen/QELEventGeneratorSM.h"


#include "Framework/Utils/Range1.h"
#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;
using namespace genie::utils::gsl;

namespace { // anonymous namespace (file only visibility)
  const double eps = std::numeric_limits<double>::epsilon();
  const double max_dbl = std::numeric_limits<double>::max();
  const double min_dbl = std::numeric_limits<double>::min();
}
//___________________________________________________________________________
QELEventGeneratorSM::QELEventGeneratorSM() :
KineGeneratorWithCache("genie::QELEventGeneratorSM")
{

}
//___________________________________________________________________________
QELEventGeneratorSM::QELEventGeneratorSM(string config) :
KineGeneratorWithCache("genie::QELEventGeneratorSM", config)
{

}
//___________________________________________________________________________
QELEventGeneratorSM::~QELEventGeneratorSM()
{

}
//___________________________________________________________________________
void QELEventGeneratorSM::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("QELEvent", pINFO) << "Generating QE event kinematics...";

  if(fGenerateUniformly) {
    LOG("QELEvent", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }
  // Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction -> InitState();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // Skip if no hit nucleon is set
  if(! evrec->HitNucleon())
  {
    LOG("QELEvent", pFATAL) << "No hit nucleon was set";
    gAbortingInErr = true;
    exit(1);
  }

  // Access the target from the interaction summary
  Target * tgt = init_state.TgtPtr();
  GHepParticle * nucleon = evrec->HitNucleon();
  // Store position of nucleon
  double hitNucPos = nucleon->X4()->Vect().Mag();
  tgt->SetHitNucPosition( hitNucPos );

  // Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  // Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  // heavy nucleus is nucleus that heavier than tritium or 3He.
  bool isHeavyNucleus = tgt->A()>3;

  sm_utils->SetInteraction(interaction);
  // phase space for heavy nucleus is different from light one
  fkps = isHeavyNucleus?kPSQ2vpfE:kPSQ2fE;
  Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
  // Try to calculate the maximum cross-section in kinematical limits
  // if not pre-computed already
  double xsec_max1  = fGenerateUniformly ? -1 : this->MaxXSec(evrec);
  // this make correct calculation of probability
  double xsec_max2  = fGenerateUniformly ? -1 : (rQ2.max<fQ2Min)? 0:this->MaxXSec(evrec, 1);
  double dvmax= isHeavyNucleus ? this->MaxXSec(evrec, 2) : 0.;


  // generate Q2, v, pF
  double Q2, v, kF, xsec;
  unsigned int iter = 0;
  bool accept = false;
  TLorentzVector q;
  while(1)
  {
     LOG("QELEvent", pINFO) << "Attempt #: " << iter;
     if(iter > 100*kRjMaxIterations)
     {
        LOG("QELEvent", pWARN)
          << "Couldn't select a valid kinematics after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

      // Pick Q2, v and pF
     double xsec_max = 0.;
     double pth = rnd->RndKine().Rndm();
     //pth < prob1/(prob1+prob2), where prob1,prob2 - probabilities to generate event in area1 (Q2<fQ2Min) and area2 (Q2>fQ2Min) which are not normalized
     if (pth <= xsec_max1*(TMath::Min(rQ2.max, fQ2Min)-rQ2.min)/(xsec_max1*(TMath::Min(rQ2.max, fQ2Min)-rQ2.min)+xsec_max2*(rQ2.max-fQ2Min)))
     {
       xsec_max = xsec_max1;
       Q2 = (rnd->RndKine().Rndm() * (TMath::Min(rQ2.max, fQ2Min)-rQ2.min)) + rQ2.min;
     }
     else
     {
        xsec_max = xsec_max2;
        Q2 = (rnd->RndKine().Rndm() * (rQ2.max-fQ2Min)) + fQ2Min;
     }
     Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
     // for nuclei heavier than deuterium generate energy transfer in defined energy interval
     double dv = 0.;
     if (isHeavyNucleus)
     {
       dv = dvmax * rnd->RndKine().Rndm(); 
       if (dv > (rv.max-rv.min))
       {
          continue;
       }
     }
     v = rv.min + dv;
     
     Range1D_t rkF = sm_utils->kFQES_SM_lim(Q2, v);
     // rkF.max = Fermi momentum
     kF = rnd->RndKine().Rndm()*sm_utils->GetFermiMomentum();
     if (kF < rkF.min)
     {
        continue;
     }

     Kinematics * kinematics = interaction->KinePtr();
     kinematics->SetKV(kKVQ2, Q2);
     kinematics->SetKV(kKVv, v);
     kinematics->SetKV(kKVPn, kF);
     xsec = fXSecModel->XSec(interaction, fkps);
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
       interaction->ResetBit(kISkipProcessChk);
       interaction->ResetBit(kISkipKinematicChk);
       break;
     }
     iter++;
  }

  // Z-frame is frame where momentum transfer is (v,0,0,qv)
  double qv = TMath::Sqrt(v*v+Q2);
  TLorentzVector transferMom(0, 0, qv, v);

  // Momentum of initial neutrino in LAB frame
  TLorentzVector * tempTLV = evrec->Probe()->GetP4();
  TLorentzVector neutrinoMom = *tempTLV;
  delete tempTLV;

  // define all angles in Z frame
  double theta = neutrinoMom.Vect().Theta();
  double phi =   neutrinoMom.Vect().Phi();
  double theta_k = sm_utils-> GetTheta_k(v, qv);
  double costheta_k = TMath::Cos(theta_k);
  double sintheta_k = TMath::Sin(theta_k);
  double E_p; //energy of initial nucleon
  double theta_p = sm_utils-> GetTheta_p(kF, v, qv, E_p);
  double costheta_p = TMath::Cos(theta_p);
  double sintheta_p = TMath::Sin(theta_p);
  double fi_p       = 2 * TMath::Pi() * rnd->RndGen().Rndm();
  double cosfi_p    = TMath::Cos(fi_p);
  double sinfi_p    = TMath::Sin(fi_p);
  double psi       = 2 * TMath::Pi() * rnd->RndGen().Rndm();

  // 4-momentum of initial neutrino in Z-frame
  TLorentzVector neutrinoMomZ(neutrinoMom.P()*sintheta_k, 0, neutrinoMom.P()*costheta_k, neutrinoMom.E());

  // 4-momentum of final lepton in Z-frame
  TLorentzVector outLeptonMom = neutrinoMomZ-transferMom;

  // 4-momentum of initial nucleon in Z-frame
  TLorentzVector inNucleonMom(kF*sintheta_p*cosfi_p, kF*sintheta_p*sinfi_p, kF*costheta_p, E_p);

  // 4-momentum of final nucleon in Z-frame
  TLorentzVector outNucleonMom = inNucleonMom+transferMom;

  // Rotate all vectors from Z frame to LAB frame
  TVector3 yvec (0.0, 1.0, 0.0);
  TVector3 zvec (0.0, 0.0, 1.0);
  TVector3 unit_nudir = neutrinoMom.Vect().Unit();

  outLeptonMom.Rotate(theta-theta_k, yvec);
  outLeptonMom.Rotate(phi, zvec);

  inNucleonMom.Rotate(theta-theta_k, yvec);
  inNucleonMom.Rotate(phi, zvec);

  outNucleonMom.Rotate(theta-theta_k, yvec);
  outNucleonMom.Rotate(phi, zvec);

  outLeptonMom.Rotate(psi, unit_nudir);
  inNucleonMom.Rotate(psi, unit_nudir);
  outNucleonMom.Rotate(psi, unit_nudir);

  // set 4-momentum of struck nucleon
  TLorentzVector * p4 = tgt->HitNucP4Ptr();
  p4->SetPx( inNucleonMom.Px() );
  p4->SetPy( inNucleonMom.Py() );
  p4->SetPz( inNucleonMom.Pz() );
  p4->SetE ( inNucleonMom.E() );

  int rpdgc = interaction->RecoilNucleonPdg();
  assert(rpdgc);
  double W = PDGLibrary::Instance()->Find(rpdgc)->Mass();
  LOG("QELEvent", pNOTICE) << "Selected: W = "<< W;
  double M = init_state.Tgt().HitNucP4().M();
  double E  = init_state.ProbeE(kRfHitNucRest);

  // (W,Q2) -> (x,y)
  double x=0, y=0;
  kinematics::WQ2toXY(E,M,W,Q2,x,y);

  // lock selected kinematics & clear running values
  interaction->KinePtr()->SetQ2(Q2, true);
  interaction->KinePtr()->SetW (W,  true);
  interaction->KinePtr()->Setx (x,  true);
  interaction->KinePtr()->Sety (y,  true);
  interaction->KinePtr()->ClearRunningValues();

  // set the cross section for the selected kinematics
  evrec->SetDiffXSec(xsec,fkps);
  if(fGenerateUniformly)
  {
     double vol     = sm_utils->PhaseSpaceVolume(fkps);
     double totxsec = evrec->XSec();
     double wght    = (vol/totxsec)*xsec;
     LOG("QELEvent", pNOTICE)  << "Kinematics wght = "<< wght;

     // apply computed weight to the current event weight
     wght *= evrec->Weight();
     LOG("QELEvent", pNOTICE) << "Current event wght = " << wght;
     evrec->SetWeight(wght);
  }
  TLorentzVector x4l(*(evrec->Probe())->X4());

  // Add the final-state lepton to the event record
  evrec->AddParticle(interaction->FSPrimLeptonPdg(), kIStStableFinalState, evrec->ProbePosition(),-1,-1,-1, outLeptonMom, x4l);

  // Set its polarization
  utils::SetPrimaryLeptonPolarization( evrec );

  GHepStatus_t ist;
  if (!fGenerateNucleonInNucleus)
     ist = kIStStableFinalState;
  else
     ist = (tgt->IsNucleus() && isHeavyNucleus) ? kIStHadronInTheNucleus : kIStStableFinalState;

  GHepParticle outNucleon(interaction->RecoilNucleonPdg(), ist, evrec->HitNucleonPosition(),-1,-1,-1, outNucleonMom , x4l);
  evrec->AddParticle(outNucleon);

  // Store struck nucleon momentum
  LOG("QELEvent",pNOTICE) << "pn: " << inNucleonMom.X() << ", " <<inNucleonMom.Y() << ", " <<inNucleonMom.Z() << ", " <<inNucleonMom.E();
  nucleon->SetMomentum(inNucleonMom);
  nucleon->SetRemovalEnergy(sm_utils->GetBindingEnergy());

  // skip if not a nuclear target
  if(evrec->Summary()->InitState().Tgt().IsNucleus())
    // add a recoiled nucleus remnant
    this->AddTargetNucleusRemnant(evrec);

  LOG("QELEvent", pINFO) << "Done generating QE event kinematics!";
}
//___________________________________________________________________________
void QELEventGeneratorSM::AddTargetNucleusRemnant(GHepRecord * evrec) const
{
// add the remnant nuclear target at the GHEP record

  LOG("QELEvent", pINFO) << "Adding final state nucleus";

  double Px = 0;
  double Py = 0;
  double Pz = 0;
  double E  = 0;

  GHepParticle * nucleus = evrec->TargetNucleus();
  int A = nucleus->A();
  int Z = nucleus->Z();

  int fd = nucleus->FirstDaughter();
  int ld = nucleus->LastDaughter();

  for(int id = fd; id <= ld; id++)
  {

     // compute A,Z for final state nucleus & get its PDG code and its mass
     GHepParticle * particle = evrec->Particle(id);
     assert(particle);
     int  pdgc = particle->Pdg();
     bool is_p  = pdg::IsProton (pdgc);
     bool is_n  = pdg::IsNeutron(pdgc);

     if (is_p) Z--;
     if (is_p || is_n) A--;

     Px += particle->Px();
     Py += particle->Py();
     Pz += particle->Pz();
     E  += particle->E();

  }//daughters

  TParticlePDG * remn = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  remn = PDGLibrary::Instance()->Find(ipdgc);
  if(!remn)
  {
    LOG("HadronicVtx", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
    assert(remn);
  }

  double Mi = nucleus->Mass();
  Px *= -1;
  Py *= -1;
  Pz *= -1;
  E = Mi-E;

  // Add the nucleus to the event record
  LOG("QELEvent", pINFO)
     << "Adding nucleus [A = " << A << ", Z = " << Z
     << ", pdgc = " << ipdgc << "]";

  int imom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
       ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);

  LOG("QELEvent", pINFO) << "Done";
  LOG("QELEvent", pINFO) << *evrec;
}
//___________________________________________________________________________
void QELEventGeneratorSM::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGeneratorSM::Configure(string config)
{
  Algorithm::Configure(config);

  Registry r( "QELEventGeneratorSM_specific", false ) ;
  r.Set( "sm_utils_algo", RgAlg("genie::SmithMonizUtils","Default") ) ;

  Algorithm::Configure(r) ;

  this->LoadConfig();
}
//____________________________________________________________________________
void QELEventGeneratorSM::LoadConfig(void)
{

  // Safety factors for the maximum differential cross section 
  fNumOfSafetyFactors = GetParamVect("Safety-Factors", vSafetyFactors, false);
  
  // Interpolator types for the maximum differential cross section 
  fNumOfInterpolatorTypes = GetParamVect("Interpolator-Types", vInterpolatorTypes, false);
  

  // Minimum energy for which max xsec would be cached, forcing explicit
  // calculation for lower eneries
  GetParamDef( "Cache-MinEnergy", fEMin, 1.00) ;
  GetParamDef( "Threshold-Q2", fQ2Min, 2.00);

  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false);

  // Generate nucleon in nucleus?
  GetParamDef( "IsNucleonInNucleus", fGenerateNucleonInNucleus, true);


  sm_utils = const_cast<genie::SmithMonizUtils *>(dynamic_cast<const genie::SmithMonizUtils *>( this -> SubAlg("sm_utils_algo") ) ) ;
}
//____________________________________________________________________________
double QELEventGeneratorSM::ComputeMaxXSec(const Interaction * interaction) const
{
    double xsec_max = -1;
    double tmp_xsec_max = -1;
    double Q20, v0;
    const int N_Q2 = 32;
    const InitialState & init_state = interaction -> InitState();
    Target * tgt = init_state.TgtPtr();
    bool isHeavyNucleus = tgt->A()>3;
    int N_v = isHeavyNucleus?32:0;

    Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
    const double logQ2min = TMath::Log(TMath::Max(rQ2.min, eps));
    const double logQ2max = TMath::Log(TMath::Min(rQ2.max, fQ2Min));
    Kinematics * kinematics = interaction->KinePtr();
    const double pFmax = sm_utils->GetFermiMomentum();
    // Now scan through kinematical variables Q2,v,kF
    for (int Q2_n=0; Q2_n <= N_Q2; Q2_n++)
    {
       // Scan around Q2
       double Q2 = TMath::Exp(Q2_n*(logQ2max-logQ2min)/N_Q2 + logQ2min);
       kinematics->SetKV(kKVQ2, Q2);
       Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
       const double logvmin = TMath::Log(TMath::Max(rv.min, eps));
       const double logvmax = TMath::Log(TMath::Max(rv.max, TMath::Max(rv.min, eps)));
       for (int v_n=0; v_n <= N_v; v_n++)
       {
          // Scan around v
          double v = TMath::Exp(v_n*(logvmax-logvmin)/N_v + logvmin);
          kinematics->SetKV(kKVv, v);
          kinematics->SetKV(kKVPn, pFmax);
          // Compute the QE cross section for the current kinematics
          double xs = fXSecModel->XSec(interaction, fkps);
          if (xs > tmp_xsec_max)
          {
             tmp_xsec_max = xs;
             Q20 = Q2;
             v0 =  v;
          }
       } // Done with v scan
    }// Done with Q2 scan
    
    const double Q2min = rQ2.min;
    const double Q2max = TMath::Min(rQ2.max, fQ2Min);
    const double vmin = sm_utils->vQES_SM_min(Q2min, Q2max);
    const double vmax = sm_utils->vQES_SM_max(Q2min, Q2max);
    
    ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
    ROOT::Math::IBaseFunctionMultiDim * f = isHeavyNucleus?static_cast<ROOT::Math::IBaseFunctionMultiDim*>(new d3XSecSM_dQ2dvdkF_E(fXSecModel, interaction, pFmax)):
                                                           static_cast<ROOT::Math::IBaseFunctionMultiDim*>(new d1XSecSM_dQ2_E(fXSecModel, interaction));
    min->SetFunction( *f );
    min->SetMaxFunctionCalls(10000);  // for Minuit/Minuit2
    min->SetMaxIterations(10000);     // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(0);
    double step = 1e-7;
    min->SetVariable(0, "Q2", Q20, step);
    min->SetVariableLimits(0, Q2min, Q2max);
    if (isHeavyNucleus)
    {
       min->SetVariable(1, "v",  v0,  step);
       min->SetVariableLimits(1, vmin,  vmax);
    }
    min->Minimize();
    xsec_max = -min->MinValue();
    if (tmp_xsec_max > xsec_max)
    {
       xsec_max = tmp_xsec_max;
    }

    return xsec_max;

}
//___________________________________________________________________________
double QELEventGeneratorSM::ComputeMaxXSec(const Interaction * interaction, const int nkey) const
{
  switch (nkey){
  case 0:
     return this->ComputeMaxXSec(interaction);
  
  case 1:
  {
     Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
     if (rQ2.max<fQ2Min)
     {
        return -1.;
     }
     double xsec_max = -1;
     double tmp_xsec_max = -1;
     double Q20, v0;
     const int N_Q2 = 32;
     const InitialState & init_state = interaction -> InitState();
     Target * tgt = init_state.TgtPtr();
     bool isHeavyNucleus = tgt->A()>3;
     int N_v = isHeavyNucleus?32:0;
     
     const double logQ2min = TMath::Log(fQ2Min);
     const double logQ2max = TMath::Log(rQ2.max);
     Kinematics * kinematics = interaction->KinePtr();
     const double pFmax = sm_utils->GetFermiMomentum();
     // Now scan through kinematical variables Q2,v,kF
     for (int Q2_n=0; Q2_n <= N_Q2; Q2_n++)
     {
        // Scan around Q2
        double Q2 = TMath::Exp(Q2_n*(logQ2max-logQ2min)/N_Q2 + logQ2min);
        kinematics->SetKV(kKVQ2, Q2);
        Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
        const double logvmin = TMath::Log(TMath::Max(rv.min, eps));
        const double logvmax = TMath::Log(TMath::Max(rv.max, TMath::Max(rv.min, eps)));
        for (int v_n=0; v_n <= N_v; v_n++)
        {
           // Scan around v
           double v = TMath::Exp(v_n*(logvmax-logvmin)/N_v + logvmin);
           kinematics->SetKV(kKVv, v);
           kinematics->SetKV(kKVPn, pFmax);
           // Compute the QE cross section for the current kinematics
           double xs = fXSecModel->XSec(interaction, fkps);
           if (xs > tmp_xsec_max)
           {
              tmp_xsec_max = xs;
              Q20 = Q2;
              v0 =  v;
           }
        } // Done with v scan
     }// Done with Q2 scan
     
     const double Q2min = fQ2Min;
     const double Q2max = rQ2.max;
     const double vmin = sm_utils->vQES_SM_min(Q2min, Q2max);
     const double vmax = sm_utils->vQES_SM_max(Q2min, Q2max);
     ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
     ROOT::Math::IBaseFunctionMultiDim * f = isHeavyNucleus?static_cast<ROOT::Math::IBaseFunctionMultiDim*>(new d3XSecSM_dQ2dvdkF_E(fXSecModel, interaction, pFmax)):
                                                            static_cast<ROOT::Math::IBaseFunctionMultiDim*>(new d1XSecSM_dQ2_E(fXSecModel, interaction));
     min->SetFunction( *f );
     min->SetMaxFunctionCalls(10000);  // for Minuit/Minuit2
     min->SetMaxIterations(10000);     // for GSL
     min->SetTolerance(0.001);
     min->SetPrintLevel(0);
     double step = 1e-7;
     min->SetVariable(0, "Q2", Q20, step);
     min->SetVariableLimits(0, Q2min, Q2max);
     if (isHeavyNucleus)
     {
        min->SetVariable(1, "v",  v0,  step);
        min->SetVariableLimits(1, vmin,  vmax);
     }
     min->Minimize();
     xsec_max = -min->MinValue();
     if (tmp_xsec_max > xsec_max)
     {
        xsec_max = tmp_xsec_max;
     }
     return xsec_max;
  }
  case 2:
  {
     double diffv_max = -1;
     double tmp_diffv_max = -1;
     const int N_Q2 = 100;
     double Q20;
     Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
     for (int Q2_n = 0; Q2_n<=N_Q2; Q2_n++) // Scan around Q2
     {
        double Q2 = rQ2.min + 1.*Q2_n * (rQ2.max-rQ2.min)/N_Q2;
        Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
        if (rv.max-rv.min > tmp_diffv_max)
        {
           tmp_diffv_max = rv.max-rv.min;
           Q20 = Q2;
        }
     } // Done with Q2 scan
     
     ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
     ROOT::Math::IBaseFunctionMultiDim * f = new dv_dQ2_E(interaction);
     min->SetFunction( *f );
     min->SetMaxFunctionCalls(10000);  // for Minuit/Minuit2
     min->SetMaxIterations(10000);     // for GSL
     min->SetTolerance(0.001);
     min->SetPrintLevel(0);
     double step = 1e-7;
     min->SetVariable(0, "Q2", Q20, step);
     min->SetVariableLimits(0, rQ2.min, rQ2.max);
     min->Minimize();
     diffv_max = -min->MinValue();
       
     if (tmp_diffv_max > diffv_max)
     {
        diffv_max = tmp_diffv_max;
     }
     return diffv_max;
  }
  default:
     return -1.;
  }
}
//___________________________________________________________________________
// GSL wrappers
//____________________________________________________________________________
d3XSecSM_dQ2dvdkF_E::d3XSecSM_dQ2dvdkF_E(
                                        const XSecAlgorithmI * m, 
                                        const Interaction    * i,
                                        double pF) :  ROOT::Math::IBaseFunctionMultiDim(), 
                                                      fModel(m), 
                                                      fInteraction(i),
                                                      fpF(pF)
{
}
d3XSecSM_dQ2dvdkF_E::~d3XSecSM_dQ2dvdkF_E()
{
}
unsigned int d3XSecSM_dQ2dvdkF_E::NDim(void) const
{
  return 2;
}
double d3XSecSM_dQ2dvdkF_E::DoEval(const double * xin) const
{
// outputs:
//   differential cross section
//
  fInteraction->KinePtr()->SetKV(kKVQ2, xin[0]);
  fInteraction->KinePtr()->SetKV(kKVv,  xin[1]);
  fInteraction->KinePtr()->SetKV(kKVPn, fpF);
  double xsec = -fModel->XSec(fInteraction, kPSQ2vpfE);
  return xsec;
}
ROOT::Math::IBaseFunctionMultiDim *
   d3XSecSM_dQ2dvdkF_E::Clone() const
{
  return new d3XSecSM_dQ2dvdkF_E(fModel, fInteraction, fpF);
}
//____________________________________________________________________________
d1XSecSM_dQ2_E::d1XSecSM_dQ2_E(
     const XSecAlgorithmI * m, 
     const Interaction    * i) :  ROOT::Math::IBaseFunctionMultiDim(), 
                                  fModel(m), 
                                  fInteraction(i)
{
}
d1XSecSM_dQ2_E::~d1XSecSM_dQ2_E()
{
}
unsigned int d1XSecSM_dQ2_E::NDim(void) const
{
  return 1;
}
double d1XSecSM_dQ2_E::DoEval(const double * xin) const
{
// outputs:
//   differential cross section
//
  fInteraction->KinePtr()->SetKV(kKVQ2, xin[0]);
  double xsec = -fModel->XSec(fInteraction, kPSQ2fE);
  return xsec;
}
ROOT::Math::IBaseFunctionMultiDim *
   d1XSecSM_dQ2_E::Clone() const
{
  return new d1XSecSM_dQ2_E(fModel, fInteraction);
}
//____________________________________________________________________________
dv_dQ2_E::dv_dQ2_E(const Interaction * i) : ROOT::Math::IBaseFunctionMultiDim(), 
                                                               fInteraction(i)
{
   AlgFactory * algf = AlgFactory::Instance();
   sm_utils = const_cast<SmithMonizUtils *>(dynamic_cast<const SmithMonizUtils *>(algf->GetAlgorithm("genie::SmithMonizUtils","Default")));
   sm_utils->SetInteraction(fInteraction);
}
dv_dQ2_E::~dv_dQ2_E()
{
}
unsigned int dv_dQ2_E::NDim(void) const
{
  return 1;
}
double dv_dQ2_E::DoEval(const double * xin) const
{
// outputs:
//   differential cross section
//
  double Q2  = xin[0];
  Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
  return rv.min-rv.max;
}
ROOT::Math::IBaseFunctionMultiDim *
   dv_dQ2_E::Clone() const
{
  return new dv_dQ2_E(fInteraction);
}
//____________________________________________________________________________
