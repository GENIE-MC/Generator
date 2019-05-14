//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Igor Kakorin <kakorin@jinr.ru>, Joint Institute for Nuclear Research
          adopted from  fortran code provided by
          Konstantin Kuzmin <kkuzmin@theor.jinr.ru>,
          Joint Institute for Nuclear Research,  Institute for Theoretical and Experimental Physics
          Vladimir Lyubushkin,
          Joint Institute for Nuclear Research
          Vadim Naumov <vnaumov@theor.jinr.ru>,
          Joint Institute for Nuclear Research
          based on code of Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.
*/
//____________________________________________________________________________

#include <TMath.h>

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
#include "Physics/QuasiElastic/XSection/SmithMonizUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Cache.h"
#include "Framework/Utils/CacheBranchFx.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

namespace { // anonymous namespace (file only visibility)
  const double eps = std::numeric_limits<double>::epsilon();
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

  // heavy nucleus is nucleus that heavier than hydrogen and deuterium
  bool isHeavyNucleus = tgt->A()>=3;

  sm_utils->SetInteraction(interaction);
  // phase space for heavy nucleus is different from light one
  fkps = isHeavyNucleus?kPSQ2vfE:kPSQ2fE;
  Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
  // Try to calculate the maximum cross-section in kinematical limits
  // if not pre-computed already
  double xsec_max1  = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);
  double xsec_max2  = (fGenerateUniformly) ? -1 : (rQ2.max<fQ2Min)? 0: this->MaxXSec2(evrec);// this make correct calculation of probability
  double vmax= isHeavyNucleus?this->MaxDiffv(evrec) : 0.;


  // generate Q2, v
  double gQ2, v, xsec;
  unsigned int iter = 0;
  bool accept = false;
  TLorentzVector q;
  while(1)
  {
     LOG("QELEvent", pINFO) << "Attempt #: " << iter;
     if(iter > kRjMaxIterations)
     {
        LOG("QELEvent", pWARN)
          << "Couldn't select a valid kinematics after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

      // Pick Q2 and v
     double xsec_max = 0.;
     double pth = rnd->RndKine().Rndm();
     //pth < prob1/(prob1+prob2), where prob1,prob2 - probabilities to generate event in area1 (Q2<fQ2Min) and area2 (Q2>fQ2Min) which are not normalized
     if (pth <= xsec_max1*(TMath::Min(rQ2.max, fQ2Min)-rQ2.min)/(xsec_max1*(TMath::Min(rQ2.max, fQ2Min)-rQ2.min)+xsec_max2*(rQ2.max-fQ2Min)))
     {
       xsec_max = xsec_max1;
       gQ2 = (rnd->RndKine().Rndm() * (TMath::Min(rQ2.max, fQ2Min)-rQ2.min)) + rQ2.min;
     }
     else
     {
        xsec_max = xsec_max2;
        gQ2 = (rnd->RndKine().Rndm() * (rQ2.max-fQ2Min)) + fQ2Min;
     }
     Range1D_t rv  = sm_utils->vQES_SM_lim(gQ2);
     // for nuclei heavier than deuterium generate energy transfer in defined energy interval
     v = 0.;
     if (isHeavyNucleus)
     {
       v = vmax * rnd->RndKine().Rndm();
       if (v > (rv.max-rv.min))
       {
          continue;
       }
     }
     v += rv.min;

     Kinematics * kinematics = interaction->KinePtr();
     kinematics->SetKV(kKVQ2, gQ2);
     kinematics->SetKV(kKVv, v);
     xsec = fXSecModel->XSec(interaction, fkps);

      //-- Decide whether to accept the current kinematics
     if(!fGenerateUniformly)
     {
       this->AssertXSecLimits(interaction, xsec, xsec_max);
       double t = xsec_max * rnd->RndKine().Rndm();

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("QELEvent", pDEBUG)<< "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
#endif
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
  double qv = TMath::Sqrt(v*v+gQ2);
  TLorentzVector transferMom(0, 0, qv, v);

  Range1D_t rkF = sm_utils->kFQES_SM_lim(gQ2, v);
  double kF = (rnd->RndKine().Rndm() * (rkF.max-rkF.min)) + rkF.min;

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
  double gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();
  LOG("QELEvent", pNOTICE) << "Selected: W = "<< gW;
  double M = init_state.Tgt().HitNucP4().M();
  double E  = init_state.ProbeE(kRfHitNucRest);

  // (W,Q2) -> (x,y)
  double gx=0, gy=0;
  kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

  // lock selected kinematics & clear running values
  interaction->KinePtr()->SetQ2(gQ2, true);
  interaction->KinePtr()->SetW (gW,  true);
  interaction->KinePtr()->Setx (gx,  true);
  interaction->KinePtr()->Sety (gy,  true);
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

  evrec->AddParticle(interaction->FSPrimLeptonPdg(), kIStStableFinalState, evrec->ProbePosition(),-1,-1,-1, outLeptonMom, x4l);

  GHepStatus_t ist;
  if (!fGenerateNucleonInNucleus)
     ist = kIStStableFinalState;
  else
     ist = (tgt->IsNucleus()) ? kIStHadronInTheNucleus : kIStStableFinalState;

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

  // Safety factor for the maximum differential cross section
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor,   1.2) ;
  GetParamDef( "MaxDiffv-SafetyFactor",fSafetyFacor_nu, 1.2);

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
    const int N_Q2 = 8;
    const int N_v = 8;

    Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
    const double logQ2min = TMath::Log(TMath::Max(rQ2.min, eps));
    const double logQ2max = TMath::Log(TMath::Min(rQ2.max, fQ2Min));

    double tmp_xsec_max = -1;
    // Now scan through kinematical variables Q2,v
    for (int Q2_n=0; Q2_n < N_Q2; Q2_n++)
    {  
       // Scan around Q2
       double Q2 = TMath::Exp(Q2_n * (logQ2max-logQ2min)/N_Q2 + logQ2min);
       Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
       const double logvmin = TMath::Log(TMath::Max(rv.min, eps));
       const double logvmax = TMath::Log(TMath::Max(rv.max, TMath::Max(rv.min, eps)));
       for (int v_n=0; v_n < N_v; v_n++)
       {  
          // Scan around v
          double v = TMath::Exp(v_n * (logvmax-logvmin)/N_v + logvmin);
          Kinematics * kinematics = interaction->KinePtr();
          kinematics->SetKV(kKVQ2, Q2);
          kinematics->SetKV(kKVv, v);
          // Compute the QE cross section for the current kinematics
          double xs = fXSecModel->XSec(interaction, fkps);
          if (xs > tmp_xsec_max)
            tmp_xsec_max = xs;
       } // Done with v scan
    }// Done with Q2 scan

    xsec_max = tmp_xsec_max;
    // Apply safety factor, since value retrieved from the cache might
    // correspond to a slightly different value
    xsec_max *= fSafetyFactor;
    return xsec_max;

}
//___________________________________________________________________________
double QELEventGeneratorSM::ComputeMaxXSec2(const Interaction * interaction) const
{
    double xsec_max = -1;
    const int N_Q2 = 8;
    const int N_v = 8;

    Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
    if (rQ2.max<fQ2Min) return xsec_max;
    const double logQ2min = TMath::Log(fQ2Min);
    const double logQ2max = TMath::Log(rQ2.max);

    double tmp_xsec_max = -1;
    // Now scan through kinematical variables Q2,v
    for (int Q2_n=0; Q2_n < N_Q2; Q2_n++)
    {  
       // Scan around Q2
       double Q2 = TMath::Exp(Q2_n * (logQ2max-logQ2min)/N_Q2 + logQ2min);
       Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
       const double logvmin = TMath::Log(TMath::Max(rv.min, eps));
       const double logvmax = TMath::Log(TMath::Max(rv.max, TMath::Max(rv.min, eps)));
       for (int v_n=0; v_n < N_v; v_n++)
       {  
          // Scan around v
          double v = TMath::Exp(v_n * (logvmax-logvmin)/N_v + logvmin);
          Kinematics * kinematics = interaction->KinePtr();
          kinematics->SetKV(kKVQ2, Q2);
          kinematics->SetKV(kKVv, v);
          // Compute the QE cross section for the current kinematics
          double xs = fXSecModel->XSec(interaction, fkps);
          if (xs > tmp_xsec_max)
            tmp_xsec_max = xs;
       } // Done with v scan
    }// Done with Q2 scan

    xsec_max = tmp_xsec_max;
    // Apply safety factor, since value retrieved from the cache might
    // correspond to a slightly different value
    xsec_max *= fSafetyFactor;
    return xsec_max;

}
//___________________________________________________________________________
double QELEventGeneratorSM::MaxXSec2(GHepRecord * event_rec) const
{
  LOG("Kinematics", pINFO)
                << "Getting max. differential xsec for the rejection method";

  double xsec_max = -1;
  Interaction * interaction = event_rec->Summary();

  LOG("Kinematics", pINFO)
                  << "Attempting to find a cached max{dxsec/dK} value";
  xsec_max = this->FindMaxXSec2(interaction);
  if(xsec_max>0) return xsec_max;

  LOG("Kinematics", pINFO)
                  << "Attempting to compute the max{dxsec/dK} value";
  xsec_max = this->ComputeMaxXSec2(interaction);
  if(xsec_max>0) {
     LOG("Kinematics", pINFO) << "max{dxsec/dK} = " << xsec_max;
     this->CacheMaxXSec2(interaction, xsec_max);
     return xsec_max;
  }

  LOG("Kinematics", pNOTICE)
            << "Can not generate event kinematics {K} (max_xsec({K};E)<=0)";
  // xsec for selected kinematics = 0
  event_rec->SetDiffXSec(0,kPSNull);
  // switch on error flag
  event_rec->EventFlags()->SetBitNumber(kKineGenErr, true);
  // reset 'trust' bits
  interaction->ResetBit(kISkipProcessChk);
  interaction->ResetBit(kISkipKinematicChk);
  // throw exception
  genie::exceptions::EVGThreadException exception;
  exception.SetReason("kinematics generation: max_xsec({K};E)<=0");
  exception.SwitchOnFastForward();
  throw exception;

  return 0;
}
//___________________________________________________________________________
double QELEventGeneratorSM::FindMaxXSec2(
                                       const Interaction * interaction) const
{
// Find a cached max xsec for the specified xsec algorithm & interaction and
// close to the specified energy

  // get neutrino energy
  double E = this->Energy(interaction);
  LOG("Kinematics", pINFO) << "E = " << E;

  if(E < fEMin) {
     LOG("Kinematics", pINFO)
         << "Below minimum energy - Forcing explicit calculation";
     return -1.;
  }

  // access the the cache branch
  CacheBranchFx * cb = this->AccessCacheBranch2(interaction);

  // if there are enough points stored in the cache buffer to build a
  // spline, then intepolate
  if( cb->Spl() ) {
     if( E >= cb->Spl()->XMin() && E <= cb->Spl()->XMax()) {
       double spl_max_xsec = cb->Spl()->Evaluate(E);
       LOG("Kinematics", pINFO)
          << "\nInterpolated: max xsec (E=" << E << ") = " << spl_max_xsec;
       return spl_max_xsec;
     }
     LOG("Kinematics", pINFO)
          << "Outside spline boundaries - Forcing explicit calculation";
     return -1.;
  }

  // if there are not enough points at the cache buffer to have a spline,
  // look whether there is another point that is sufficiently close
  double dE = TMath::Min(0.25, 0.05*E);
  const map<double,double> & fmap = cb->Map();
  map<double,double>::const_iterator iter = fmap.lower_bound(E);
  if(iter != fmap.end()) {
     if(TMath::Abs(E - iter->first) < dE) return iter->second;
  }

  return -1;

}
//___________________________________________________________________________
void QELEventGeneratorSM::CacheMaxXSec2(
                     const Interaction * interaction, double max_xsec) const
{
  LOG("Kinematics", pINFO)
                       << "Adding the computed max{dxsec/dK} value to cache";
  CacheBranchFx * cb = this->AccessCacheBranch2(interaction);

  double E = this->Energy(interaction);
  if(max_xsec>0) cb->AddValues(E,max_xsec);

  if(! cb->Spl() ) {
    if( cb->Map().size() > 40 ) cb->CreateSpline();
  }

  if( cb->Spl() ) {
     if( E < cb->Spl()->XMin() || E > cb->Spl()->XMax() ) {
        cb->CreateSpline();
     }
  }
}
//___________________________________________________________________________
CacheBranchFx * QELEventGeneratorSM::AccessCacheBranch2(
                                      const Interaction * interaction) const
{
// Returns the cache branch for this algorithm and this interaction. If no
// branch is found then one is created.

  Cache * cache = Cache::Instance();

  // build the cache branch key as: namespace::algorithm/config/interaction
  string algkey = this->Id().Key();
  string intkey = interaction->AsString();
  string key    = cache->CacheBranchKey(algkey, intkey, "2nd");

  CacheBranchFx * cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
  if(!cache_branch) {
    //-- create the cache branch at the first pass
    LOG("Kinematics", pINFO) << "No Max d^nXSec/d{K}^n cache branch found";
    LOG("Kinematics", pINFO) << "Creating cache branch - key = " << key;

    cache_branch = new CacheBranchFx("max[d^nXSec/d^n{K}] over phase space");
    cache->AddCacheBranch(key, cache_branch);
  }
  assert(cache_branch);

  return cache_branch;
}
//___________________________________________________________________________
double QELEventGeneratorSM::ComputeMaxDiffv(const Interaction *) const
{
  double max_diffv = -1;
  const int N_Q2 = 10;

  Range1D_t rQ2 = sm_utils->Q2QES_SM_lim();
  for (int Q2_n = 0; Q2_n<N_Q2; Q2_n++) // Scan around Q2
  {
     double Q2 = rQ2.min + 1.*Q2_n * (rQ2.max-rQ2.min)/N_Q2;
     Range1D_t rv  = sm_utils->vQES_SM_lim(Q2);
     if (rv.max-rv.min > max_diffv)
        max_diffv = rv.max-rv.min;
  } // Done with Q2 scan
  max_diffv *= fSafetyFactor;
  return max_diffv;

}
//___________________________________________________________________________
double QELEventGeneratorSM::MaxDiffv(GHepRecord * event_rec) const
{
  LOG("Kinematics", pINFO)
                << "Getting max. vmax(Q2)-vmin(Q2) for the rejection method";

  double max_diffv = -1;
  Interaction * interaction = event_rec->Summary();

  LOG("Kinematics", pINFO)
                  << "Attempting to find a cached max{vmax(Q2)-vmin(Q2)} value";
  max_diffv = this->FindMaxDiffv(interaction);
  if(max_diffv>0) return max_diffv;

  LOG("Kinematics", pINFO)
                  << "Attempting to compute the max{vmax(Q2)-vmin(Q2)} value";
  max_diffv = this->ComputeMaxDiffv(interaction);
  if(max_diffv>0) {
     LOG("Kinematics", pINFO) << "max{vmax(Q2)-vmin(Q2)} = " << max_diffv;
     this->CacheMaxDiffv(interaction, max_diffv);
     return max_diffv;
  }

  LOG("Kinematics", pNOTICE)
            << "Can not generate event kinematics (max{vmax(Q2)-vmin(Q2);E}<=0)";
  // xsec for selected kinematics = 0
  event_rec->SetDiffXSec(0,kPSNull);
  // switch on error flag
  event_rec->EventFlags()->SetBitNumber(kKineGenErr, true);
  // reset 'trust' bits
  interaction->ResetBit(kISkipProcessChk);
  interaction->ResetBit(kISkipKinematicChk);
  // throw exception
  genie::exceptions::EVGThreadException exception;
  exception.SetReason("kinematics generation: max{vmax(Q2)-vmin(Q2);E}<=0");
  exception.SwitchOnFastForward();
  throw exception;

  return 0;
}
//___________________________________________________________________________
double QELEventGeneratorSM::FindMaxDiffv(const Interaction * interaction) const
{
// Find a cached maximum of vmax(Q2)-vmin(Q2) for xsec algorithm & interaction and
// close to the specified energy

  // get neutrino energy
  double E = this->Energy(interaction);
  LOG("Kinematics", pINFO) << "E = " << E;

  if(E < fEMin) {
     LOG("Kinematics", pINFO)
         << "Below minimum energy - Forcing explicit calculation";
     return -1.;
  }

  // access the the cache branch
  CacheBranchFx * cb = this->AccessCacheBranchDiffv(interaction);

  // if there are enough points stored in the cache buffer to build a
  // spline, then intepolate
  if( cb->Spl() ) {
     if( E >= cb->Spl()->XMin() && E <= cb->Spl()->XMax()) 
     {
       double spl_maxdiffv = cb->Spl()->Evaluate(E);
       LOG("Kinematics", pINFO)
          << "\nInterpolated: max vmax(Q2)-vmin(Q2) (E=" << E << ") = " << spl_maxdiffv;
       return spl_maxdiffv;
     }
     LOG("Kinematics", pINFO)
          << "Outside spline boundaries - Forcing explicit calculation";
     return -1.;
  }

  // if there are not enough points at the cache buffer to have a spline,
  // look whether there is another point that is sufficiently close
  double dE = TMath::Min(0.25, 0.05*E);
  const map<double,double> & fmap = cb->Map();
  map<double,double>::const_iterator iter = fmap.lower_bound(E);
  if(iter != fmap.end()) 
  {
     if(TMath::Abs(E - iter->first) < dE) return iter->second;
  }

  return -1;

}
//___________________________________________________________________________
void QELEventGeneratorSM::CacheMaxDiffv(const Interaction * interaction, double max_diffv) const
{
  LOG("Kinematics", pINFO)
                       << "Adding the computed max{vmax(Q2)-vmin(Q2)} value to cache";
  CacheBranchFx * cb = this->AccessCacheBranchDiffv(interaction);

  double E = this->Energy(interaction);
  if(max_diffv>0) cb->AddValues(E,max_diffv);

  if(! cb->Spl() )
  {
    if( cb->Map().size() > 40 ) cb->CreateSpline();
  }

  if( cb->Spl() )
  {
     if( E < cb->Spl()->XMin() || E > cb->Spl()->XMax() ) 
     {
        cb->CreateSpline();
     }
  }
}
//___________________________________________________________________________
CacheBranchFx * QELEventGeneratorSM::AccessCacheBranchDiffv(const Interaction * interaction) const
{
// Returns the cache branch for this algorithm and this interaction. If no
// branch is found then one is created.

  Cache * cache = Cache::Instance();

  // build the cache branch key as: namespace::algorithm/config/interaction
  string algkey = this->Id().Key();
  string intkey = interaction->AsString();
  string key    = cache->CacheBranchKey(algkey, intkey, "diffv");

  CacheBranchFx * cache_branch =
              dynamic_cast<CacheBranchFx *> (cache->FindCacheBranch(key));
  if(!cache_branch)
  {
    //-- create the cache branch at the first pass
    LOG("Kinematics", pINFO) << "No Max vmax(Q2)-vmin(Q2) cache branch found";
    LOG("Kinematics", pINFO) << "Creating cache branch - key = " << key;

    cache_branch = new CacheBranchFx("max[vmax(Q2)-vmin(Q2)] over phase space");
    cache->AddCacheBranch(key, cache_branch);
  }
  assert(cache_branch);

  return cache_branch;
}
//___________________________________________________________________________
