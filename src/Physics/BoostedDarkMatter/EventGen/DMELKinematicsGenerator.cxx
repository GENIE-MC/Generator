//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 
*/
//____________________________________________________________________________

#include <TMath.h>

#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/BoostedDarkMatter/EventGen/DMELKinematicsGenerator.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::constants;
using namespace genie::utils;

//___________________________________________________________________________
DMELKinematicsGenerator::DMELKinematicsGenerator() :
KineGeneratorWithCache("genie::DMELKinematicsGenerator")
{

}
//___________________________________________________________________________
DMELKinematicsGenerator::DMELKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::DMELKinematicsGenerator", config)
{

}
//___________________________________________________________________________
DMELKinematicsGenerator::~DMELKinematicsGenerator()
{

}
//___________________________________________________________________________
void DMELKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("DMELKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // store the struck nucleon position for use by the xsec method
  double hitNucPos = evrec->HitNucleon()->X4()->Vect().Mag();
  interaction->InitStatePtr()->TgtPtr()->SetHitNucPosition(hitNucPos);

  //-- Note: The kinematic generator would be using the free nucleon cross
  //   section (even for nuclear targets) so as not to double-count nuclear
  //   suppression. This assumes that a) the nuclear suppression was turned
  //   on when computing the cross sections for selecting the current event 
  //   and that b) if the event turns out to be unphysical (Pauli-blocked) 
  //   the next attempted event will be forced to DM again.
  //   (discussion with Hugh - GENIE/NeuGEN integration workshop - 07APR2006
  interaction->SetBit(kIAssumeFreeNucleon);

  //-- Get the limits for the generated Q2
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t Q2 = kps.Limits(kKVQ2);

  if(Q2.max <=0 || Q2.min>=Q2.max) {
     LOG("DMELKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  LOG("DMELKinematics",pNOTICE) << "Setting max XSec";
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);
  LOG("DMELKinematics",pNOTICE) << "Set max XSec to " << xsec_max;

  //-- Try to select a valid Q2 using the rejection method

  // kinematical limits
  double Q2min  = Q2.min+kASmallNum;
  double Q2max  = Q2.max-kASmallNum;
//double QD2min = utils::kinematics::Q2toQD2(Q2min);
//double QD2max = utils::kinematics::Q2toQD2(Q2max);
  double xsec   = -1.;
  double gQ2    =  0.;

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("DMELKinematics", pWARN)
          << "Couldn't select a valid Q^2 after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     //-- Generate a Q2 value within the allowed phase space
/*
     if(fGenerateUniformly) {
         gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     } else {
         // In unweighted mode - use transform that takes out the dipole form
         double gQD2 = QD2min + (QD2max-QD2min) * rnd->RndKine().Rndm();
         gQ2  = utils::kinematics::QD2toQ2(gQD2);
     }
*/
     LOG("DMELKinematics",pNOTICE) << "Attempting to sample";
     gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     interaction->KinePtr()->SetQ2(gQ2);
     LOG("DMELKinematics", pINFO) << "Trying: Q^2 = " << gQ2;

     //-- Computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSQ2fE);

     //-- Decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();
     //double J = kinematics::Jacobian(interaction,kPSQ2fE,kPSQD2fE);
        double J = 1.;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("DMELKinematics", pDEBUG)
            << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
#endif
        accept = (t < J*xsec);
     } else {
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("DMELKinematics", pINFO) << "Selected: Q^2 = " << gQ2;

        // reset bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);
        interaction->ResetBit(kIAssumeFreeNucleon);

        // compute the rest of the kinematical variables

        // get neutrino energy at struck nucleon rest frame and the
        // struck nucleon mass (can be off the mass shell)
        const InitialState & init_state = interaction->InitState();
        double E  = init_state.ProbeE(kRfHitNucRest);
        double M = init_state.Tgt().HitNucP4().M();

        LOG("DMELKinematics", pNOTICE) << "E = " << E << ", M = "<< M;

        // The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
        // For DM/Charm events it is set to be equal to the on-shell mass of
        // the generated charm baryon (Lamda_c+, Sigma_c+ or Sigma_c++)
        // Similarly for strange baryons
        //
        const XclsTag & xcls = interaction->ExclTag();
        int rpdgc = 0;
        if(xcls.IsCharmEvent()) { rpdgc = xcls.CharmHadronPdg();           }
        else if(xcls.IsStrangeEvent()) { rpdgc = xcls.StrangeHadronPdg();           }
        else                    { rpdgc = interaction->RecoilNucleonPdg(); }
        assert(rpdgc);
        double gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();

        LOG("DMELKinematics", pNOTICE) << "Selected: W = "<< gW;

        // (W,Q2) -> (x,y)
        double gx=0, gy=0;
        kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSQ2fE);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
          double vol     = kinematics::PhaseSpaceVolume(interaction,kPSQ2fE);
          double totxsec = evrec->XSec();
          double wght    = (vol/totxsec)*xsec;
          LOG("DMELKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
          wght *= evrec->Weight();
          LOG("DMELKinematics", pNOTICE) << "Current event wght = " << wght;
          evrec->SetWeight(wght);
        }

        // lock selected kinematics & clear running values
        interaction->KinePtr()->SetQ2(gQ2, true);
        interaction->KinePtr()->SetW (gW,  true);
        interaction->KinePtr()->Setx (gx,  true);
        interaction->KinePtr()->Sety (gy,  true);
        interaction->KinePtr()->ClearRunningValues();

        return;
     }
  }// iterations
}
//___________________________________________________________________________
void DMELKinematicsGenerator::SpectralFuncExperimentalCode(
  GHepRecord * evrec) const
{
  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = new Interaction(*evrec->Summary());
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  // store the struck nucleon position for use by the xsec method
  double hitNucPos = evrec->HitNucleon()->X4()->Vect().Mag();
  interaction->InitStatePtr()->TgtPtr()->SetHitNucPosition(hitNucPos);

  //-- Note: The kinematic generator would be using the free nucleon cross
  //   section (even for nuclear targets) so as not to double-count nuclear
  //   suppression. This assumes that a) the nuclear suppression was turned
  //   on when computing the cross sections for selecting the current event 
  //   and that b) if the event turns out to be unphysical (Pauli-blocked) 
  //   the next attempted event will be forced to DM again.
  //   (discussion with Hugh - GENIE/NeuGEN integration workshop - 07APR2006
  interaction->SetBit(kIAssumeFreeNucleon);

  //-- Assume scattering off a nucleon on the mass shell (PWIA prescription)
  double Mn  = interaction->InitState().Tgt().HitNucMass(); // PDG mass, take it to be on-shell
  double pxn = interaction->InitState().Tgt().HitNucP4().Px();
  double pyn = interaction->InitState().Tgt().HitNucP4().Py();
  double pzn = interaction->InitState().Tgt().HitNucP4().Pz();
  double En  = interaction->InitState().Tgt().HitNucP4().Energy();
  double En0 = TMath::Sqrt(pxn*pxn + pyn*pyn + pzn*pzn + Mn*Mn);
  double Eb  = En0 - En;
  interaction->InitStatePtr()->TgtPtr()->HitNucP4Ptr()->SetE(En0);
  double ml  = interaction->FSPrimLepton()->Mass();
  
  //-- Get the limits for the generated Q2
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t Q2 = kps.Limits(kKVQ2);

  if(Q2.max <=0 || Q2.min>=Q2.max) {
     LOG("DMELKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
//  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);
  double xsec_max = this->MaxXSec(evrec);

  // get neutrino energy at struck nucleon rest frame and the
  // struck nucleon mass (can be off the mass shell)
  const InitialState & init_state = interaction->InitState();
  double E  = init_state.ProbeE(kRfHitNucRest);

  LOG("DMELKinematics", pNOTICE) << "E = " << E << ", M = "<< Mn;

  //-- Try to select a valid Q2 using the rejection method

  // kinematical limits
  double Q2min  = Q2.min+kASmallNum;
  double Q2max  = Q2.max-kASmallNum;
  double xsec   = -1.;
  double gQ2    =  0.;
  double gW     =  0.;
  double gx     =  0.;
  double gy     =  0.;

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("DMELKinematics", pWARN)
          << "Couldn't select a valid Q^2 after " << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     //-- Generate a Q2 value within the allowed phase space
     gQ2 = Q2min + (Q2max-Q2min) * rnd->RndKine().Rndm();
     LOG("DMELKinematics", pNOTICE) << "Trying: Q^2 = " << gQ2;

     // The hadronic inv. mass is equal to the recoil nucleon on-shell mass.
     // For DM/Charm events it is set to be equal to the on-shell mass of
     // the generated charm baryon (Lamda_c+, Sigma_c+ or Sigma_c++)
     //
     const XclsTag & xcls = interaction->ExclTag();
     int rpdgc = 0;
     if(xcls.IsCharmEvent()) { rpdgc = xcls.CharmHadronPdg();           }
     else                    { rpdgc = interaction->RecoilNucleonPdg(); }
     assert(rpdgc);
     gW = PDGLibrary::Instance()->Find(rpdgc)->Mass();

     // (W,Q2) -> (x,y)
     kinematics::WQ2toXY(E,Mn,gW,gQ2,gx,gy);

     LOG("DMELKinematics", pNOTICE) << "W = "<< gW;
     LOG("DMELKinematics", pNOTICE) << "x = "<< gx;
     LOG("DMELKinematics", pNOTICE) << "y = "<< gy;

     // v
     double gv  = gy * E;
     double gv2 = gv*gv;

     LOG("DMELKinematics", pNOTICE) << "v = "<< gv;

     // v -> v~
     double gvtilde  = gv + Mn - Eb - TMath::Sqrt(Mn*Mn+pxn*pxn+pyn*pyn+pzn*pzn);
     double gvtilde2 = gvtilde*gvtilde;

     LOG("DMELKinematics", pNOTICE) << "v~ = "<< gvtilde;

     // Q~^2
     double gQ2tilde = gQ2 - gv2 + gvtilde2;

     LOG("DMELKinematics", pNOTICE) << "Q~^2 = "<< gQ2tilde;

     // Set updated Q2
     interaction->KinePtr()->SetQ2(gQ2tilde);

     //-- Computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSQ2fE);

     //-- Decide whether to accept the current kinematics
//     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->RndKine().Rndm();
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("DMELKinematics", pDEBUG)
            << "xsec= " << xsec << ", Rnd= " << t;
#endif
        accept = (t < xsec);
//     } else {
//        accept = (xsec>0);
//     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("DMELKinematics", pNOTICE) << "Selected: Q^2 = " << gQ2;

        // reset bits
//        interaction->ResetBit(kISkipProcessChk);
//        interaction->ResetBit(kISkipKinematicChk);
//        interaction->ResetBit(kIAssumeFreeNucleon);

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSQ2fE);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
//        if(fGenerateUniformly) {
//          double vol     = kinematics::PhaseSpaceVolume(interaction,kPSQ2fE);
//          double totxsec = evrec->XSec();
//          double wght    = (vol/totxsec)*xsec;
//          LOG("DMELKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
//          wght *= evrec->Weight();
//          LOG("DMELKinematics", pNOTICE) << "Current event wght = " << wght;
//          evrec->SetWeight(wght);
//        }

        // lock selected kinematics & clear running values
//        interaction->KinePtr()->SetQ2(gQ2, true);
//        interaction->KinePtr()->SetW (gW,  true);
//        interaction->KinePtr()->Setx (gx,  true);
//        interaction->KinePtr()->Sety (gy,  true);
//        interaction->KinePtr()->ClearRunningValues();

        evrec->Summary()->KinePtr()->SetQ2(gQ2, true);
        evrec->Summary()->KinePtr()->SetW (gW,  true);
        evrec->Summary()->KinePtr()->Setx (gx,  true);
        evrec->Summary()->KinePtr()->Sety (gy,  true);
        evrec->Summary()->KinePtr()->ClearRunningValues();
	delete interaction;

        return;
     }
  }// iterations
}
//___________________________________________________________________________
void DMELKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMELKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DMELKinematicsGenerator::LoadConfig(void)
{
// Load sub-algorithms and config data to reduce the number of registry
// lookups

  //-- Safety factor for the maximum differential cross section
	GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor , 1.25 ) ;

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
	GetParamDef( "Cache-MinEnergy", fEMin, 1.00 ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
	GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

}
//____________________________________________________________________________
double DMELKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small dQ2 step.

  double max_xsec = 0.0;

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t rQ2 = kps.Limits(kKVQ2);
  LOG("DMELKinematics", pDEBUG) << "Range of Q^2: " << rQ2.min << " to " << rQ2.max;
  if(rQ2.min <=0 || rQ2.max <= rQ2.min) return 0.;

  const double logQ2min = TMath::Log(rQ2.min + kASmallNum);
  const double logQ2max = TMath::Log(rQ2.max - kASmallNum);

  const int N  = 15;
  const int Nb = 10;

  double dlogQ2   = (logQ2max - logQ2min) /(N-1);
  double xseclast = -1;
  bool   increasing;

  for(int i=0; i<N; i++) {
     double Q2 = TMath::Exp(logQ2min + i * dlogQ2);
     interaction->KinePtr()->SetQ2(Q2);
     double xsec = fXSecModel->XSec(interaction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DMELKinematics", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
     max_xsec = TMath::Max(xsec, max_xsec);
     increasing = xsec-xseclast>=0;
     xseclast   = xsec;

     // once the cross section stops increasing, I reduce the step size and
     // step backwards a little bit to handle cases that the max cross section
     // is grossly underestimated (very peaky distribution & large step)
     if(!increasing) {
       dlogQ2/=(Nb+1);
       for(int ib=0; ib<Nb; ib++) {
	 Q2 = TMath::Exp(TMath::Log(Q2) - dlogQ2);
         if(Q2 < rQ2.min) continue;
         interaction->KinePtr()->SetQ2(Q2);
         xsec = fXSecModel->XSec(interaction, kPSQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
         LOG("DMELKinematics", pDEBUG)  << "xsec(Q2= " << Q2 << ") = " << xsec;
#endif
         max_xsec = TMath::Max(xsec, max_xsec);
       }
       break;
     }
  }//Q^2

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("DMELKinematics", pDEBUG) << interaction->AsString();
  SLOG("DMELKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DMELKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;
#endif

  return max_xsec;
}
//___________________________________________________________________________
