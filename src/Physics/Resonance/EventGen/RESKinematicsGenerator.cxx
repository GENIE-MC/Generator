//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 03, 2009 - CA
   Moved into the new RES package from its previous location (EVGModules)
 @ Sep 21, 2009 - CA
   In ComputeMaxXSec() protect against NW=1 when calculating the step size dW
   as pointed out by Robert.
 @ Oct 20, 2009 - CA
   For charged lepton scattering do a more brute force kinematical selection
   and do not use the importance sampling envelope used for v scattering
   (was very inefficient)
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.
 @ Jul 26, 2018 - IL (Afroditi Papadopoulou, Adi Ashkenazi - Massachusetts Institute of Technology)
   Included importance sampling envelop both for neutrino and electron scattering
*/
//____________________________________________________________________________

#include <TMath.h>
#include <TF2.h>
#include <TROOT.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen//RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Resonance/EventGen/RESKinematicsGenerator.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
RESKinematicsGenerator::RESKinematicsGenerator() :
KineGeneratorWithCache("genie::RESKinematicsGenerator")
{
  fEnvelope = 0;
}
//___________________________________________________________________________
RESKinematicsGenerator::RESKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::RESKinematicsGenerator", config)
{
  fEnvelope = 0;
}
//___________________________________________________________________________
RESKinematicsGenerator::~RESKinematicsGenerator()
{
  if(fEnvelope) delete fEnvelope;
}
//___________________________________________________________________________
void RESKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("RESKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction from the GHEP record
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);

  //-- EM or CC/NC process (different importance sampling methods)
  bool is_em = interaction->ProcInfo().IsEM();

  //-- Compute the W limits
  //  (the physically allowed W's, unless an external cut is imposed)
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t W = kps.Limits(kKVW);
  //costas says this is bad  if(W.max>1.7) W.max=1.7;
//assert(W.min>=0. && W.min<W.max);

  if(W.max <=0 || W.min>=W.max) {
     LOG("RESKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  const InitialState & init_state = interaction -> InitState();
  double E = init_state.ProbeE(kRfHitNucRest);
  //  double M = init_state.Tgt().HitNucP4().M();
  //  double ml  = interaction->FSPrimLepton()->Mass();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Try to select a valid W, Q2 pair using the rejection method
  double dW   = W.max - W.min;
  double xsec = -1;

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
         LOG("RESKinematics", pWARN)
              << "*** Could not select a valid (W,Q^2) pair after "
                                                    << iter << " iterations";
         evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
         genie::exceptions::EVGThreadException exception;
         exception.SetReason("Couldn't select kinematics");
         exception.SwitchOnFastForward();
         throw exception;
     }

     double gW   = 0; // current hadronic invariant mass
     double gQ2  = 0; // current momentum transfer
     double gQD2 = 0; // tranformed Q2 to take out dipole form

     if(fGenerateUniformly) {
       //-- Generate a W uniformly in the kinematically allowed range.
       //   For the generated W, compute the Q2 range and generate a value
       //   uniformly over that range
       gW  = W.min + dW  * rnd->RndKine().Rndm();
       Range1D_t Q2 = kps.Q2Lim_W();
       if(Q2.max<=0. || Q2.min>=Q2.max) continue;
       gQ2 = Q2.min + (Q2.max-Q2.min) * rnd->RndKine().Rndm();

       interaction->SetBit(kISkipKinematicChk);

     } else {


       // > neutrino scattering
       // Selecting unweighted event kinematics using an importance sampling
       // method. Q2 with be transformed to QD2 to take out the dipole form.
       // An importance sampling envelope will be constructed for W.
         // first pass, configure the sampling envelope
         if(iter==1) {
            LOG("RESKinematics", pINFO) << "Initializing the sampling envelope";
            if(!fEnvelope) {
               LOG("RESKinematics", pFATAL) << "Null sampling envelope!";
               exit(1);
            }
            interaction->KinePtr()->SetW(W.min);
            Range1D_t Q2 = kps.Q2Lim_W();
	    double Q2min  = -99.;
            if (is_em) { Q2min  = Q2.min + kASmallNum; }
            else { Q2min  = 0 + kASmallNum; }
            double Q2max  = Q2.max - kASmallNum;
            
	    // In unweighted mode - use transform that takes out the dipole form
            double QD2min = utils::kinematics::Q2toQD2(Q2max);
            double QD2max = utils::kinematics::Q2toQD2(Q2min);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
            LOG("RESKinematics", pDEBUG)
                <<  "Q^2: [" << Q2min  << ", " << Q2max  << "] => "
                << "QD^2: [" << QD2min << ", " << QD2max << "]";
#endif
            double mR, gR;
            if(!interaction->ExclTag().KnownResonance()) {
               mR = 1.2;
               gR = 0.6;
            } else {
               Resonance_t res = interaction->ExclTag().Resonance();
               mR = res::Mass(res);
               gR = (E>mR) ? 0.220 : 0.400;
            }
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
            LOG("RESKinematics", pDEBUG)
               <<  "(m,g) = (" << mR << ", " << gR
               << "), max(xsec,W) = (" << xsec_max << ", " << W.max << ")";
#endif
            fEnvelope->SetRange(QD2min,W.min,QD2max,W.max); // range
            fEnvelope->SetParameter(0,  mR);                // resonance mass
            fEnvelope->SetParameter(1,  gR);                // resonance width
            fEnvelope->SetParameter(2,  xsec_max);          // max differential xsec
            fEnvelope->SetParameter(3,  W.max);             // kinematically allowed Wmax
         }// first pass

         // Generate W,QD2 using the 2-D envelope as PDF
         fEnvelope->GetRandom2(gQD2,gW);

         // QD2 -> Q2
         gQ2 = utils::kinematics::QD2toQ2(gQD2);
     } // uniformly over phase space?

     LOG("RESKinematics", pINFO) << "Trying: W = " << gW << ", Q2 = " << gQ2;

     //-- Set kinematics for current trial
     interaction->KinePtr()->SetW(gW);
     interaction->KinePtr()->SetQ2(gQ2);

     //-- Computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSWQ2fE);

     //-- Decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        // > charged lepton scattering
        if(is_em) {
          this->AssertXSecLimits(interaction, xsec, xsec_max);
          double t  = xsec_max * rnd->RndKine().Rndm();
          accept = (t < xsec);
	  LOG("RESKinematics", pINFO) << "xsec = " << xsec << ", ran*max = " << t << ", accept= " << accept;
      }
        // > neutrino scattering (using importance sampling envelope)
        else {
          double max = fEnvelope->Eval(gQD2, gW);
          double t   = max * rnd->RndKine().Rndm();
          double J   = kinematics::Jacobian(interaction,kPSWQ2fE,kPSWQD2fE);

          this->AssertXSecLimits(interaction, xsec, max);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("RESKinematics", pDEBUG)
                     << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
#endif
          accept = (t < J*xsec);
        } // charged lepton or neutrino scattering?
     }
     else {
        accept = (xsec>0);
     } // uniformly over phase space

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("RESKinematics", pINFO)
                            << "Selected: W = " << gW << ", Q2 = " << gQ2;
        // reset 'trust' bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);

        // compute x,y for selected W,Q2
        // note: hit nucleon can be off the mass-shell
        double gx=-1, gy=-1;
        double M = init_state.Tgt().HitNucP4().M();
        //double M = init_state.Tgt().HitNucMass();
        kinematics::WQ2toXY(E,M,gW,gQ2,gx,gy);

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec,kPSWQ2fE);

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
          double vol     = kinematics::PhaseSpaceVolume(interaction,kPSWQ2fE);
          double totxsec = evrec->XSec();
          double wght    = (vol/totxsec)*xsec;
          LOG("RESKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
          wght *= evrec->Weight();
          LOG("RESKinematics", pNOTICE) << "Current event wght = " << wght;
          evrec->SetWeight(wght);
        }

        // lock selected kinematics & clear running values
        interaction->KinePtr()->SetQ2(gQ2, true);
        interaction->KinePtr()->SetW (gW,  true);
        interaction->KinePtr()->Setx (gx,  true);
        interaction->KinePtr()->Sety (gy,  true);
        interaction->KinePtr()->ClearRunningValues();

        return;
     } // accept
  } // iterations
}
//___________________________________________________________________________
void RESKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void RESKinematicsGenerator::LoadConfig(void)
{
  // Safety factor for the maximum differential cross section
  this->GetParamDef("MaxXSec-SafetyFactor", fSafetyFactor, 1.25);

  // Minimum energy for which max xsec would be cached, forcing explicit
  // calculation for lower eneries
  this->GetParamDef("Cache-MinEnergy", fEMin, 0.5);

  // Load Wcut used in DIS/RES join scheme
  this->GetParam("Wcut", fWcut);

  // Maximum allowed fractional cross section deviation from maxim cross
  // section used in rejection method
  this->GetParamDef("MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999.);
  assert(fMaxXSecDiffTolerance>=0);

  // Generate kinematics uniformly over allowed phase space and compute
  // an event weight?
  this->GetParamDef("UniformOverPhaseSpace", fGenerateUniformly, false);

  // Envelope employed when importance sampling is used
  // (initialize with dummy range)
  if(fEnvelope) delete fEnvelope;
  fEnvelope = new TF2("res-envelope",
        kinematics::RESImportanceSamplingEnvelope,0.01,1,0.01,1,4);
  // stop ROOT from deleting this object of its own volition
  gROOT->GetListOfFunctions()->Remove(fEnvelope);
}
//____________________________________________________________________________
double RESKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But this needs to be fast - do not use a very fine grid.

  double max_xsec = 0.;

  const InitialState & init_state = interaction -> InitState();
  double E = init_state.ProbeE(kRfHitNucRest);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("RESKinematics", pDEBUG) << "Scanning phase space for E= " << E;
#endif
  //bool scan1d = (E>1.0);
  bool scan1d = false;

  double md;
  if(!interaction->ExclTag().KnownResonance()) md=1.23;
  else {
     Resonance_t res = interaction->ExclTag().Resonance();
     md=res::Mass(res);
  }

  if(scan1d) {
    // ** 1-D Scan
    //
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
    LOG("RESKinematics", pDEBUG)
              << "Will search for max{xsec} along W(=MRes) = " << md;
#endif
    // Set W around the value where d^2xsec/dWdQ^2 peaks
    interaction->KinePtr()->SetW(md);

    const KPhaseSpace & kps = interaction->PhaseSpace();
    Range1D_t rQ2 = kps.Q2Lim_W();
    if( rQ2.max < kMinQ2Limit || rQ2.min <=0 ) return 0.;

    int    NQ2      = 25;
    int    NQ2b     = 5;
    double logQ2min = TMath::Log(rQ2.min+kASmallNum);
    double logQ2max = TMath::Log(rQ2.max-kASmallNum);
    double dlogQ2 = (logQ2max - logQ2min) /(NQ2-1);

    double xseclast   = -1;
    bool   increasing = true;

    for(int iq2=0; iq2<NQ2; iq2++) {
      double Q2 = TMath::Exp(logQ2min + iq2 * dlogQ2);
      interaction->KinePtr()->SetQ2(Q2);
      double xsec = fXSecModel->XSec(interaction, kPSWQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
      LOG("RESKinematics", pDEBUG)
                << "xsec(W= " << md << ", Q2= " << Q2 << ") = " << xsec;
#endif
      max_xsec = TMath::Max(xsec, max_xsec);
      increasing = xsec-xseclast>=0;
      xseclast=xsec;

      // once the cross section stops increasing, I reduce the step size and
      // step backwards a little bit to handle cases that the max cross section
      // is grossly underestimated (very peaky distribution & large step)
      if(!increasing) {
        dlogQ2/=NQ2b;
        for(int iq2b=0; iq2b<NQ2b; iq2b++) {
	  Q2 = TMath::Exp(TMath::Log(Q2) - dlogQ2);
          if(Q2 < rQ2.min) continue;
          interaction->KinePtr()->SetQ2(Q2);
          xsec = fXSecModel->XSec(interaction, kPSWQ2fE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("RESKinematics", pDEBUG)
                 << "xsec(W= " << md << ", Q2= " << Q2 << ") = " << xsec;
#endif
          max_xsec = TMath::Max(xsec, max_xsec);
        }
        break;
      }
    } // Q2

  } else {

    // ** 2-D Scan
    //
    const KPhaseSpace & kps = interaction->PhaseSpace();
    Range1D_t rW = kps.WLim();

    int    NW   = 20;
    double Wmin = rW.min + kASmallNum;
    double Wmax = rW.max - kASmallNum;

    Wmax = TMath::Min(Wmax,fWcut);

    Wmin = TMath::Max(Wmin, md-.3);
    Wmax = TMath::Min(Wmax, md+.3);

    if(Wmax-Wmin<0.05) { NW=1; Wmin=Wmax; }

    double dW = (NW>1) ? (Wmax-Wmin)/(NW-1) : 0.;

    for(int iw=0; iw<NW; iw++) {
      double W = Wmin + iw*dW;
      interaction->KinePtr()->SetW(W);

      int NQ2  = 25;
      int NQ2b =  4;

      Range1D_t rQ2 = kps.Q2Lim_W();

      if( rQ2.max < kMinQ2Limit || rQ2.min <=0 ) continue;
      if( rQ2.max-rQ2.min<0.02 ) {NQ2=5; NQ2b=3;}

      double logQ2min   = TMath::Log(rQ2.min+kASmallNum);
      double logQ2max   = TMath::Log(rQ2.max-kASmallNum);
      double dlogQ2     = (logQ2max - logQ2min) /(NQ2-1);
      double xseclast   = -1;
      bool   increasing = true;

      for(int iq2=0; iq2<NQ2; iq2++) {
        double Q2 = TMath::Exp(logQ2min + iq2 * dlogQ2);
        interaction->KinePtr()->SetQ2(Q2);
        double xsec = fXSecModel->XSec(interaction, kPSWQ2fE);
        LOG("RESKinematics", pDEBUG)
                << "xsec(W= " << W << ", Q2= " << Q2 << ") = " << xsec;
        max_xsec = TMath::Max(xsec, max_xsec);
        increasing = xsec-xseclast>=0;
        xseclast=xsec;

        // once the cross section stops increasing, I reduce the step size and
        // step backwards a little bit to handle cases that the max cross section
        // is grossly underestimated (very peaky distribution & large step)
        if(!increasing) {
         dlogQ2/=NQ2b;
         for(int iq2b=0; iq2b<NQ2b; iq2b++) {
	   Q2 = TMath::Exp(TMath::Log(Q2) - dlogQ2);
           if(Q2 < rQ2.min) continue;
           interaction->KinePtr()->SetQ2(Q2);
           xsec = fXSecModel->XSec(interaction, kPSWQ2fE);
           LOG("RESKinematics", pDEBUG)
                 << "xsec(W= " << W << ", Q2= " << Q2 << ") = " << xsec;
           max_xsec = TMath::Max(xsec, max_xsec);
         }
         break;
	}
      } // Q2
    }//W
  }//2d scan

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  // Apply larger safety factor for smaller energies.
  max_xsec *= ( (E<md) ? 2. : fSafetyFactor);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("RESKinematics", pDEBUG) << interaction->AsString();
  LOG("RESKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  LOG("RESKinematics", pDEBUG) << "Computed using " << fXSecModel->Id();
#endif

  return max_xsec;
}
//___________________________________________________________________________
