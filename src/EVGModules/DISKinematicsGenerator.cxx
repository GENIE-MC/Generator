//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Controls.h"
#include "Conventions/KineVar.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/DISKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator() :
KineGeneratorWithCache("genie::DISKinematicsGenerator")
{

}
//___________________________________________________________________________
DISKinematicsGenerator::DISKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::DISKinematicsGenerator", config)
{

}
//___________________________________________________________________________
DISKinematicsGenerator::~DISKinematicsGenerator()
{

}
//___________________________________________________________________________
void DISKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  if(fGenerateUniformly) {
    LOG("DISKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Get the interaction and set the 'trust' bits
  Interaction * interaction = evrec->GetInteraction();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the physical W range taking into account any user cuts
  Range1D_t W  = this->WRange(interaction);

  if(W.max <=0 || W.min>=W.max) {
     LOG("DISKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kNoAvailablePhaseSpace, true);
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
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- get energy and 'nucleon mass' needed for W,Q2->x,y conversion
  const InitialState & init_state = interaction->GetInitialState();
  double Ev  = init_state.GetProbeE(kRfStruckNucAtRest);
  double M   = init_state.GetTarget().StruckNucleonP4()->M(); // can be off m-shell

  //-- Try to select a valid (x,y) pair using the rejection method

  double xmin  = kMinX;
  double xmax  = kMaxX;
  double ymin  = kMinY;
  double ymax  = kMaxY;
  double dx    = xmax - xmin;
  double dy    = ymax - ymin;

  double gx=-1, gy=-1, gW=-1, gQ2=-1, xsec=-1;

  register unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
       LOG("DISKinematics", pWARN)
        << " Couldn't select kinematics after " << iter << " iterations";

       evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select kinematics");
       exception.SwitchOnFastForward();
       throw exception;
     }

     //-- random x,y
     gx = xmin  + dx  * rnd->Random1().Rndm();
     gy = ymin  + dy  * rnd->Random1().Rndm();

     //-- (x,y) => (W,Q2)
     kinematics::XYtoWQ2(Ev,M,gW,gQ2,gx,gy);

     interaction->GetKinematicsPtr()->Setx(gx);
     interaction->GetKinematicsPtr()->Sety(gy);
     interaction->GetKinematicsPtr()->SetW(gW);
     interaction->GetKinematicsPtr()->SetQ2(gQ2);

     LOG("DISKinematics", pINFO) 
             << "Trying: x=" << gx << ", y=" << gy 
                             << " (-> W=" << gW << ", Q2=" << gQ2 << ")";

     //-- check whether we are in the allowed phase space
     Range1D_t Q2 = this->Q2Range(interaction);
     bool inW  = math::IsWithinLimits(gW, W );
     bool inQ2 = math::IsWithinLimits(gQ2,Q2);
     bool in   = inW && inQ2;

     if(!in) {
       LOG("DISKinematics", pINFO) << "*** Not in allowed phase space ";
       continue;
     }

     //-- compute the cross section for current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);

        double t = xsec_max * rnd->Random1().Rndm();
//      double J = kinematics::Jacobian(interaction,kPSxyfE,kPSlogxlogyfE);
        double J = 1;

        LOG("DISKinematics", pDEBUG)
              << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
        accept = (t < J*xsec);
     } 
     else {
        accept = (xsec>0);
     }

     if(accept) {
         // -------------- KINEMATICAL SELECTION DONE ----------------         
         LOG("DISKinematics", pNOTICE) 
                << "Selected: x=" << gx << ", y=" << gy 
                             << " (-> W=" << gW << ", Q2=" << gQ2 << ")";

         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec);

         // for uniform kinematics, compute an event weight as
         // wght = (phase space volume)*(differential xsec)/(event total xsec)
         if(fGenerateUniformly) {
            double vol     = kinematics::PhaseSpaceVolume(interaction,kPSxyfE);
            double totxsec = evrec->GetXSec();
            double wght    = (vol/totxsec)*xsec;
            LOG("DISKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

            // apply computed weight to the current event weight
            wght *= evrec->GetWeight();
            LOG("DISKinematics", pNOTICE) << "Current event wght = " << wght;
            evrec->SetWeight(wght);
         }

         // lock selected kinematics & clear running values
         interaction->GetKinematicsPtr()->SetW (gW,  true);
         interaction->GetKinematicsPtr()->SetQ2(gQ2, true);
         interaction->GetKinematicsPtr()->Setx (gx,  true);
         interaction->GetKinematicsPtr()->Sety (gy,  true);
         interaction->GetKinematicsPtr()->ClearRunningValues();

         return;
     }
  } // iterations
}
//___________________________________________________________________________
void DISKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadSubAlg(void)
{
// Reads its configuration from its Registry and loads all the sub-algorithms
// needed
  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadConfigData(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Get the user kinematical limits on W
  fWmin = fConfig->GetDoubleDef("W-min", -999999);
  fWmax = fConfig->GetDoubleDef("W-max",  999999);

  //-- Get the user kinematical limits on Q2
  fQ2min = fConfig->GetDoubleDef("Q2-min", -999999);
  fQ2max = fConfig->GetDoubleDef("Q2-max",  999999);

  //-- Safety factor for the maximum differential cross section
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
  fEMin = fConfig->GetDoubleDef("min-energy-cached", 1.0);

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  fMaxXSecDiffTolerance = fConfig->GetDoubleDef("max-xsec-diff-tolerance",0.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("uniform-over-phase-space", false);
}
//____________________________________________________________________________
Range1D_t DISKinematicsGenerator::WRange(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t W = kinematics::KineRange(interaction,kKVW);
  LOG("DISKinematics", pDEBUG)
             << "Physical W range: " << "[" << W.min << ", " << W.max << "]";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if(fWmin>0 && fWmax>0) 
         kinematics::ApplyCutsToKineLimits(W, fWmin, fWmax);
  LOG("DISKinematics", pDEBUG)
    << "W range (including cuts): " << "[" << W.min << ", " << W.max << "]";

  return W;
}
//___________________________________________________________________________
Range1D_t DISKinematicsGenerator::Q2Range(
                                       const Interaction * interaction) const
{
  //-- Get the physically allowed kinematical region for this interaction
  Range1D_t Q2 = kinematics::KineRange(interaction, kKVQ2);
  LOG("DISKinematics", pDEBUG)
         << "Physical Q2 range: " << "[" << Q2.min << ", " << Q2.max << "]";

  //-- Define the W range: the user selection (if any) is not allowed to
  //   extend it to an unphysical region but is allowed to narrow it down.
  if(fQ2min>0 && fQ2max>0)
         kinematics::ApplyCutsToKineLimits(Q2, fQ2min, fQ2max);
  LOG("DISKinematics", pDEBUG)
     << "Q2 range (including cuts): " << "["<< Q2.min<< ", "<< Q2.max<< "]";

  return Q2;
}
//___________________________________________________________________________
double DISKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But this needs to be fast - do not use a very fine grid.

  LOG("DISKinematics", pDEBUG)<< "Computing max xsec in allowed phase space";

  double max_xsec = 0.0;

  const int NW   = 20;
  const int NQ2  = 15;
  const int NQ2b =  5;

  const InitialState & init_state = interaction -> GetInitialState();
  double Ev = init_state.GetProbeE(kRfStruckNucAtRest);
  double M  = init_state.GetTarget().StruckNucleonMass(); 

  //-- Get the physical W range taking into account any user cuts
  Range1D_t W = this->WRange(interaction);
  assert(W.min>0 && W.max>W.min);

  double Wmin = W.min + kASmallNum;
  double Wmax = W.max - kASmallNum;
  double dW   = (Wmax-Wmin)/(NW-1);
  double gx=-1, gy=-1;

  double xseclast_w = -1;
  bool increasing_w;

  for(int i=0; i<NW; i++) {
     double gW = Wmin + i*dW;
     interaction->GetKinematicsPtr()->SetW(gW);

     LOG("DISKinematics", pDEBUG) << "W = " << gW;

     //-- get allowed Q2 range taking into account any user cuts
     Range1D_t Q2 = this->Q2Range(interaction);
     if(Q2.min<0 || Q2.max<Q2.min) continue;
     double logQ2min = TMath::Log(Q2.min+kASmallNum);
     double logQ2max = TMath::Log(Q2.max-kASmallNum);
     double dlogQ2   = (logQ2max-logQ2min)/(NQ2-1);

     double xseclast_q2 = -1;
     bool increasing_q2;

     for(int j=0; j<NQ2; j++) {
        double gQ2 = TMath::Exp(logQ2min + j*dlogQ2);
        interaction->GetKinematicsPtr()->SetQ2(gQ2);

        //-- (W,Q2) => (x,y)
        kinematics::WQ2toXY(Ev,M,gW,gQ2,gx,gy);
        interaction->GetKinematicsPtr()->Setx(gx);
        interaction->GetKinematicsPtr()->Sety(gy);

        double xsec = fXSecModel->XSec(interaction, kPSxyfE);
        LOG("DISKinematics", pINFO) 
                << "xsec(W=" << gW << ",Q2=" << gQ2 
   	                  << ",x=" << gx << ",y=" << gy << ")=" << xsec;

        // update maximum xsec
        max_xsec = TMath::Max(xsec, max_xsec);

        increasing_q2 = xsec-xseclast_q2>=0;
        xseclast_q2   = xsec;

        // once the cross section stops increasing, I reduce the step size and
        // step backwards a little bit to handle cases that the max cross section
        // is grossly underestimated (very peaky distribution & large step)
        if(!increasing_q2) {
          LOG("DISKinematics", pDEBUG) 
           << "d2xsec/dxdy|W stopped increasing. Stepping back & exiting Q2 loop";
          dlogQ2/=(NQ2b+1);
          for(int ik=0; ik<NQ2b; ik++) {
   	     gQ2 = TMath::Exp(TMath::Log(gQ2) - dlogQ2);
             if(gQ2 < Q2.min) continue;
             interaction->GetKinematicsPtr()->SetQ2(gQ2);
             kinematics::WQ2toXY(Ev,M,gW,gQ2,gx,gy);
             interaction->GetKinematicsPtr()->Setx(gx);
             interaction->GetKinematicsPtr()->Sety(gy);

             xsec = fXSecModel->XSec(interaction, kPSxyfE);
             LOG("DISKinematics", pINFO) 
                << "xsec(W=" << gW << ",Q2=" << gQ2 
   	                      << ",x=" << gx << ",y=" << gy << ")=" << xsec;
	  }
          break;
        } // stepping back within last bin
     } // Q2
     increasing_w = max_xsec-xseclast_w>=0;
     xseclast_w   = max_xsec;
     if(!increasing_w) {
       LOG("DISKinematics", pDEBUG) 
                          << "d2xsec/dxdy stopped increasing. Exiting W loop";
       break;
     }
  }// W

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  //  max_xsec *= fSafetyFactor;
  max_xsec *= ( (Ev<2.5) ? 2. : fSafetyFactor);

  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________

