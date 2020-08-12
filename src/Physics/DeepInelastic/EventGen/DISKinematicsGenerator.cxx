//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <cfloat>

#include <TMath.h>

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Physics/DeepInelastic/EventGen/DISKinematicsGenerator.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Physics/Common/RadiativeCorrector.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;
//using namespace genie::constants;
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
  if(fGenerateUniformly) {
    LOG("DISKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);

  //-- Get neutrino energy and hit 'nucleon mass'
  const InitialState & init_state = interaction->InitState();
  double Ev  = init_state.ProbeE(kRfHitNucRest);
  double M   = init_state.Tgt().HitNucP4().M(); // can be off m-shell

  //-- Get the physical W range
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t W  = kps.Limits(kKVW);
  if(W.max <=0 || W.min>=W.max) {
     LOG("DISKinematics", pWARN) << "No available phase space";
     evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("No available phase space");
     exception.SwitchOnFastForward();
     throw exception;
  }

  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  LOG("DISKinematics", pNOTICE) << "x: [" << xl.min << ", " << xl.max << "]";
  LOG("DISKinematics", pNOTICE) << "y: [" << yl.min << ", " << yl.max << "]";

  assert(xl.min>0 && yl.min>0);

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Try to select a valid (x,y) pair using the rejection method

  double dx = xl.max - xl.min;
  double dy = yl.max - yl.min;
  double gx=-1, gy=-1, gW=-1, gQ2=-1, xsec=-1;

  unsigned int iter = 0;
  bool accept  = false;
  bool doRad   = false;
  bool doneISR = false;
  if (fDoRadiativeCorrection) {
      doRad = true;
      doneISR = false;
  }

  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
       LOG("DISKinematics", pWARN)
         << " Couldn't select kinematics after " << iter << " iterations";
       evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select kinematics");
       exception.SwitchOnFastForward();
       throw exception;
     }

     //-- random x,y
     gx = xl.min + dx * rnd->RndKine().Rndm();
     gy = yl.min + dy * rnd->RndKine().Rndm();
     interaction->KinePtr()->Setx(gx);
     interaction->KinePtr()->Sety(gy);
     kinematics::UpdateWQ2FromXY(interaction);

     LOG("DISKinematics", pNOTICE)
        << "Trying: x = " << gx << ", y = " << gy
        << " (W  = " << interaction->KinePtr()->W()  << ","
        << " (Q2 = " << interaction->KinePtr()->Q2() << ")";

     //-- compute the cross section for current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        this->AssertXSecLimits(interaction, xsec, xsec_max);
        double t = xsec_max * rnd->RndKine().Rndm();
	double J = 1;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("DISKinematics", pDEBUG)
              << "xsec= " << xsec << ", J= " << J << ", Rnd= " << t;
#endif
        accept = (t < J*xsec);
     }
     else {
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
         LOG("DISKinematics", pNOTICE)
            << "Selected:  x = " << gx << ", y = " << gy
            << " (W  = " << interaction->KinePtr()->W()  << ","
            << " (Q2 = " << interaction->KinePtr()->Q2() << ")";

         // reset trust bits
         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec,kPSxyfE);

         // for uniform kinematics, compute an event weight as
         // wght = (phase space volume)*(differential xsec)/(event total xsec)
         if(fGenerateUniformly) {
            double vol     = kinematics::PhaseSpaceVolume(interaction,kPSxyfE);
            double totxsec = evrec->XSec();
            double wght    = (vol/totxsec)*xsec;
            LOG("DISKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

            // apply computed weight to the current event weight
            wght *= evrec->Weight();
            LOG("DISKinematics", pNOTICE) << "Current event wght = " << wght;
            evrec->SetWeight(wght);
         }

         // compute W,Q2 for selected x,y
         //bool is_em = interaction->ProcInfo().IsEM();
         kinematics::XYtoWQ2(Ev,M,gW,gQ2,gx,gy);

         LOG("DISKinematics", pNOTICE)
                        << "Selected x,y => W = " << gW << ", Q2 = " << gQ2;

         //-- In the case of Radiative correction, the radiated photon energy
         //   depends on the final state lepton kinematics.
         if(doRad && !doneISR) {
          LOG("DISKinematics", pINFO)  << "Doing ISR "<< fModel <<" cutoff "<<fCutoff<<" thickness "<<fThickness;
          TLorentzVector p4l = this->GetFinalStateLeptonKinematic(evrec,Ev,gy,gQ2);
          RadiativeCorrector * fISRCorrector = new RadiativeCorrector();
          fISRCorrector->SetISR(true);
	  fISRCorrector->SetDoInternalRad(fDoInternal);
          fISRCorrector->SetModel(fModel);
          fISRCorrector->SetCutoff(fCutoff);
          fISRCorrector->SetThickness(fThickness);
          fISRCorrector->SetQ2(gQ2);
          fISRCorrector->SetP4l(p4l);
          fISRCorrector->ProcessEventRecord(evrec);
          doneISR = true;
          continue;
        }

	 // lock selected kinematics & clear running values
         interaction->KinePtr()->SetW (gW,  true);
         interaction->KinePtr()->SetQ2(gQ2, true);
         interaction->KinePtr()->Setx (gx,  true);
         interaction->KinePtr()->Sety (gy,  true);
         interaction->KinePtr()->ClearRunningValues();
         return;
     }
  } // iterations
}
//___________________________________________________________________________
TLorentzVector DISKinematicsGenerator::GetFinalStateLeptonKinematic(GHepRecord * evrec, double E, double gy, double gQ2) const
{
          Interaction * interaction = evrec->Summary();
          const InitialState & init_state = interaction->InitState();
          LOG("DISKinematics", pNOTICE) << "Event selected for radiative correction Selected in kRfHitNucRest: Q^2 = " << gQ2<<" gy = ";
          const TLorentzVector & pnuc4 = init_state.Tgt().HitNucP4(); //[@LAB]
          TVector3 beta = pnuc4.BoostVector();
          double El  = (1-gy)*E;
          double ml  = interaction->FSPrimLepton()->Mass();
          double ml2 = TMath::Power(ml,2);
          double plp = El - 0.5*(gQ2+ml2)/E;                          // p(//)
          double plt = TMath::Sqrt(TMath::Max(0.,El*El-plp*plp-ml2)); // p(-|)
          LOG("DISKinematics", pNOTICE) << "Calulating temp final state for radiative correction fsl @ Nucleon rest frame: E = " << El << ", |p//| = " << plp << ", [pT] = " << plt;
	  // Randomize transverse components
	  RandomGen * rnd_temp = RandomGen::Instance();
          //double phi  = 2*kPi * rnd_temp->RndLep().Rndm();
          double phi  = 2* TMath::Pi() * rnd_temp->RndLep().Rndm();
          double pltx = plt * TMath::Cos(phi);
          double plty = plt * TMath::Sin(phi);
          TLorentzVector * p4v = evrec->CorrectProbe()->GetP4(); // v 4p @ LAB
          p4v->Boost(-1.*beta);                           // v 4p @ Nucleon rest frame
          // Take a unit vector along the neutrino direction @ the nucleon rest frame
          TVector3 unit_nudir = p4v->Vect().Unit();
          // Rotate lepton momentum vector from the reference frame (x'y'z') where 
	  // {z':(neutrino direction), z'x':(theta plane)} to the nucleon rest frame
	  TVector3 p3l(pltx,plty,plp);
          p3l.RotateUz(unit_nudir);
          // Lepton 4-momentum in the nucleon rest frame
          TLorentzVector p4l(p3l,El);
          p4l.Boost(beta); // active Lorentz transform
          LOG("DISKinematics", pNOTICE) << "Calulating temp final state for radiative correction fsl @ LAB: E " <<p4l.E() << " px "<<p4l.Px() << " py "<<p4l.Py() << " pz "<<p4l.Pz() ;
          return p4l;
}
//___________________________________________________________________________
void DISKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void DISKinematicsGenerator::LoadConfig(void)
{
// Reads its configuration data from its configuration Registry and loads them
// in private data members to avoid looking up at the Registry all the time.

  //-- Safety factor for the maximum differential cross section
	GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor,  1.25 ) ;

  //-- Minimum energy for which max xsec would be cached, forcing explicit
  //   calculation for lower eneries
	GetParamDef( "Cache-MinEnergy", fEMin, 0.8 ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
	GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
    assert(fMaxXSecDiffTolerance>=0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
    GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

    GetParam("doRadiativeCorrection", fDoRadiativeCorrection, false) ;
    if (fDoRadiativeCorrection) {
       GetParam( "RadiativeCorrectionModel" , fModel);
       GetParam( "RadiativeCorrectionCutoff",fCutoff);
       GetParam( "RadiativeCorrectionThickness",fThickness);
       GetParam( "RadiativeCorrectionDoInternal",fDoInternal);
    }

}
//____________________________________________________________________________
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

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISKinematics", pDEBUG)<< "Computing max xsec in allowed phase space";
#endif
  double max_xsec = 0.0;

  const InitialState & init_state = interaction->InitState();
  //double Ev = init_state.ProbeE(kRfHitNucRest);
  //const ProcessInfo & proc = interaction->ProcInfo();
  const Target & tgt = init_state.Tgt();

  int Ny  = 20;
  int Nx  = 40;
  int Nxb = 3;

  double xpeak    = .2;
  double xwindow  = .2;
  double ypeak    = .5;
  double ywindow  = .5;

  if(tgt.HitQrkIsSet()) {
    if(tgt.HitSeaQrk()) {
       xpeak    = .1;
       xwindow  = .1;
       ypeak    = .7;
       ywindow  = .3;
    } else {
       xpeak    = .2;
       xwindow  = .2;
       ypeak    = .7;
       ywindow  = .3;
    }
  }

  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t xl = kps.Limits(kKVx);
  Range1D_t yl = kps.Limits(kKVy);

  double xmin    = TMath::Max(xpeak-xwindow, TMath::Max(xl.min, 5E-3));
  double xmax    = TMath::Min(xpeak+xwindow, xl.max);
  double ymin    = TMath::Max(ypeak-ywindow, TMath::Max(yl.min, 2E-3));
  double ymax    = TMath::Min(ypeak+ywindow, yl.max);
  double dx      = (xmax-xmin)/(Nx-1);
  double dy      = (ymax-ymin)/(Ny-1);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  LOG("DISKinematics", pDEBUG)
    << "Searching max. in x [" << xmin << ", " << xmax << "], y [" << ymin << ", " << ymax << "]";
#endif
  double xseclast_y = -1;
  bool increasing_y;

  for(int i=0; i<Ny; i++) {
     double gy = ymin + i*dy;
     //double gy = TMath::Power(10., logymin + i*dlogy);
     interaction->KinePtr()->Sety(gy);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("DISKinematics", pDEBUG) << "y = " << gy;
#endif
     double xseclast_x = -1;
     bool increasing_x;

     for(int j=0; j<Nx; j++) {
        double gx = xmin + j*dx;
	//double gx = TMath::Power(10., logxmin + j*dlogx);
        interaction->KinePtr()->Setx(gx);
        kinematics::UpdateWQ2FromXY(interaction);

        double xsec = fXSecModel->XSec(interaction, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("DISKinematics", pINFO)
                << "xsec(y=" << gy << ", x=" << gx << ") = " << xsec;
#endif
        // update maximum xsec
        max_xsec = TMath::Max(xsec, max_xsec);

        increasing_x = xsec-xseclast_x>=0;
        xseclast_x   = xsec;

        // once the cross section stops increasing, I reduce the step size and
        // step backwards a little bit to handle cases that the max cross section
        // is grossly underestimated (very peaky distribution & large step)
        if(!increasing_x) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
          LOG("DISKinematics", pDEBUG)
           << "d2xsec/dxdy|x stopped increasing. Stepping back & exiting x loop";
#endif
          //double dlogxn = dlogx/(Nxb+1);
          double dxn = dx/(Nxb+1);
          for(int ik=0; ik<Nxb; ik++) {
	     //gx = TMath::Exp(TMath::Log(gx) - dlogxn);
   	     gx = gx - dxn;
             interaction->KinePtr()->Setx(gx);
             kinematics::UpdateWQ2FromXY(interaction);
             xsec = fXSecModel->XSec(interaction, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
             LOG("DISKinematics", pINFO)
                << "xsec(y=" << gy << ", x=" << gx << ") = " << xsec;
#endif
	  }
          break;
        } // stepping back within last bin
     } // x
     increasing_y = max_xsec-xseclast_y>=0;
     xseclast_y   = max_xsec;
     if(!increasing_y) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       LOG("DISKinematics", pDEBUG)
           << "d2xsec/dxdy stopped increasing. Exiting y loop";
#endif
       break;
     }
  }// y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy
  //  max_xsec *= fSafetyFactor;
  //max_xsec *= ( (Ev<3.0) ? 2.5 : fSafetyFactor);
  max_xsec *= 3;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("DISKinematics", pDEBUG) << interaction->AsString();
  SLOG("DISKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("DISKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;
#endif

  return max_xsec;
}
//___________________________________________________________________________
