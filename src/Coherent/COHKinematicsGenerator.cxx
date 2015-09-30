//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 09, 2009 - CA
   Moved into the new Coherent package from its previous location  (EVGModules 
   package)
 @ Mar 03, 2009 - CA
   Renamed COHPiKinematicsGenerator -> COHKinematicsGenerator in
   anticipation of reusing the code for simulating coherent production of
   vector mesons.
 @ May 06, 2009 - CA
   Fix a problem with the search for the max cross section over the allowed
   phase space which prevented kinematics to be generated for events near the 
   energy threshold.
 @ Feb 06, 2013 - CA
   When the value of the differential cross-section for the selected kinematics
   is set to the event, set the corresponding KinePhaseSpace_t value too.

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>
#include <TF2.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/GBuild.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/Units.h"
#include "Coherent/COHKinematicsGenerator.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
COHKinematicsGenerator::COHKinematicsGenerator() :
KineGeneratorWithCache("genie::COHKinematicsGenerator")
{
  fEnvelope = 0;
}
//___________________________________________________________________________
COHKinematicsGenerator::COHKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::COHKinematicsGenerator", config)
{
  fEnvelope = 0;
}
//___________________________________________________________________________
COHKinematicsGenerator::~COHKinematicsGenerator()
{
  if(fEnvelope) delete fEnvelope;
}
//___________________________________________________________________________
void COHKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("COHKinematics", pNOTICE)
          << "Generating kinematics uniformly over the allowed phase space";
  }

  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);
  interaction->SetBit(kISkipKinematicChk);

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec);

  //-- Get the kinematical limits for the generated x,y
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t y = kps.YLim();
  assert(y.min>0. && y.max>0. && y.min<1. && y.max<1. && y.min<y.max);

  const double xmin = kASmallNum;
  const double xmax = 1.- kASmallNum;
  const double ymin = y.min + kASmallNum;
  const double ymax = y.max - kASmallNum;
  const double dx   = xmax - xmin;
  const double dy   = ymax - ymin;

  //------ Try to select a valid x,y pair

  unsigned int iter = 0;
  bool accept=false;
  double xsec=-1, gx=-1, gy=-1;

  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
        LOG("COHKinematics", pWARN)
             << "*** Could not select a valid (x,y) pair after "
                                               << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     if(fGenerateUniformly) {
        //-- Generate a x,y pair uniformly in the kinematically allowed range.
        gx = xmin + dx * rnd->RndKine().Rndm();
        gy = ymin + dy * rnd->RndKine().Rndm();

     } else {
        //-- Select unweighted kinematics using importance sampling method. 

        if(iter==1) {
         LOG("COHKinematics", pNOTICE) << "Initializing the sampling envelope";
         double Ev = interaction->InitState().ProbeE(kRfLab);
         fEnvelope->SetRange(xmin,ymin,xmax,ymax);
         fEnvelope->SetParameter(0, xsec_max);  
         fEnvelope->SetParameter(1, Ev);        
       }

       // Generate W,QD2 using the 2-D envelope as PDF
       fEnvelope->GetRandom2(gx,gy);
     }

     LOG("COHKinematics", pINFO) << "Trying: x = " << gx << ", y = " << gy;

     interaction->KinePtr()->Setx(gx);
     interaction->KinePtr()->Sety(gy);

     // computing cross section for the current kinematics
     xsec = fXSecModel->XSec(interaction, kPSxyfE);

     //-- decide whether to accept the current kinematics
     if(!fGenerateUniformly) {
        double max = fEnvelope->Eval(gx, gy);
        double t   = max * rnd->RndKine().Rndm();

        this->AssertXSecLimits(interaction, xsec, max);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
        LOG("COHKinematics", pDEBUG) 
            << "xsec= " << xsec << ", J= 1, Rnd= " << t;
#endif
        accept = (t<xsec);
     }
     else { 
        accept = (xsec>0);
     }

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
        LOG("COHKinematics", pNOTICE) << "Selected: x = "<< gx << ", y = "<< gy;

        // the COH cross section should be a triple differential cross section
        // d^2xsec/dxdydt where t is the the square of the 4p transfer to the
        // nucleus. The cross section used for kinematical selection should have
        // the t-dependence integrated out. The t-dependence is of the form
        // ~exp(-bt). Now that the x,y kinematical variables have been selected
        // we can generate a t using the t-dependence as a PDF.
        const InitialState & init_state = interaction->InitState();
        double Ev    = init_state.ProbeE(kRfLab);
        double Epi   = gy*Ev; // pion energy
        double Epi2  = TMath::Power(Epi,2);
        double pme2  = kPionMass2/Epi2;   
        double xME   = kNucleonMass*gx/Epi;
        double tA    = 1. + xME - 0.5*pme2;
        double tB    = TMath::Sqrt(1.+ 2*xME) * TMath::Sqrt(1.-pme2);
        double tmin  = 2*Epi2 * (tA-tB);
        double tmax  = 2*Epi2 * (tA+tB);
        double A     = (double) init_state.Tgt().A(); 
        double A13   = TMath::Power(A,1./3.);
        double R     = fRo * A13 * units::fermi; // nuclear radius
        double R2    = TMath::Power(R,2.);
        double b     = 0.33333 * R2;
        double tsum  = (TMath::Exp(-b*tmin) - TMath::Exp(-b*tmax))/b; 
        double rt    = tsum * rnd->RndKine().Rndm();
        double gt    = -1.*TMath::Log(-1.*b*rt + TMath::Exp(-1.*b*tmin))/b;

        LOG("COHKinematics", pNOTICE)
          << "Selected: t = "<< gt << ", from ["<< tmin << ", "<< tmax << "]";

        // for uniform kinematics, compute an event weight as
        // wght = (phase space volume)*(differential xsec)/(event total xsec)
        if(fGenerateUniformly) {
          double vol     = y.max-y.min; // dx=1, dt: irrelevant
          double totxsec = evrec->XSec();
          double wght    = (vol/totxsec)*xsec;
          LOG("COHKinematics", pNOTICE)  << "Kinematics wght = "<< wght;

          // apply computed weight to the current event weight
          wght *= evrec->Weight();
          LOG("COHKinematics", pNOTICE) << "Current event wght = " << wght;
          evrec->SetWeight(wght);
        }

        // reset bits
        interaction->ResetBit(kISkipProcessChk);
        interaction->ResetBit(kISkipKinematicChk);

        // lock selected kinematics & clear running values
        interaction->KinePtr()->Setx(gx, true);
        interaction->KinePtr()->Sety(gy, true);
        interaction->KinePtr()->Sett(gt, true);
        interaction->KinePtr()->SetW(kPionMass, true);
        interaction->KinePtr()->SetQ2(2*kNucleonMass*gx*gy*Ev, true);
        interaction->KinePtr()->ClearRunningValues();

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec*TMath::Exp(-b*gt)/tsum,kPSxytfE);

        return;
     }
  }// iterations
}
//___________________________________________________________________________
double COHKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG)
          << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";
#endif
  double max_xsec = 0.;

  double Ev = in->InitState().ProbeE(kRfLab);

  const int Nx = 50;
  const int Ny = 50;

  const KPhaseSpace & kps = in->PhaseSpace();
  Range1D_t y = kps.YLim();

  const double logxmin = TMath::Log10(1E-5);
  const double logxmax = TMath::Log10(1.0);
  const double logymin = TMath::Log10(y.min);
  const double logymax = TMath::Log10(y.max);
/*
  double dy=0;
  double log10Ev = TMath::Log10(Ev);
  double yc = TMath::Power(10,-0.5813-0.8492*log10Ev);
  const double logymin = TMath::Log10( TMath::Max(y.min,yc-dy) );
  const double logymax = TMath::Log10( TMath::Min(y.max,yc+dy) );
*/
  const double dlogx   = (logxmax - logxmin) /(Nx-1);
  const double dlogy   = (logymax - logymin) /(Ny-1);

  for(int i=0; i<Nx; i++) {
   double gx = TMath::Power(10, logxmin + i * dlogx);
   for(int j=0; j<Ny; j++) {
     double gy = TMath::Power(10, logymin + j * dlogy);

     double Q2 = 2*kNucleonMass*gx*gy*Ev;
     if(Ev>1.0 && Q2>0.01) continue;

     in->KinePtr()->Setx(gx);
     in->KinePtr()->Sety(gy);

     double xsec = fXSecModel->XSec(in, kPSxyfE);
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
     LOG("COHKinematics", pDEBUG)  
     	 << "xsec(x= " << gx << ", y= " << gy << ") = " << xsec;
#endif
     max_xsec = TMath::Max(max_xsec, xsec);

   }//y
  }//x

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG) << in->AsString();
  SLOG("COHKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("COHKinematics", pDEBUG) << "Computed using alg = " << fXSecModel->Id();
#endif

  return max_xsec;
}
//___________________________________________________________________________
double COHKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void COHKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHKinematicsGenerator::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  //-- COH model parameter Ro
  fRo = fConfig->GetDoubleDef("Ro", gc->GetDouble("COH-Ro"));

  //-- max xsec safety factor (for rejection method) and min cached energy
  fSafetyFactor = fConfig->GetDoubleDef("MaxXSec-SafetyFactor", 1.6);
  fEMin         = fConfig->GetDoubleDef("Cache-MinEnergy",     -1.0);

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  fGenerateUniformly = fConfig->GetBoolDef("UniformOverPhaseSpace", false);

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  fMaxXSecDiffTolerance = 
         fConfig->GetDoubleDef("MaxXSec-DiffTolerance",999999.);
  assert(fMaxXSecDiffTolerance>=0);

  //-- Envelope employed when importance sampling is used 
  //   (initialize with dummy range)
  if(fEnvelope) delete fEnvelope;
  fEnvelope = new TF2("envelope",
    	  kinematics::COHImportanceSamplingEnvelope,0.,1,0.,1,2);
}
//____________________________________________________________________________

