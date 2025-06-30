//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Alfonso Garcia <alfonsog \at nikhef.nl>
         NIKHEF

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Physics/HEDIS/EventGen/HEDISKinematicsGenerator.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/KineVar.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/Range1.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

#include <TMath.h>

#include <string>

using namespace genie;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
HEDISKinematicsGenerator::HEDISKinematicsGenerator() :
KineGeneratorWithCache("genie::HEDISKinematicsGenerator")
{

}
//___________________________________________________________________________
HEDISKinematicsGenerator::HEDISKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::HEDISKinematicsGenerator", config)
{

}
//___________________________________________________________________________
HEDISKinematicsGenerator::~HEDISKinematicsGenerator()
{

}
//___________________________________________________________________________
void HEDISKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

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
  double Ev  = init_state.ProbeE(kRfLab);
  double M   = init_state.Tgt().HitNucP4().M(); // can be off m-shell

  //-- Get the physical W range 
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t xl = kps.XLim();
  Range1D_t Q2l = kps.Q2Lim();

  //-- x and y lower limit restrict by limits in SF tables
  Q2l.min = TMath::Max(Q2l.min,fSFQ2min);
  Q2l.max = TMath::Min(Q2l.max,fSFQ2max);
  xl.min  = TMath::Max(TMath::Max(xl.min,Q2l.min/2./M/Ev),fSFXmin);

  LOG("HEDISKinematics", pNOTICE) << "x: [" << xl.min << ", " << xl.max << "]"; 
  LOG("HEDISKinematics", pNOTICE) << "log10Q2: [" << Q2l.min << ", " << Q2l.max << "]"; 

  //-- For the subsequent kinematic selection with the rejection method:
  
  //Scan through a wide region to find the maximum
  Range1D_t xrange_wide(xl.min*fWideRange,xl.max/fWideRange); 
  Range1D_t Q2range_wide(Q2l.min*fWideRange,Q2l.max/fWideRange); 
  double x_wide    = 0.;
  double Q2_wide   = 0.;
  double xsec_wide = this->Scan(interaction,xrange_wide,Q2range_wide,fWideNKnotsX,fWideNKnotsQ2,2*M*Ev,x_wide,Q2_wide);

  //Scan through a fine region to find the maximum
  Range1D_t xrange_fine(TMath::Max(x_wide/fFineRange,xrange_wide.min),TMath::Min(x_wide*fFineRange,xrange_wide.max)); 
  Range1D_t Q2range_fine(TMath::Max(Q2_wide/fFineRange,Q2range_wide.min),TMath::Min(Q2_wide*fFineRange,Q2range_wide.max)); 
  double x_fine    = 0.;
  double Q2_fine   = 0.;
  double xsec_fine = this->Scan(interaction,xrange_fine,Q2range_fine,fFineNKnotsX,fFineNKnotsQ2,2*M*Ev,x_fine,Q2_fine);

  //Apply safety factor
  double xsec_max = fSafetyFactor * TMath::Max(xsec_wide,xsec_fine);

  //-- Try to select a valid (x,y) pair using the rejection method
  double log10xmin  = TMath::Log10(xl.min);  
  double log10xmax  = TMath::Log10(xl.max); 
  double log10Q2min = TMath::Log10(Q2l.min);  
  double log10Q2max = TMath::Log10(Q2l.max); 

  double dlog10x = log10xmax - log10xmin; 
  double dlog10Q2 = log10Q2max - log10Q2min; 
  
  double gx=-1, gQ2=-1, xsec=-1;

  unsigned int iter = 0;
  bool accept = false;
  while(1) {
     iter++;
     if(iter > kRjMaxIterations) {
       LOG("HEDISKinematics", pWARN)
         << " Couldn't select kinematics after " << iter << " iterations";
       evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select kinematics");
       exception.SwitchOnFastForward();
       throw exception;
     }
    
     gx = TMath::Power( 10., log10xmin + dlog10x * rnd->RndKine().Rndm() ); 
     gQ2 = TMath::Power( 10., log10Q2min + dlog10Q2 * rnd->RndKine().Rndm() ); 

     interaction->KinePtr()->Setx(gx);
     interaction->KinePtr()->SetQ2(gQ2);
     kinematics::UpdateWYFromXQ2(interaction);

     LOG("HEDISKinematics", pDEBUG) 
        << "Trying: x = " << gx << ", Q2 = " << gQ2 
        << " (W  = " << interaction->KinePtr()->W()  << ","
        << "  y = " << interaction->KinePtr()->y() << ")";

     //-- compute the cross section for current kinematics
     xsec = fXSecModel->XSec(interaction, kPSlog10xlog10Q2fE);

     //-- decide whether to accept the current kinematics
     this->AssertXSecLimits(interaction, xsec, xsec_max);

     double t = xsec_max * rnd->RndKine().Rndm();

     LOG("HEDISKinematics", pDEBUG) << "xsec= " << xsec << ", Rnd= " << t;

     accept = (t < xsec);

     //-- If the generated kinematics are accepted, finish-up module's job
     if(accept) {
         LOG("HEDISKinematics", pNOTICE) 
            << "Selected:  x = " << gx << ", Q2 = " << gQ2
            << " (W  = " << interaction->KinePtr()->W()  << ","
            << " (Y = " << interaction->KinePtr()->y() << ")";

         // reset trust bits
         interaction->ResetBit(kISkipProcessChk);
         interaction->ResetBit(kISkipKinematicChk);

         // set the cross section for the selected kinematics
         evrec->SetDiffXSec(xsec,kPSxQ2fE);

         // lock selected kinematics & clear running values
         interaction->KinePtr()->SetW (interaction->KinePtr()->W(),  true);
         interaction->KinePtr()->Sety (interaction->KinePtr()->y(),  true);
         interaction->KinePtr()->SetQ2(gQ2, true);
         interaction->KinePtr()->Setx (gx,  true);
         interaction->KinePtr()->ClearRunningValues();
         return;
     }
  } // iterations
}
//___________________________________________________________________________
double HEDISKinematicsGenerator::Scan(Interaction * interaction, Range1D_t xrange,Range1D_t Q2range, int NKnotsQ2, int NKnotsX, double ME2, double & x_scan, double & Q2_scan) const
{

  double xsec_scan = 0.;
  Q2_scan   = 0.;
  x_scan    = 0.;

  double stepQ2 = TMath::Power(Q2range.max/Q2range.min,1./(NKnotsQ2-1));

  for ( int iq=0; iq<NKnotsQ2; iq++) {
    double Q2_aux = Q2range.min*TMath::Power(stepQ2,iq);
    xrange.min = TMath::Max(xrange.min,Q2_aux/ME2);
    double stepx = TMath::Power(xrange.max/xrange.min,1./(NKnotsX-1));
    for ( int ix=0; ix<NKnotsX; ix++) {
      double x_aux = xrange.min*TMath::Power(stepx,ix);
      interaction->KinePtr()->Setx(x_aux);
      interaction->KinePtr()->SetQ2(Q2_aux);
      kinematics::UpdateWYFromXQ2(interaction);      
      double xsec_aux = fXSecModel->XSec(interaction, kPSlog10xlog10Q2fE);
      LOG("HEDISKinematics", pDEBUG) << "x = " << x_aux << " , Q2 = " << Q2_aux << ", xsec = " << xsec_aux; 
      if (xsec_aux>xsec_scan) {
        xsec_scan = xsec_aux;
        Q2_scan   = Q2_aux;
        x_scan    = x_aux;
      }
    }
  }

  LOG("HEDISKinematics", pNOTICE) << "scan -> x = " << x_scan << " , Q2 = " << Q2_scan << ", xsec = " << xsec_scan; 

  return xsec_scan;

}
//___________________________________________________________________________
double HEDISKinematicsGenerator::ComputeMaxXSec(const Interaction * /* interaction */ ) const
{
  return 0;
}
//___________________________________________________________________________
void HEDISKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void HEDISKinematicsGenerator::LoadConfig(void)
{

  //-- Safety factor for the maximum differential cross section
  GetParamDef("MaxXSec-SafetyFactor", fSafetyFactor, 2. ) ;
  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef("MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
  assert(fMaxXSecDiffTolerance>=0);

  GetParamDef("ScanWide-NKnotsX",   fWideNKnotsX,   10 ) ;
  GetParamDef("ScanWide-NKnotsQ2", fWideNKnotsQ2,   10 ) ;
  GetParamDef("ScanWide-Range",       fWideRange,  1.1 ) ;
  GetParamDef("ScanFine-NKnotsX",   fFineNKnotsX,   10 ) ;
  GetParamDef("ScanFine-NKnotsQ2", fFineNKnotsQ2,   10 ) ;
  GetParamDef("ScanFine-Range",       fFineRange,  10. ) ;

  // Limits from the SF tables that are useful to reduce computation 
  // time of the max cross section
  GetParam("XGrid-Min",  fSFXmin ) ;
  GetParam("Q2Grid-Min", fSFQ2min ) ;
  GetParam("Q2Grid-Max", fSFQ2max ) ;

}
