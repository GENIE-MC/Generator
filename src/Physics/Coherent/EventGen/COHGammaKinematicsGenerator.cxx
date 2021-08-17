//____________________________________________________________________________
/*
  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda < mroda \at liverpool.ac.uk>
  University of Liverpool 

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
#include <iostream>
#include <functional>   
#include <numeric>      // std::accumulate


#include <TROOT.h>
#include <TMath.h>
#include <TF2.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Physics/Coherent/EventGen/COHGammaKinematicsGenerator.h"
#include "Framework/Conventions/KinePhaseSpace.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/XSectionIntegration/GSLXSecFunc.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace genie::utils;

//___________________________________________________________________________
COHGammaKinematicsGenerator::COHGammaKinematicsGenerator() :
  KineGeneratorWithCache("genie::COHGammaKinematicsGenerator"), 
  fGammaLimits( nullptr ), 
  ftPhaseSpace( true )
{ 

}
//___________________________________________________________________________
COHGammaKinematicsGenerator::COHGammaKinematicsGenerator(string config) :
  KineGeneratorWithCache("genie::COHGammaKinematicsGenerator", config), 
  fGammaLimits( nullptr ), 
  ftPhaseSpace( true )
{

}
//___________________________________________________________________________
COHGammaKinematicsGenerator::~COHGammaKinematicsGenerator()
{

}
//___________________________________________________________________________
void COHGammaKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  if(fGenerateUniformly) {
    LOG("COHGammaKinematicsGenerator", pNOTICE)
      << "Generating kinematics uniformly over the allowed phase space";
  }

  Interaction * in = evrec -> Summary() ;

  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();
  
  if ( ! fXSecModel -> ValidProcess( in )  ) {
    std::stringstream message = "Cannot calculate kinematics for " ;
    message << fXSecModel->Id().Name();
 
    LOG("COHGammaKinematicsGenerator",pFATAL) << message.str() ;

    exceptions::EVGThreadException ex;
    ex.SetReason( message.str() );
    ex.SwitchOnFastForward() ;
    throw ex;
  }

  in -> SetBit(kISkipProcessChk);
  //in -> SetBit(kISkipKinematicChk);

  double xsec_max = (fGenerateUniformly) ? -1 : this->MaxXSec(evrec); 

  // Initialise a random number generator 
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant

  std::array<Range1D_t,4> ranges = { fGammaLimits -> EGamma( *in ), 
                                     ftPhaseSpace ? fGammaLimits -> t( *in ) : fGammaLimits -> ThetaLepton( *in ), 
                                     fGammaLimits -> ThetaGamma( *in ), 
                                     fGammaLimits -> PhiGamma( *in ) } ;
  
  std::array<double, 4> widths ;
  for ( unsigned int i = 0 ; i < ranges.size() ; ++i ) {
    widths[i] = ranges[i].max - ranges[i].min ; 
  }
  
  //------ Try to select a valid set of kinematics
  unsigned int iter = 0;
  bool accept=false;
  
  double xsec=-1 ; 
  std::array<double,4> point ; 
  
  Interaction local_interaction( * in ) ;

  ROOT::Math::IBaseFunctionMultiDim * func = nullptr ;
  if ( ftPhaseSpace ) func = new utils::gsl::d4Xsec_dEgdtdThetagdPhig( fXSecModel, & local_interaction ) ;
  else func = new utils::gsl::d4Xsec_dEgdThetaldThetagdPhig( fXSecModel, & local_interaction ) ;
  
  while(1) {
    iter++;
    
    if( iter > fMaxIterations )  this->throwOnTooManyIterations(iter,evrec);
    
    // generate the kinematic point
    for ( unsigned int i = 0 ; i < ranges.size() ; ++i ) {
      point[i] = ranges[i].min + widths[i] * rnd->RndKine().Rndm();
    }

    // LOG("COHKinematics", pINFO) << "Trying: Gamma(" << g_E_g << ", "
    // 				<< g_theta_g << ", " << g_phi_g << "),   Lep(" 
    // 				<< g_theta_l << ")";
    
    xsec = (*func)( point.data() ) ; 
    
    //realized this method assumes E_l 
    
    if (!fGenerateUniformly) {
      //-- decide whether to accept the current kinematics
      double t   = xsec_max * rnd->RndKine().Rndm();

      LOG("COHKinematics", pINFO) << "Got: xsec = " << xsec << ", t = " << 
        t << " (max_xsec = " << xsec_max << ")";
      
      this->AssertXSecLimits( in, xsec, xsec_max );
      
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

      // update kinematic variables
      Kinematics * kine = in -> KinePtr() ;
      const Kinematics & local_kine = local_interaction.Kine() ;
      
      kine -> Setx ( local_kine.x() , true );
      kine -> Sety ( local_kine.y() , true );
      kine -> SetQ2( local_kine.Q2(), true );
      kine -> SetW ( local_kine.W() , true );
      kine -> Sett ( local_kine.t() , true );
      
      // generate the phase for the lepton
      double phi_l =  2*constants::kPi * rnd->RndKine().Rndm(); // final overall roation added

      // add the phase to both lepton and gamma so that the relative is the same
      TLorentzVector lep ( local_kine.FSLeptonP4() ) ;
      lep.SetPhi( phi_l ) ;
      kine -> SetFSLeptonP4( lep ) ;

      TLorentzVector gamma( local_kine.HadSystP4() ) ;
      double phi_g = gamma.Phi() ;
      phi_g += phi_l ;
      gamma.SetPhi( phi_g ) ;
      kine -> SetHadSystP4( gamma ) ;
      
      // for uniform kinematics, compute an event weight as
      // wght = (phase space volume)*(differential xsec)/(event total xsec)
      if(fGenerateUniformly) {
        // Phase space volume needs checking
	double vol = std::accumulate( widths.begin(), widths.end(), 1., std::multiplies<double>() ) ;
        double totxsec = evrec->XSec();
        double wght    = (vol/totxsec)*xsec;
        LOG("COHKinematics", pNOTICE)  << "Kinematics wght = "<< wght;
	
        // apply computed weight to the current event weight
        wght *= evrec->Weight();
        LOG("COHKinematics", pNOTICE) << "Current event wght = " << wght;
        evrec->SetWeight(wght);
      }

      evrec->SetDiffXSec(xsec,kPSEgTlTgPgfE);
      
      // reset bits
      in -> ResetBit(kISkipProcessChk);
      in -> ResetBit(kISkipKinematicChk);
      
      break ; 
    }//if accept
		  
  }//while still throwing events
  
  delete func ;   

  return;
}
//___________________________________________________________________________
double COHGammaKinematicsGenerator::ComputeMaxXSec(const Interaction * in) const
{
  // Computes the maximum differential cross section in the requested phase
  // space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
  // method and the value is cached at a circular cache branch for retrieval
  // during subsequent event generation.

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
  SLOG("COHKinematics", pDEBUG)
    << "Scanning the allowed phase space {K} for the max(dxsec/d{K})";
#endif
  if ( ! fXSecModel -> ValidProcess( in )  ) {
    LOG("COHGammaKinematicsGenerator",pFATAL) << "Cannot calculate max-xsec for " 
					      << fXSecModel->Id().Name();
    exit(0) ;
  }
  
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer("Minuit2" );
  //  min -> SetPrintLevel(3) ;
  
  gsl::d4Xsec_dEgdThetaldThetagdPhig f(fXSecModel,in);
  f.SetFactor(-1.); // Make it return negative of cross-section so we can minimize
  
  min->SetFunction( f );
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(0.05);//not sure what this does....
  

  //need min, max, n, and d for all the xsec varables
  //for COH gamma they are Eg, theta_l theta_gamma, phi_gamma 
  //then set variables and get min
  
  std::array<string, 4> names = { "E_g", 
				  ftPhaseSpace ? "t" : "theta_l", 
				  "theta_g", 
				  "phi_g" } ;

  std::array<Range1D_t,4> ranges = { fGammaLimits -> EGamma( *in ), 
				     ftPhaseSpace ? fGammaLimits -> t( *in ) : fGammaLimits -> ThetaLepton( *in ), 
				     fGammaLimits -> ThetaGamma( *in ), 
				     fGammaLimits -> PhiGamma( *in ) } ;
  
  std::array<double,4> start, steps, temp_point ;
  // Please not that if Minuit2 minimizer is used, the steps are not used
  // but for consistency we are evaluating it
  
  for ( unsigned int i = 0 ; i < ranges.size() ; ++i ) {
    double width = ranges[i].max - ranges[i].min ;
    steps[i] = width / ( fMinimScanPoints[i] +1 ) ;
  }
  
  double xsec = 0; 
  
  // preliimnary scan 
  for ( unsigned int i = 1 ; i <= fMinimScanPoints[0] ; ++i ) {
    temp_point[0] = ranges[0].min + steps[0]*i ;

    for ( unsigned int j = 1 ; j <= fMinimScanPoints[1] ; ++j ) {
      temp_point[1] = ranges[1].min + steps[1]*j ;
      
      for ( unsigned int k = 1 ; k <= fMinimScanPoints[2] ; ++k ) {
	temp_point[2] = ranges[2].min + steps[2]*k ;

	for ( unsigned int l = 1 ; l <= fMinimScanPoints[3] ; ++l ) {
	  temp_point[3] = ranges[3].min + steps[3]*l ;

	  double temp_xsec = - f( temp_point.data() ) ;
	  if ( temp_xsec > xsec ) {
	    start = temp_point ;
	    xsec = temp_xsec ;
	  }
	  
	}
      }
    }
  }

 
  for ( unsigned int i = 0 ; i < ranges.size() ; ++i ) {
    min -> SetLimitedVariable( i, names[i], start[i], steps[i], ranges[i].min, ranges[i].max ) ;
  }
  
  min->Minimize();
  
  double max_xsec = -min->MinValue(); //back to positive xsec

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
double COHGammaKinematicsGenerator::Energy(const Interaction * interaction) const
{
  // Override the base class Energy() method to cache the max xsec for the
  // neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void COHGammaKinematicsGenerator::throwOnTooManyIterations(unsigned int iters,
                                                      GHepRecord* evrec) const
{
  LOG("COHKinematics", pWARN)
    << "*** Could not select valid kinematics after "
    << iters << " iterations";
  evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
  genie::exceptions::EVGThreadException exception;
  exception.SetReason("Couldn't select kinematics");
  exception.SwitchOnFastForward();
  throw exception;
}
//___________________________________________________________________________
void COHGammaKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHGammaKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHGammaKinematicsGenerator::LoadConfig(void)
{

  bool error = false ;

  GetParamDef( "tPhaseSpace", ftPhaseSpace, false ) ;
  
  
  //-- max xsec safety factor (for rejection method) and min cached energy
  GetParamDef( "MaxXSec-SafetyFactor", fSafetyFactor, 1.6 ) ;
  GetParamDef( "Cache-MinEnergy", fEMin,  -1.0 ) ;

  //-- Generate kinematics uniformly over allowed phase space and compute
  //   an event weight?
  GetParamDef( "UniformOverPhaseSpace", fGenerateUniformly, false ) ;

  //-- Maximum allowed fractional cross section deviation from maxim cross
  //   section used in rejection method
  GetParamDef( "MaxXSec-DiffTolerance", fMaxXSecDiffTolerance, 999999. ) ;
  assert(fMaxXSecDiffTolerance>=0);

  int max_iter ;
  GetParamDef( "MaxIterations", max_iter, (int) kRjMaxIterations ) ;
  if ( max_iter > 0 ) {
    fMaxIterations = max_iter ;
  }

  const Algorithm * temp = SubAlg( "IntegrationLimits" ) ;
  fGammaLimits = dynamic_cast<const COHGammaIntegrationLimits *>( temp ) ;
  if (! fGammaLimits ) {
    LOG( "COHGammaKinematicsGenerator", pERROR ) << "Gamma integration limits subalgo failed to load" ;
    error = true ;
  }

  std::vector<int> scan_points ;
  GetParamVect( "MinimScanPoint", scan_points ) ;
  if ( scan_points.size() < fMinimScanPoints.size() ) {
    LOG( "COHGammaKinematicsGenerator", pERROR ) << "Not enough information for phase space scan" ;
    error = true ;
  }
		
  for ( unsigned int i = 0 ; i < fMinimScanPoints.size() ; ++i ) {
    fMinimScanPoints[i] = std::max( 1, scan_points[i] ) ;
  }
  
  if ( error ) {
    LOG( "COHGammaKinematicsGenerator", pFATAL ) << "Invalid configuration. Exiting" ;
    exit( 78 ) ;
  } 

}
//____________________________________________________________________________

