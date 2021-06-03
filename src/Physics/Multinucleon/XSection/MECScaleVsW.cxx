//_________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J.Tena and M.Roda
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Multinucleon/XSection/MECScaleVsW.h"
#include "Framework/Utils/StringUtils.h" 
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;

//_________________________________________________________________________
MECScaleVsW::MECScaleVsW() : 
  XSecScaleMap()
{
  
}
//_________________________________________________________________________
MECScaleVsW::MECScaleVsW(string config) : 
  XSecScaleMap(config) 
{
  
}
//_________________________________________________________________________
MECScaleVsW::~MECScaleVsW()
{

}
//_________________________________________________________________________
void MECScaleVsW::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void MECScaleVsW::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________
double MECScaleVsW::GetScaling( const Interaction & interaction ) const
{
  double Q0 = interaction.Kine().GetKV(kKVQ0) ; 
  double Q3 = interaction.Kine().GetKV(kKVQ3) ; 

  return GetScaling( Q0, Q3 ) ; 
}
//_________________________________________________________________________
double MECScaleVsW::GetScaling( const double Q0, const double Q3 ) const
{
  // Get the vectors that include the kinematic limits of W for a given event:
  std::vector<double> W_limits, weights ; 
  GetVectorWithLimits( W_limits, weights,Q0, Q3, fDefaultWeight ) ;

  // The Scaling is done using the "experimenter's W", which assumes a single nucleon
  // See motivation in : https://arxiv.org/pdf/1601.02038.pdf
  // Calculate event W:
  double Mn = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  double W = sqrt( pow(Mn,2) + 2*Mn*Q0 - pow(Q0,2) + pow(Q0,2) ) ;

  // Calculate scaling:

  for ( unsigned int i = 0 ; i < weights.size() - 1 ; ++i ) {
    if ( W >= W_limits[i] && W < W_limits[i+1] ) {
      return ScaleFunction( W, W_limits[i], W_limits[i+1], weights[i], weights[i+1] ) ; 
    } 
  }

  return 1. ; 
} 

//_________________________________________________________________________
void MECScaleVsW::GetVectorWithLimits( std::vector<double> & W_limits, std::vector<double> & weights,
				       const double Q0, const double Q3 , const double weight ) const {
  // This function is responsible to add the phase space limits in the WValues vector in case they are not included
  // in the configuration setup.

  double Mn = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  double W_min = sqrt( pow(Mn,2) + 2*Mn*fW1_Q0Q3_limits->Eval(Q3) - pow(Q3,2) + pow(fW1_Q0Q3_limits->Eval(Q3),2) ) ;
  double W_max = sqrt( pow(Mn,2) + 2*Mn*Q0 ) ; // Imposing Q2 = 0 

  W_limits = fWValues ; 
  weights = fWeights ; 

  if( W_limits[W_limits.size()-1] < W_max ){
    // Add upper limit
    W_limits.push_back( W_max ) ; 
    weights.push_back( weight ) ; 
  } 

  if( W_limits[0] > W_min ) {
    // Add lower limit from spline:
    W_limits.insert( W_limits.begin(), W_min ) ; 
    weights.insert( fWeights.begin(), weight ) ; 
  }
                                                
} 
//_________________________________________________________________________
    
double MECScaleVsW::ScaleFunction( double W, double W_min, double W_max, double scale_min, double scale_max ) const 
{
  // This function is responsible to calculate the scale at a given W
  // It interpolates the value between scale_min (W_min) and scale_max (W_max) linearly 

  double scale = ( W_min * scale_max - W_max * scale_min ) + ( scale_min - scale_max ) * W ;
  scale /= W_min - W_max ; 

  return scale ; 

} 

//_________________________________________________________________________

void MECScaleVsW::LoadConfig(void)
{
  bool good_config = true ; 
  GetParam( "MECScaleVsW-Weight-Default", fDefaultWeight ) ;
  GetParamVect( "MECScaleVsW-Weights", fWeights ) ; 
  GetParamVect( "MECScaleVsW-WValues", fWValues ) ; 

  if( fWeights.size() != fWValues.size() ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Entries don't match" ;
    LOG("MECScaleVsW", pERROR) << "Weights size: " << fWeights.size() ;
    LOG("MECScaleVsW", pERROR) << "WValues size: " << fWValues.size() ;
  }

  std::vector<double> limit_Q0, limit_Q3 ;
  if( GetConfig().Exists("MECScaleVsW-LowerLimitQ0") ) {
    GetParamVect( "MECScaleVsW-LowerLimitQ0", limit_Q0 ) ;
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Lower limit of the phase space for Q0 is not specified" ;
  }

  if( GetConfig().Exists("MECScaleVsW-LowerLimitQ3") ) {
    GetParamVect( "MECScaleVsW-LowerLimitQ0", limit_Q3 ) ;
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Lower limit of the phase space for Q3 is not specified" ;
  }

  if( limit_Q0.size() != limit_Q3.size() ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Entries don't match" ;
    LOG("MECScaleVsW", pERROR) << "Lower limit for Q0 size: " << limit_Q0.size() ;
    LOG("MECScaleVsW", pERROR) << "Lower limit for Q3 size: " << limit_Q3.size() ;
  }

  fW1_Q0Q3_limits = new TSpline3("fW1_Q0Q3_limits",&limit_Q0[0],&limit_Q3[0],limit_Q0.size());
  
  if( ! good_config ) {
    LOG("MECScaleVsW", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

}

//_________________________________________________________________________
