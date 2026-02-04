//_________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J.Tena and M.Roda
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Physics/Multinucleon/XSection/MECScaleVsW.h"
#include "Framework/Utils/StringUtils.h" 
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;

//_________________________________________________________________________
MECScaleVsW::MECScaleVsW() : 
  XSecScaleI("genie::MECScaleVsW")
{
  
}
//_________________________________________________________________________
MECScaleVsW::MECScaleVsW(string config) : 
  XSecScaleI("genie::MECScaleVsW",config) 
{
  
}
//_________________________________________________________________________
MECScaleVsW::~MECScaleVsW()
{

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
  MECScaleVsW::weight_type_map weight_map = GetMapWithLimits( Q0, Q3 ) ;
    
  // The Scaling is done using the "experimenter's W", which assumes a single nucleon
  // See motivation in : https://arxiv.org/pdf/1601.02038.pdf
  // Calculate event W:
  static double Mn = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  double W = pow(Mn,2) + 2*Mn*Q0 - pow(Q3,2) + pow(Q0,2) ;
  // Do not scale if W<0. This can happen while we try to get the correct kinematics.
  // If the kinematics is not correct, W can be negative, and we scale with a nan.
  // To avoid this we do this check.
  if ( W < 0 ) return 1. ; 
  W = sqrt( W ) ; 

  // Calculate scaling:
  MECScaleVsW::weight_type_map::iterator it_min = weight_map.begin() ; 
  MECScaleVsW::weight_type_map::iterator it_max = std::next( weight_map.begin(), weight_map.size() -1 ) ; 

  if ( W < it_min->first || W > it_max->first ) return fDefaultWeight ; 

  while ( std::distance( it_min, it_max ) > 1 ) {
    unsigned int step = std::distance( weight_map.begin(), it_min ) + std::distance( it_min, it_max ) / 2 ; 
    MECScaleVsW::weight_type_map::iterator it_middle = std::next( weight_map.begin(), step ) ; 
    if ( W < it_middle->first ) it_max = it_middle ; 
    else it_min = it_middle ; 
  }

  return ScaleFunction( W, *it_min, *it_max ) ; 
} 

//_________________________________________________________________________
MECScaleVsW::weight_type_map MECScaleVsW::GetMapWithLimits( const double Q0, const double Q3 ) const {
  // This function is responsible to add the phase space limits in the WValues vector in case they are not included
  // in the configuration setup.

  static double Mn = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  double W_min = sqrt( pow(Mn,2) + 2*Mn*fW1_Q0Q3_limits.Eval(Q3) - pow(Q3,2) + pow(fW1_Q0Q3_limits.Eval(Q3),2) ) ;
  double W_max = sqrt( pow(Mn,2) + 2*Mn*Q0 ) ; // Imposing Q2 = 0 

  // Insert phase space limits:
  MECScaleVsW::weight_type_map w_map = fWeightsMap ; 
  w_map.insert( weight_type_pair( W_max, fUpperLimitWeight ) ) ;
  w_map.insert( weight_type_pair( W_min, fLowLimitWeight ) ) ;

  return w_map ; 
} 
//_________________________________________________________________________
    
double MECScaleVsW::ScaleFunction( const double W, const weight_type_pair min, const weight_type_pair max ) const 
{
  // This function is responsible to calculate the scale at a given W
  // It interpolates the value between scale_min (W_min) and scale_max (W_max) linearly 
  return ( max.second - min.second ) * ( W - min.first ) / ( max.first - min.first ) + min.second ; 
} 

//_________________________________________________________________________

void MECScaleVsW::LoadConfig(void)
{
  bool good_config = true ; 
  // Reset members
  fWeightsMap.clear(); 

  if( GetConfig().Exists("MECScaleVsW-Default-Weight") ) {
    GetParam( "MECScaleVsW-Default-Weight", fDefaultWeight ) ;
  } else {
      good_config = false ; 
      LOG("MECScaleVsW", pERROR) << "Default weight is not specified." ;
  }

  std::vector<double> Weights, WValues ;
  GetParamVect( "MECScaleVsW-Weights", Weights, false ) ; 
  GetParamVect( "MECScaleVsW-WValues", WValues, false ) ; 
  
  if( Weights.size() != WValues.size() ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Entries don't match" ;
    LOG("MECScaleVsW", pERROR) << "Weights size: " << Weights.size() ;
    LOG("MECScaleVsW", pERROR) << "WValues size: " << WValues.size() ;
  }

  // Store weights and WValues in map:
  for( unsigned int i = 0 ; i<Weights.size() ; ++i ) {
    fWeightsMap.insert( weight_type_pair( WValues[i], Weights[i] ) ) ;
  }
  
  std::vector<double> limit_Q0, limit_Q3 ;
  if( GetParamVect( "MECScaleVsW-LowerLimitQ0", limit_Q0 ) == 0 ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "MECScaleVsW-LowerLimitQ0 is empty" ;
  }

  if( GetParamVect( "MECScaleVsW-LowerLimitQ3", limit_Q3 ) == 0 ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "MECScaleVsW-LowerLimitQ3 is empty" ;
  }
  
  if( limit_Q0.size() != limit_Q3.size() ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Entries don't match" ;
    LOG("MECScaleVsW", pERROR) << "Lower limit for Q0 size: " << limit_Q0.size() ;
    LOG("MECScaleVsW", pERROR) << "Lower limit for Q3 size: " << limit_Q3.size() ;
  }

  GetParamDef("MECScleVsW-LowerLimit-Weight", fLowLimitWeight, fDefaultWeight ) ; 
  GetParamDef("MECScleVsW-UpperLimit-Weight", fUpperLimitWeight, fDefaultWeight ) ; 

  if( ! good_config ) {
    LOG("MECScaleVsW", pERROR) << "Configuration has failed.";
    exit(78) ;
  }

  fW1_Q0Q3_limits = TSpline3("fW1_Q0Q3_limits",limit_Q3.data(),limit_Q0.data(),limit_Q3.size()); 
  
}

//_________________________________________________________________________
