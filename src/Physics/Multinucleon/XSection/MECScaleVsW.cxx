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

  return 1 ; 
}
//_________________________________________________________________________
double MECScaleVsW::GetScaling( const double q0, const double q3 ) const
{
  // Get the vectors that include the kinematic limits of W for a given event:
  std::vector<double> W_limits, weights ; 
  GetVectorWithLimits( W_limits, weights,q0, q3, fDefaultWeight ) ;

  // The Scaling is done using the "experimenter's W", which assumes a single nucleon
  // See motivation in : https://arxiv.org/pdf/1601.02038.pdf
  // Calculate event W:
  double Mn = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  double W = sqrt( pow(Mn,2) + 2*Mn*q0 - pow(q0,2) + pow(q0,2) ) ;

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
				       const double q0, const double q3 , const double weight ) const {

  double Mn = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  double W_min = sqrt( pow(Mn,2) + 2*Mn*fW1_q0q3_limits->Eval(q3) - pow(q3,2) + pow(fW1_q0q3_limits->Eval(q3),2) ) ;
  double W_max = sqrt( pow(Mn,2) + 2*Mn*q0 ) ; // Imposing Q2 = 0 

  W_limits = fWValues ; 
  weights = fWeights ; 

  if( W_limits[0] > W_min ) {
    // Add lower limit from spline:
    W_limits.insert( W_limits.begin(), W_min ) ; 
    weights.insert( fWeights.begin(), weight ) ; 
  }

  if( W_limits[W_limits.size()-1] < W_max ){
    // Add upper limit
    W_limits.push_back( W_max ) ; 
    weights.push_back( weight ) ; 
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
  GetParamDef( "MECScaleVsW-Weight-Default", fDefaultWeight, 1. ) ;
  GetParamVect( "MECScaleVsW-Weights", fWeights ) ; 
  GetParamVect( "MECScaleVsW-WValues", fWValues ) ; 

  if( fWeights.size() != fWValues.size() ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Entries don't match" ;
    LOG("MECScaleVsW", pERROR) << "Weights size: " << fWeights.size() ;
    LOG("MECScaleVsW", pERROR) << "WValues size: " << fWValues.size() ;
  }

  std::vector<double> limit_q0, limit_q3 ;
  GetParamVect( "MECScaleVsW-LowerLimitQ0", limit_q0 ) ;
  GetParamVect( "MECScaleVsW-LowerLimitQ3", limit_q3 ) ;

  if( limit_q0.size() != limit_q3.size() ) {
    good_config = false ; 
    LOG("MECScaleVsW", pERROR) << "Entries don't match" ;
    LOG("MECScaleVsW", pERROR) << "Lower limit for q0 size: " << limit_q0.size() ;
    LOG("MECScaleVsW", pERROR) << "Lower limit for q3 size: " << limit_q3.size() ;
  }

  fW1_q0q3_limits = new TSpline3("fW1_q0q3_limits",&limit_q0[0],&limit_q3[0],limit_q0.size());
  
  /*
  // Store spline : ( include this in the xml? ) 
  // Digitalized from a high stat sample. This is Enu independent.
  double limit_q0[] = { 0.0479, 0.0288, 0.0192, 0.0288, 0.029, 0.0415, 0.0511, 0.0735, 0.1054, 0.1374, 0.1821, 0.2300, 0.2843, 0.3546, 0.4313, 0.4824};
  double limit_q3[] = { 0.054, 0.0963, 0.2139, 0.299, 0.3719, 0.4451, 0.5125, 0.5973, 0.6802, 0.7630, 0.8401, 0.9114, 1.0058, 1.0790, 1.1676, 1.1946};
  */

  /*
  fWValues.assign( 5 , 0. ) ;
  // fWValues[0] corresponds to W1. This is obtained from fW1_q0q3_limits 
  fWValues[1] = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  fWValues[2] = 1.12; // GeV. Value from https://arxiv.org/pdf/1601.02038.pdf // W_dip 
  fWValues[3] = PDGLibrary::Instance()->Find(kPdgP33m1232_DeltaP)->Mass(); // This corresponds to the Delta mass. 
  // fWValues[4] satisfies Q2 = 0. 
  */



}

//_________________________________________________________________________
