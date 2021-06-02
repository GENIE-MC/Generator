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
#include "Physics/Multinucleon/XSection/ScaleMECW.h"
#include "Framework/Utils/StringUtils.h" 
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;

//_________________________________________________________________________
ScaleMECW::ScaleMECW() : 
  Algorithm("genie::ScaleMECW") 
{
  
}
//_________________________________________________________________________
ScaleMECW::ScaleMECW(string config) : 
  Algorithm("genie::ScaleMECW",config) 
{
  
}
//_________________________________________________________________________
ScaleMECW::~ScaleMECW()
{

}
//_________________________________________________________________________
void ScaleMECW::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void ScaleMECW::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________
double ScaleMECW::GetScaling( const double q0, const double q3 ) const
{
  // The Scaling is done using the "experimenter's W", which assumes a single nucleon
  // See motivation in : https://arxiv.org/pdf/1601.02038.pdf
  double Mn = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass

  std::vector<double> limits = fLimits ; 
  // The W1 and W2 limits vary with Q2. They must be updated:
  limits[0] = sqrt( pow(Mn,2) + 2*Mn*fW1_q0q3_limits->Eval(q3) - pow(q3,2) + pow(fW1_q0q3_limits->Eval(q3),2) ) ;
  limits[4] = sqrt( pow(Mn,2) + 2*Mn*q0 ) ; // Imposing Q2 = 0
   
  // Calculate event W:
  double W = sqrt( pow(Mn,2) + 2*Mn*q0 - pow(q0,2) + pow(q0,2) ) ;

  // Calculate scaling:

  for ( unsigned int i = 0 ; i < fWeights.size() - 1 ; ++i ) {
    if ( W >= limits[i] && W < limits[i+1] ) {
      return ScaleFunction( W, limits[i], limits[i+1], fWeights[i], fWeights[i+1] ) ; 
    } 
  }

  return 1. ; 
} 

//_________________________________________________________________________
double ScaleMECW::ScaleFunction( double W, double W_min, double W_max, double scale_min, double scale_max ) const 
{
  // This function is responsible to calculate the scale at a given W
  // It interpolates the value between scale_min (W_min) and scale_max (W_max) linearly 

  double scale = ( W_min * scale_max - W_max * scale_min ) + ( scale_min - scale_max ) * W ;
  scale /= W_min - W_max ; 

  return scale ; 

} 

//_________________________________________________________________________

void ScaleMECW::LoadConfig(void)
{
  double default_weight ; 
  GetParamDef( "ScaleMECW-Weight-Default", default_weight, 1. ) ;
  fWeights.assign( 5 , default_weight ) ;

  // Store weights in vector 
  GetParam( "ScaleMECW-Weight-QELRegion", fWeights[1] ) ;
  GetParam( "ScaleMECW-Weight-RESRegion", fWeights[3] ) ;

  fLimits.assign( 5 , 0. ) ;
  // fLimits[0] corresponds to W1. This is obtained from fW1_q0q3_limits 
  fLimits[1] = ( PDGLibrary::Instance()->Find(kPdgProton)->Mass() + PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() ) * 0.5 ;  // Nucleon mass
  fLimits[2] = 1.12; // GeV. Value from https://arxiv.org/pdf/1601.02038.pdf // W_dip 
  fLimits[3] = PDGLibrary::Instance()->Find(kPdgP33m1232_DeltaP)->Mass(); // This corresponds to the Delta mass. 
  // fLimits[4] satisfies Q2 = 0. 

  // Store spline : ( include this in the xml? ) 
  // Digitalized from a high stat sample. This is Enu independent.
  double q0_min[] = { 0.0479, 0.0288, 0.0192, 0.0288, 0.029, 0.0415, 0.0511, 0.0735, 0.1054, 0.1374, 0.1821, 0.2300, 0.2843, 0.3546, 0.4313, 0.4824};
  double q3_min[] = { 0.054, 0.0963, 0.2139, 0.299, 0.3719, 0.4451, 0.5125, 0.5973, 0.6802, 0.7630, 0.8401, 0.9114, 1.0058, 1.0790, 1.1676, 1.1946};

  fW1_q0q3_limits = new TSpline3("fW1_q0q3_limits",q3_min,q0_min,16);

}

//_________________________________________________________________________
