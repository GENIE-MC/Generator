//_________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J. Schwehr, D. Cherdack, R. Gran
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/LinearCombination/XSecLinearCombinations.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"

using namespace genie;
using namespace genie::constants;

//_________________________________________________________________________
XSecLinearCombinations::XSecLinearCombinations() :
  XSecAlgorithmI("genie::XSecLinearCombinations")
{

}
//_________________________________________________________________________
XSecLinearCombinations::XSecLinearCombinations(string config) :
  XSecAlgorithmI("genie::XSecLinearCombinations", config)
{

}
//_________________________________________________________________________
XSecLinearCombinations::~XSecLinearCombinations()
{

}
//_________________________________________________________________________
double XSecLinearCombinations::XSec(
  const Interaction * interaction, KinePhaseSpace_t kps) const
{
  double xsec = 0 ; 
  std::vector<double> vect_xsec ;

  // Store xsec prediction for each entry :
  for( unsigned int i = 0 ; i < fXSections.size() ; ++ i ) {
    vect_xsec.push_back( fXSections[i] -> XSec( interaction, kps ) ) ;
  }

  // Check if we need to normalize the weights: 
  double total_w = 1 ; 
  if( fNormalise ) {
    double sum = 0 ; 
    for( unsigned int j = 0 ; j < fLinearCoefficients.size() ; ++j ) {
      sum += fLinearCoefficients[j] ;
    }
    total_w = sum ; 
  }

  for( unsigned int k = 0; k < vect_xsec.size() ; ++k ) {
    xsec += vect_xsec[k] * fLinearCoefficients[k] ; 
    xsec *= total_w ; 
  }  

  return xsec;
}
//_________________________________________________________________________
double XSecLinearCombinations::Integral( const Interaction * interaction ) const
{

  double xsec = 0 ; 
  std::vector<double> vect_xsec ;

  // Store xsec prediction for each entry :
  for( unsigned int i = 0 ; i < fXSections.size() ; ++ i ) {
    vect_xsec.push_back( fXSections[i] -> Integral( interaction ) ) ; 
  }

  // Check if we need to normalize the weights: 
  double total_w = 1 ; 
  if( fNormalise ) {
    double sum = 0 ; 
    for( unsigned int j = 0 ; j < fLinearCoefficients.size() ; ++j ) {
      sum += fLinearCoefficients[j] ;
    }
    total_w = sum ; 
  }

  for( unsigned int k = 0; k < vect_xsec.size() ; ++k ) {
    xsec += vect_xsec[k] * fLinearCoefficients[k] ; 
    xsec *= total_w ; 
  }  
  
  return xsec;
}

//_________________________________________________________________________

bool XSecLinearCombinations::ValidProcess(
        const Interaction * interaction) const
{
  bool is_valid = true ;

  // CHECK THAT THE XSEC USED CONSIDER A VALID PROCESS
  for( unsigned int i = 0 ; i < fXSections.size() ; ++ i ) {
    if( ! fXSections[i] -> ValidProcess( interaction ) ) {
      is_valid = false ;
      break ; 
    }
  }
  
  return is_valid ;
}

//_________________________________________________________________________
void XSecLinearCombinations::Configure(const Registry & config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//____________________________________________________________________________
void XSecLinearCombinations::Configure(string config)
{
    Algorithm::Configure(config);
    this->LoadConfig();
}
//_________________________________________________________________________
void XSecLinearCombinations::LoadConfig(void)
{
  GetParamVect( "LinearCoefficients", fLinearCoefficients ) ; 
  GetParam( "Normalise", fNormalise ) ;

  GetParamVectKeys( "CrossSection", fXSectionsKeys ) ;

  if( fXSectionsKeys.size() != fLinearCoefficients.size() ) {
    LOG("XSecLinearCombinations", pFATAL) << "Entries don't match" ;
    LOG("XSecLinearCombinations", pFATAL) << "Cross Sections size: " << fXSectionsKeys.size() ;
    LOG("XSecLinearCombinations", pFATAL) << "Linear Coefficients size: " << fLinearCoefficients.size() ;
  }

  for( unsigned int i = 0 ; i < fXSectionsKeys.size() ; ++i ) {
    fXSections.push_back( dynamic_cast<const XSecAlgorithmI *> ( this->SubAlg( fXSectionsKeys[i] ) ) );
    assert( fXSections[i] ) ;
  }

}
