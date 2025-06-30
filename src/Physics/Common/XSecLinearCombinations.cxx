//_________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Original code contributed by J. Schwehr, D. Cherdack, R. Gran
 Substantial code refactorizations by the core GENIE group.
*/
//_________________________________________________________________________

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/Common/XSecLinearCombinations.h"

using namespace genie;

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
  // Store xsec prediction for each entry :
  for( unsigned int i = 0 ; i < fXSections.size() ; ++ i ) {
    if( fXSections[i] -> ValidKinematics( interaction ) ) { 
      xsec += ( fXSections[i] -> XSec( interaction, kps ) ) * fLinearCoefficients[i] ;
    }
  }

  return xsec;
}

//_________________________________________________________________________

double XSecLinearCombinations::Integral( const Interaction * interaction ) const
{

  double xsec = 0 ; 
  // Store xsec prediction for each entry :
  for( unsigned int i = 0 ; i < fXSections.size() ; ++ i ) {
    xsec += ( fXSections[i] -> Integral( interaction ) ) * fLinearCoefficients[i] ; 
  }

  return xsec;
}

//_________________________________________________________________________

bool XSecLinearCombinations::ValidProcess( const Interaction * interaction) const
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
  bool good_config = true ;
  GetParamVect( "LinearCoefficients", fLinearCoefficients ) ; 

  bool Normalise = false ;
  GetParam( "Normalise", Normalise ) ;

  if( Normalise ) {
    double sum = 0 ; 
    for( unsigned int i = 0 ; i < fLinearCoefficients.size() ; ++i ) {
      sum += fLinearCoefficients[i] ;
    }

    if( sum == 0 ) {
      good_config = false ;
      LOG("XSecLinearCombinations", pERROR) << "Sum of linear coefficients is zero. Cannot normalize " ;
    }

    for( unsigned int j = 0 ; j < fLinearCoefficients.size() ; ++j ) {
      fLinearCoefficients[j] /= sum ;
    }
  }

  std::vector<RgKey> XSectionsKeys ;
  GetParamVectKeys( "CrossSection", XSectionsKeys ) ;

  if( XSectionsKeys.size() != fLinearCoefficients.size() ) {
    good_config = false ; 
    LOG("XSecLinearCombinations", pERROR) << "Entries don't match" ;
    LOG("XSecLinearCombinations", pERROR) << "Cross Sections size: " << XSectionsKeys.size() ;
    LOG("XSecLinearCombinations", pERROR) << "Linear Coefficients size: " << fLinearCoefficients.size() ;
  }

  for( unsigned int i = 0 ; i < XSectionsKeys.size() ; ++i ) {
    fXSections.push_back( dynamic_cast<const XSecAlgorithmI *> ( this->SubAlg( XSectionsKeys[i] ) ) );
    if ( ! fXSections[i] ) {
      good_config = false ;
      LOG("XSecLinearCombinations", pERROR) << "The SubAlg is not a XSecAlgorithmI";
    }
  }

  if( ! good_config ) {
    LOG("XSecLinearCombinations", pERROR) << "XSecLinearCombilation Configuration has failed.";
    exit(78) ;
  }

}
