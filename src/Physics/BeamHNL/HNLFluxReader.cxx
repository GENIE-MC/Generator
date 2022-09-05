//----------------------------------------------------------------------------
// Implementation file for HNLFluxReader.h
// 10/Jan/2021
//----------------------------------------------------------------------------
// TODO: figure out naming scheme for flux files
//       exceptions + error messages
//       what to do about K0L, mu parents & polarisation??
//       figure out implementation of generateVtxT
//----------------------------------------------------------------------------

#include "HNLFluxReader.h"

using namespace genie::HNL;
using namespace genie::HNL::HNLenums;

// definition of extern variables
std::string HNLFluxReader::fPath;

int HNLFluxReader::selectMass( const double mN ){
    /// move mass to closest mass hypothesis

    const massHyp_t massesHyp[] = {
	kLight0Hyp,  kLight1Hyp,  kLight2Hyp,  kLight3Hyp,  kLight4Hyp,
	kLight5Hyp,  kLight6Hyp,  kLight7Hyp,  kLight8Hyp,  kLight9Hyp,
	kLightAHyp,  kLightBHyp,  kLightCHyp,  kLightDHyp,  kLightEHyp,
	kLightFHyp,  kLightGHyp,  kLightHHyp,  kLightIHyp,  kLightJHyp,
	kMedium0Hyp, kMedium1Hyp, kMedium2Hyp, kMedium3Hyp, kMedium4Hyp,
	kMedium5Hyp, kMedium6Hyp, kMedium7Hyp, kMedium8Hyp, kMedium9Hyp,
	kHeavy0Hyp,  kHeavy1Hyp,  kHeavy2Hyp,  kHeavy3Hyp,  kHeavy4Hyp,
	kHeavy5Hyp,  kHeavy6Hyp,  kHeavy7Hyp,  kHeavy8Hyp,  kHeavy9Hyp,
	kHeavyAHyp,  kHeavyBHyp,  kHeavyCHyp,  kHeavyDHyp,  kHeavyEHyp,
	kHeavyFHyp,  kHeavyGHyp,  kHeavyHHyp,  kHeavyIHyp,  kHeavyJHyp };
    const int nMasses = sizeof(massesHyp)/sizeof(massesHyp[0]) - 1;

    // because masses are kept in a *map*, gotta build array of the second elements!
    double masses[ nMasses + 1 ];
    for( int i = 0; i <= nMasses; i++ ){
      massHyp_t thisHyp = massesHyp[i];
      auto pos = massHypMap.find( thisHyp );
      masses[i] = pos->second;
      LOG("HNL", pINFO) 
	<< "At position " << i << " the mass hypothesis is " << masses[i];
    }

    if( mN < 0.0 ){
      LOG("HNL", pERROR) 
	<< "Illegal mass mN = " << mN; 
      exit(1); 
    }

    int mp = -1; double fmN = 0.0; // mass of the selected masspoint
    if( mN >= masses[ nMasses ] ){ mp = nMasses; fmN = masses[ nMasses ]; }
    while( masses[ mp + 1 ] < mN && mp < nMasses ){ mp++; } // decide interval
    
    // generally decide mass + point by closest endpoint in interval
    const double dLeft  = std::abs( mN - masses[ mp ] );
    const double dRight = masses[ mp + 1 ] - mN;

    LOG("HNL", pDEBUG)
     << "Stats:"
     << "\n Input mass: " << mN 
     << "\n Choice interval: [ " << masses[mp] << ", " << masses[mp+1] << " ] " 
     << "\n Left and right distance: " << dLeft << ", " << dRight
     << "\n Chosen point: " << ( ( dLeft < dRight ) ? mp : mp + 1 ) << " ( "
     << ( ( dLeft < dRight ) ? "LEFT )" : "RIGHT )" );

    fmN = ( dLeft < dRight ) ? masses[ mp ] : masses[ mp + 1 ];
    mp  = ( dLeft < dRight ) ? mp : mp + 1;

    return mp;
    
}

void HNLFluxReader::selectFile( std::string fin, const double mN ){
  std::string filePath = fin;
  filePath.append( "/FHC/" ); // RETHERE make this configurable!!!
  filePath.append( "EqualCouplings/" );
  const int mp = selectMass( mN );
  filePath.append( Form( "mp%02d", mp ) );
  filePath.append( ".root" );
  fPath = filePath;
}

TH1F * HNLFluxReader::getFluxHist1F( std::string fin, std::string hName, int HType ){
    // flux file contains 4 dirs: numu, numubar, nue, nuebar. Each has flux + helper hists
    TFile *f = TFile::Open( fin.c_str() );

    // descend into dir first!
    TDirectory *baseDir = f->GetDirectory( "" );
    std::string strType;
    switch( HType ){
	case kNumu:    strType = std::string( "numu" ); break;
	case kNumubar: strType = std::string( "numubar" ); break;
	case kNue:     strType = std::string( "nue" ); break;
	case kNuebar:  strType = std::string( "nuebar" ); break;
    }
    TDirectory *deepDir = baseDir->GetDirectory( strType.c_str( ) );
    assert( deepDir );

    LOG( "HNL", pDEBUG )
      << "Getting flux from histo with base name " << hName.c_str();

    // construct standardised name
    std::string strHist = std::string( hName.c_str( ) );
    strHist.append( "_all" );
    assert( deepDir->GetListOfKeys()->Contains( strHist.c_str() ) );
    
    TH1F * histPtr = dynamic_cast< TH1F* >( deepDir->Get( strHist.c_str() ) );

    return histPtr;
}

// overload to new interface
TH1F * HNLFluxReader::getFluxHist1F( std::string fin, int masspoint, bool isParticle ){
  
  TFile *f = TFile::Open( fin.c_str() );
  assert( f );
  
  std::string histName( Form("hHNLFluxCenterAcc_%d", masspoint) );
  if( !isParticle ) histName.append("_bar");
  TObject * histObj = f->Get( histName.c_str() );
  
  LOG( "HNL", pDEBUG )
    << "Getting flux from histo with name " << histObj->GetName()
    << " and title " << histObj->GetTitle();

  TH1F * histPtr = dynamic_cast< TH1F* >( histObj );

  return histPtr;
}

TH3D * HNLFluxReader::getFluxHist3D( std::string fin, std::string dirName, std::string hName ){
  TFile * f = TFile::Open( fin.c_str() );

  TDirectory * baseDir = f->GetDirectory( "" );
  TDirectory * deepDir = baseDir->GetDirectory( dirName.c_str() );
  assert( deepDir );
  assert( deepDir->GetListOfKeys()->Contains( hName.c_str() ) );

  TH3D * histPtr = dynamic_cast< TH3D* >( deepDir->Get( hName.c_str() ) );
  return histPtr;
}

std::vector< double > * HNLFluxReader::generateVtx3X( TH3D * prodVtxHist )
{
  if( prodVtxHist->GetEntries() == 0 ){ // no production vtx? ok, set HNL to be 1km upstream
    std::vector< double > * vtxDir = new std::vector< double >();
    vtxDir->emplace_back( 0.0 );
    vtxDir->emplace_back( 0.0 );
    vtxDir->emplace_back( -100000.0 );
    vtxDir->emplace_back( 0.0 );
    return vtxDir;
  }

  double ux = 0.0, uy = 0.0, uz = 0.0;
  prodVtxHist->GetRandom3( ux, uy, uz );

  std::vector< double > * vtxDir = new std::vector< double >();
  vtxDir->emplace_back( ux );
  vtxDir->emplace_back( uy );
  vtxDir->emplace_back( uz );
  vtxDir->emplace_back( 0.0 ); // RETHERE: add some timing information?

  return vtxDir;
}
