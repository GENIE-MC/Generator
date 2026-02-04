
//----------------------------------------------------------------------------
/*!
  Copyright (c) 2003-2025, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  
  Author: John Plows <komninos-john.plows \at physics.ox.ac.uk>
          University of Oxford
 */
//----------------------------------------------------------------------------

#include "Physics/BeamHNL/HNLBRCalculator.h"

using namespace genie;
using namespace genie::hnl;

//----------------------------------------------------------------------------
BRCalculator::BRCalculator() :
  ChannelCalculatorI("genie::hnl::BRCalculator")
{

}
//----------------------------------------------------------------------------
BRCalculator::BRCalculator(string name) :
 ChannelCalculatorI(name)
{

}
//----------------------------------------------------------------------------
BRCalculator::BRCalculator(string name, string config) :
  ChannelCalculatorI(name, config)
{

}
//----------------------------------------------------------------------------
BRCalculator::~BRCalculator()
{

}
//----------------------------------------------------------------------------
void BRCalculator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//----------------------------------------------------------------------------
void BRCalculator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//----------------------------------------------------------------------------
void BRCalculator::LoadConfig(void)
{
  if( fIsConfigLoaded ) return;

  this->GetParam( "HNL-Mass", fMass );
  this->GetParamVect( "HNL-LeptonMixing", fCouplings );
  fUe42 = fCouplings.at(0);
  fUm42 = fCouplings.at(1);
  fUt42 = fCouplings.at(2);
  this->GetParam( "HNL-Majorana", fMajorana );

  mPi0 = PDGLibrary::Instance()->Find(genie::kPdgPi0)->Mass();     
  mPi  = PDGLibrary::Instance()->Find(genie::kPdgPiP)->Mass();     
  mMu  = PDGLibrary::Instance()->Find(genie::kPdgMuon)->Mass();    
  mK   = PDGLibrary::Instance()->Find(genie::kPdgKP)->Mass();	     
  mK0  = PDGLibrary::Instance()->Find(genie::kPdgK0)->Mass();	     
  mE   = PDGLibrary::Instance()->Find(genie::kPdgElectron)->Mass();

  this->GetParam( "WeinbergAngle", wAng );
  s2w = std::pow( std::sin( wAng ), 2.0 );

  this->GetParam( "CKM-Vud", Vud );
  Vud2 = Vud * Vud;

  this->GetParam( "Pion-FFactor", fpi );
  fpi2 = fpi * fpi;

  BR_C1 = 1./4. * ( 1. - 4. * s2w + 8. * s2w * s2w );
  BR_C2 = 1./2. * ( -s2w + 2. * s2w * s2w );

  this->GetParam( "PMNS-Ue1", Ue1 );
  this->GetParam( "PMNS-Ue2", Ue2 );
  this->GetParam( "PMNS-Ue3", Ue3 );
  this->GetParam( "PMNS-Um1", Um1 );
  this->GetParam( "PMNS-Um2", Um2 );
  this->GetParam( "PMNS-Um3", Um3 );
  this->GetParam( "PMNS-Ut1", Ut1 );
  this->GetParam( "PMNS-Ut2", Ut2 );
  this->GetParam( "PMNS-Ut3", Ut3 );

  kscale_K3e = { 
    { 0.0, 1.0 }, { 0.01, (2.0 + 0.968309)/3.0 },
    { 0.019970, 0.968309 }, { 0.029963, 0.952842 }, { 0.040037, 0.922646 }, { 0.049839, 0.907908 },
    { 0.060075, 0.879136 }, { 0.070117, 0.824297 }, { 0.080272, 0.772880 }, { 0.090139, 0.736432 },
    { 0.099924, 0.679466 }, { 0.110059, 0.607039 }, { 0.120445, 0.560082 }, { 0.130123, 0.508503 },
    { 0.140579, 0.461674 }, { 0.150900, 0.405874 }, { 0.159906, 0.356819 }, { 0.170544, 0.303751 },
    { 0.180722, 0.267038 }, { 0.190278, 0.227323 }, { 0.200340, 0.193514 }, { 0.209579, 0.159513 },
    { 0.219244, 0.129386 }, { 0.230837, 0.101623 }, { 0.239932, 0.081113 }, { 0.249386, 0.060704 },
    { 0.260887, 0.043990 }, { 0.269425, 0.031878 }, { 0.280041, 0.021660 }, { 0.289206, 0.014251 },
    { 0.300601, 0.007729 }, { 0.310440, 0.004059 }, { 0.320600, 0.001935 }, { 0.330028, 0.000713 },
    { 0.339733, 0.000184 }, { 0.350852, 0.00000953 }
  };

  kscale_K3mu = {
    { 0.0, 1.0 }, { 0.01, (2.0 + 0.968309)/3.0 },
    { 0.019970, 0.968309 }, { 0.029963, 0.937622 }, { 0.040037, 0.907908 }, { 0.049839, 0.879136 },
    { 0.060075, 0.824297 }, { 0.070117, 0.772880 }, { 0.079757, 0.713094 }, { 0.090139, 0.647424 },
    { 0.100570, 0.578413 }, { 0.110059, 0.525146 }, { 0.120045, 0.454300 }, { 0.130964, 0.393012 },
    { 0.140127, 0.334561 }, { 0.149932, 0.275778 }, { 0.159906, 0.227323 }, { 0.170544, 0.181443 },
    { 0.180722, 0.137994 }, { 0.190278, 0.101623 }, { 0.200340, 0.071309 }, { 0.210933, 0.046917 },
    { 0.220661, 0.025445 }, { 0.229355, 0.013362 }, { 0.239932, 0.003930 }, { 0.250997, 0.000037 }
  };

  kscale_mu3e = {
    { 0.0, 1.0 }, { 0.01, (2.0 + 0.772880)/3.0 },
    { 0.020099, 0.772880 }, { 0.029963, 0.560082 }, { 0.040037, 0.356819 }, { 0.050161, 0.193514 },
    { 0.060075, 0.089341 }, { 0.069667, 0.032922 }, { 0.080272, 0.007247 }, { 0.090139, 0.000713 },
    { 0.100570, 0.00000363 }
  };

  fIsConfigLoaded = true;
}
//----------------------------------------------------------------------------
double BRCalculator::KinematicScaling( HNLProd_t hnlprod ) const
{
  return this->KScale_Global( hnlprod, fMass );
}
//----------------------------------------------------------------------------
double BRCalculator::DecayWidth( HNLDecayMode_t hnldm ) const
{
  return this->DWidth_Global( hnldm, fMass );
}
//----------------------------------------------------------------------------
double BRCalculator::GetFormfactorF1( double x ) const {
  if( x < 0. || x > 0.5 ) { LOG( "HNL", pERROR ) << "BRCalculator::GetFormfactorF1:: Illegal x = " << x; exit( 3 ); }
  if( x == 0.5 ) return 0.;
  int i = x/selector::PARTWIDTH;
  if( x - i*selector::PARTWIDTH ==0 ) return selector::FormfactorF1[i];
  return 1./2. * ( selector::FormfactorF1[i] + selector::FormfactorF1[i+1] );
}
//----------------------------------------------------------------------------
double BRCalculator::GetFormfactorF2( double x ) const {
  if( x < 0. || x > 0.5 ) { LOG( "HNL", pERROR ) << "BRCalculator::GetFormfactorF2:: Illegal x = " << x; exit( 3 ); }
  if( x == 0.5 ) return 0.;
  int i = x/selector::PARTWIDTH;
  if( x - i*selector::PARTWIDTH==0 ) return selector::FormfactorF2[i];
  return 1./2. * ( selector::FormfactorF2[i] + selector::FormfactorF2[i+1] );
}
//----------------------------------------------------------------------------
// interface to scale factors
double BRCalculator::KScale_Global( HNLProd_t hnldm, const double M ) const {
  if( !utils::hnl::IsProdKinematicallyAllowed( hnldm ) ){
    return 0.0;
  }
  
  switch( hnldm ){
  case kHNLProdPion2Muon: return KScale_PseudoscalarToLepton( mPi, M, mMu );
  case kHNLProdPion2Electron: return KScale_PseudoscalarToLepton( mPi, M, mE );
  case kHNLProdKaon2Muon: return KScale_PseudoscalarToLepton( mK, M, mMu );
  case kHNLProdKaon2Electron: return KScale_PseudoscalarToLepton( mK, M, mE );
  case kHNLProdKaon3Muon: return KScale_PseudoscalarToPiLepton( mK, M, mMu );
  case kHNLProdKaon3Electron: return KScale_PseudoscalarToPiLepton( mK, M, mE );
  case kHNLProdNeuk3Muon: return KScale_PseudoscalarToPiLepton( mK0, M, mMu );
  case kHNLProdNeuk3Electron: return KScale_PseudoscalarToPiLepton( mK0, M, mE );
  case kHNLProdMuon3Numu:
  case kHNLProdMuon3Nue:
  case kHNLProdMuon3Nutau:
    return KScale_MuonToNuAndElectron( M );
  default: return 0.0;
  }

  return 0.0;
}
//----------------------------------------------------------------------------
// HNL production widths
double BRCalculator::KScale_PseudoscalarToLepton( const double mP, const double M, const double ma ) const {
  double da = std::pow( utils::hnl::MassX( ma, mP ) , 2.0 );
  double di = std::pow( utils::hnl::MassX( M,  mP ) , 2.0 );
  double num = utils::hnl::RhoFunc( da, di );
  double den = da * std::pow( (1.0 - da), 2.0 );
  return num/den;
}
//----------------------------------------------------------------------------
double BRCalculator::DWidth_PseudoscalarToLepton( const double mP, const double M, const double Ua42, const double ma ) const {
  assert( M + ma <= mP );

  double KScale = KScale_PseudoscalarToLepton( mP, M, ma );
  return Ua42 * KScale;
}
//----------------------------------------------------------------------------
double BRCalculator::KScale_PseudoscalarToPiLepton( const double mP, const double M, const double ma ) const {
  assert( mP == mK || mP == mK0 ); // RETHERE remove this when/if heavier pseudoscalars are considered
  assert( ma == mE || ma == mMu );
  
  std::map< double, double > scaleMap = ( ma == mE ) ? kscale_K3e : kscale_K3mu;

  std::map< double, double >::iterator scmit = scaleMap.begin();
  // iterate until we know between which two map points M is
  // if we're very lucky, M will coincide with a map point
  while( (*scmit).first <= M && scmit != scaleMap.end() ){ ++scmit; }
  std::map< double, double >::iterator scpit = std::prev( scmit, 1 );
  //LOG( "HNL", pDEBUG )
  //  << "Requested map for M = " << M << ": iter at ( " << (*scpit).first << ", " << (*scmit).first << " ]";
  assert( scmit != scaleMap.end() );
  // if coincide then return scale there
  if( scaleMap.find( M ) != scaleMap.end() ) return (*scmit).second;
  // otherwise transform scmit-1 and scmit second to log, do a linear extrapolation and return
  double l1 = TMath::Log( (*scpit).second );
  double l2 = TMath::Log( (*scmit).second );
  double t  = ( M - (*scpit).first ) / ( (*scmit).first - (*scpit).first );
  return TMath::Exp( l1 + ( l2 - l1 ) * t );
}
//----------------------------------------------------------------------------
double BRCalculator::DWidth_PseudoscalarToPiLepton( const double mP, const double M, const double Ua42, const double ma ) const {
  assert( M + ma + mPi0 <= mP );

  double KScale = KScale_PseudoscalarToPiLepton( mP, M, ma );
  return Ua42 * KScale;
}
//----------------------------------------------------------------------------
double BRCalculator::KScale_MuonToNuAndElectron( const double M ) const {
  std::map< double, double > scaleMap = kscale_mu3e;
  std::map< double, double >::iterator scmit = scaleMap.begin();
  while( (*scmit).first <= M && scmit != scaleMap.end() ){ ++scmit; }
  std::map< double, double >::iterator scpit = std::prev( scmit, 1 );
  //LOG( "HNL", pDEBUG )
  //  << "Requested map for M = " << M << ": iter at ( " << (*scpit).first << ", " << (*scmit).first << " ]";
  assert( scmit != scaleMap.end() );

  if( scaleMap.find( M ) != scaleMap.end() ) return (*scmit).second;

  double l1 = TMath::Log( (*scpit).second );
  double l2 = TMath::Log( (*scmit).second );
  double t  = ( M - (*scpit).first ) / ( (*scmit).first - (*scpit).first );
  return TMath::Exp( l1 + ( l2 - l1 ) * t );
}
//----------------------------------------------------------------------------
double BRCalculator::DWidth_MuonToNuAndElectron( const double M, const double Ue42, const double Umu42, const double Ut42 ) const {
  assert( M + mE <= mMu );

  double KScale = KScale_MuonToNuAndElectron( M );
  return ( Ue42 + Umu42 + Ut42 ) * KScale;
}
//----------------------------------------------------------------------------
// interface to decay widths
double BRCalculator::DWidth_Global( HNLDecayMode_t hnldm, const double M ) const {
  if( !utils::hnl::IsKinematicallyAllowed( hnldm, M ) ){
    return 0.0;
  }
  
  switch( hnldm ){
  case kHNLDcyNuNuNu   : return this->DWidth_Invisible( M, fUe42, fUm42, fUt42 );
  case kHNLDcyNuEE     : return this->DWidth_SameLepton( M, fUe42, fUm42, fUt42, mE, false );
  case kHNLDcyNuMuE    : return (this->DWidth_DiffLepton( M, fUe42, fUm42, fMajorana ) + this->DWidth_DiffLepton( M, fUm42, fUe42, fMajorana ));
  case kHNLDcyPi0Nu    : return this->DWidth_PiZeroAndNu( M, fUe42, fUm42, fUt42 );
  case kHNLDcyPiE      : return this->DWidth_PiAndLepton( M, fUe42, mE );
  case kHNLDcyNuMuMu   : return this->DWidth_SameLepton( M, fUe42, fUm42, fUt42, mMu, true );
  case kHNLDcyPiMu     : return this->DWidth_PiAndLepton( M, fUm42, mMu );
  case kHNLDcyPi0Pi0Nu : return this->DWidth_Pi0Pi0Nu( M, fUe42, fUm42, fUt42 );
  case kHNLDcyPiPi0E   : return this->DWidth_PiPi0Ell( M, mE, fUe42, fUm42, fUt42, true );
  case kHNLDcyPiPi0Mu  : return this->DWidth_PiPi0Ell( M, mMu, fUe42, fUm42, fUt42, false );
  default: return 0.0;
  }

  return 0.0;
}
//----------------------------------------------------------------------------
// total decay widths, various channels
double BRCalculator::DWidth_PiZeroAndNu( const double M, const double Ue42, const double Umu42, const double Ut42 ) const {
  const double x       = genie::utils::hnl::MassX( mPi0, M );
  const double preFac  = GF2 * M*M*M / ( 32. * pi );
  const double kinPart = ( 1. - x*x ) * ( 1. - x*x );
  return preFac * ( Ue42 + Umu42 + Ut42 ) * fpi2 * kinPart;
}
//----------------------------------------------------------------------------
double BRCalculator::DWidth_PiAndLepton( const double M, const double Ua42, const double ma ) const {
  const double xPi     = genie::utils::hnl::MassX( mPi, M );
  const double xLep    = genie::utils::hnl::MassX( ma, M );
  const double preFac  = GF2 * M*M*M / ( 16. * pi );
  const double kalPart = TMath::Sqrt( genie::utils::hnl::Kallen( 1, xPi*xPi, xLep*xLep ) );
  const double othPart = 1. - xPi*xPi - xLep*xLep * ( 2. + xPi*xPi - xLep*xLep );
  return preFac * fpi2 * Ua42 * Vud2 * kalPart * othPart;
}
//----------------------------------------------------------------------------
double BRCalculator::DWidth_Invisible( const double M, const double Ue42, const double Umu42, const double Ut42 ) const {
  const double preFac = GF2 * TMath::Power( M, 5. ) / ( 192. * pi*pi*pi );
  return preFac * ( Ue42 + Umu42 + Ut42 );
}
//----------------------------------------------------------------------------
double BRCalculator::DWidth_SameLepton( const double M, const double Ue42, const double Umu42, const double Ut42, const double mb, bool bIsMu ) const {
  const double preFac = GF2 * TMath::Power( M, 5. ) / ( 192. * pi*pi*pi );
  const double x      = genie::utils::hnl::MassX( mb, M );
  const double f1     = GetFormfactorF1( x );
  const double f2     = GetFormfactorF2( x );
  const double C1Part = ( Ue42 + Umu42 + Ut42 ) * f1 * BR_C1;
  const double C2Part = ( Ue42 + Umu42 + Ut42 ) * f2 * BR_C2;
  const double D1Part = bIsMu ? 2. * s2w * Umu42 * f1 : 2. * s2w * Ue42 * f1;
  const double D2Part = bIsMu ? s2w * Umu42 * f2 : s2w * Ue42 * f2;
  return preFac * ( C1Part + C2Part + D1Part + D2Part );
}
//----------------------------------------------------------------------------
double BRCalculator::DWidth_DiffLepton( const double M, const double Ua42, const double Ub42, const int IsMajorana ) const {
  const double preFac = GF2 * TMath::Power( M, 5. ) / ( 192. * pi*pi*pi );
  const double x = genie::utils::hnl::MassX( mMu, M );
  const double kinPol = 1. - 8. * x*x + 8. * TMath::Power( x, 6. ) - TMath::Power( x, 8. );
  const double kinLn  = -12. * TMath::Power( x, 4. ) * TMath::Log( x*x );
  const double kinPart = kinPol + kinLn;
  const double coupPart = (!IsMajorana) ? Ua42 : Ua42 + Ub42; // 2nd diagram in Majorana case!
  return preFac * kinPart * coupPart;
}
//----------------------------------------------------------------------------
// note that these BR are very very tiny.
double BRCalculator::DWidth_PiPi0Ell( const double M, const double ml,
					      const double Ue42, const double Umu42, const double Ut42,
					      const bool isElectron) const
{
  // because the actual decay width is very hard to integrate onto a full DWidth,
  // build 2Differential and then integrate numerically
  // using Simpson's method for 2D.

  const double preFac = fpi2 * fpi2 * GF2 * GF2 * Vud2 * M / ( 32.0 * pi*pi*pi );
  const double Ua1 = isElectron ? Ue1 : Um1;
  const double Ua2 = isElectron ? Ue2 : Um2;
  const double Ua3 = isElectron ? Ue3 : Um3;
  __attribute__((unused)) const double Ua4 = isElectron ? std::sqrt( Ue42 ) : std::sqrt( Umu42 );

  const double Ue4 = std::sqrt( Ue42 );
  const double Um4 = std::sqrt( Umu42 );
  const double Ut4 = std::sqrt( Ut42 );
  // assume all these to be real
  const double bigMats =
    Ua1 * ( Ue4 * Ue1 + Um4 * Um1 + Ut4 * Ut1 ) +
    Ua2 * ( Ue4 * Ue2 + Um4 * Um2 + Ut4 * Ut2 ) +
    Ua3 * ( Ue4 * Ue3 + Um4 * Um3 + Ut4 * Ut3 );

  // now limits
  const double maxMu =
    ( ( M - mPi0 ) * ( M - mPi0 ) - mPi*mPi + ml*ml ) / ( 2.0 * ( M - mPi0 ) );
  const double maxPi =
    ( ( M - mPi0 ) * ( M - mPi0 ) + mPi*mPi - ml*ml ) / ( 2.0 * ( M - mPi0 ) );

  // gotta put in the formula
  TF2 * f = new TF2( "fPiPi0Ell", PiPi0EllForm, mPi, maxPi, ml, maxMu, 4 );
  f->SetParameter( 0, M );
  f->SetParameter( 1, ml );
  f->SetParameter( 2, mPi );
  f->SetParameter( 3, mPi0 );

  // now we can use composite Simpson, iterating on both axes simultaneously
  // This is like using Fubini over and over again for sampled Emu ==> integrate
  // out Epi ==> Simpson again for Emu. Can see more at
  // https://math.stackexchange.com/questions/1319892/simpsons-rule-for-double-integrals.

  const int nSteps = 10000 + 1;
  const double hEMu = ( maxMu - ml ) / ( nSteps - 1 );
  const double hEPi = ( maxPi - mPi ) / ( nSteps - 1 );
  const double preSimp = hEMu * hEPi / ( 9.0 * ( nSteps - 1 ) * ( nSteps - 1 ) );

  double intNow = 0.0;
  for( int i = 0; i < nSteps; i++ ){
    for( int j = 0; j < nSteps; j++ ){
      double midW = 0.0;
      //determine midpoint coefficient for this step
      if( i % (nSteps - 1) == 0 ){ // edge case i
	if( j % (nSteps - 1) == 0 ){ midW = 1.0; } // edge case j
	else if( j % 2 == 0 ){ midW = 2.0; } // even j
	else{ midW = 4.0; } // odd j
      }
      else if( i % 2 == 0 ){ // even i
	if( j % (nSteps - 1) == 0 ){ midW = 2.0; } // edge case j
	else if( j % 2 == 0 ){ midW = 4.0; } // even j
	else{ midW = 8.0; } // odd j
      }
      else{ // odd i
	if( j % (nSteps - 1) == 0 ){ midW = 4.0; } // edge case j
	else if( j % 2 == 0 ){ midW = 8.0; } // even j
	else{ midW = 16.0; } // odd j
      }
      // finally, evaluate f at this point
      const double xev  = mPi + i * hEPi;
      const double yev  = ml + j * hEMu;
      const double fev  = f->Eval( xev, yev );

      // and add to integral
      intNow += std::abs( preSimp * midW * fev );
    }
  }

  delete f;
    
  intNow *= preFac * bigMats;

  return intNow;
	    
}
//----------------------------------------------------------------------------
// *especially* this channel, there's N4 in the propagator so it emits *both* the pi-zeros!!!
// It is subleading in |U_\ell 4|^2, therefore not important to get this exactly right
double BRCalculator::DWidth_Pi0Pi0Nu( const double M,
					      const double Ue42, const double Umu42, const double Ut42 ) const
{ 
  const double preFac = fpi2 * fpi2 * GF2 * GF2 * std::pow( M, 5.0 ) / ( 64.0 * pi*pi*pi );

  const double Ue4 = std::sqrt( Ue42 );
  const double Um4 = std::sqrt( Umu42 );
  const double Ut4 = std::sqrt( Ut42 );

  // once again, assume all PMNS matrix elements real
  const double bigMats = std::pow( Ue4 * ( Ue1 + Ue2 + Ue3 ) +
				   Um4 * ( Um1 + Um2 + Um3 ) +
				   Ut4 * ( Ut1 + Ut2 + Ut3 ), 2.0 );
  const double smallMats = std::pow( Ue42 + Umu42 + Ut42 , 2.0 );

  // let's make the limits
  const double maxNu = 
    ( ( M - mPi0 ) * ( M - mPi0 ) - mPi0*mPi0 ) / ( 2.0 * ( M - mPi0 ) );
  const double maxPi = 
    ( ( M - mPi0 ) * ( M - mPi0 ) + mPi0*mPi0 ) / ( 2.0 * ( M - mPi0 ) );

  // gotta put in the formula
  TF2 * f = new TF2( "fPi0Pi0Nu", Pi0Pi0NuForm, mPi0, maxPi, 0.0, maxNu, 2 );
  f->SetParameter( 0, M );
  f->SetParameter( 1, mPi0 );

  // using composite Simpson to evaluate
  
  const int nSteps = 10000 + 1;
  const double hENu = ( maxNu - 0.0 ) / ( nSteps - 1 );
  const double hEPi = ( maxPi - mPi0 ) / ( nSteps - 1 );
  const double preSimp = hENu * hEPi / ( 9.0 * ( nSteps - 1 ) * ( nSteps - 1 ) );

  double intNow = 0.0;
  for( int i = 0; i < nSteps; i++ ){
    for( int j = 0; j < nSteps; j++ ){
      double midW = 0.0;
      //determine midpoint coefficient for this step
      if( i % (nSteps - 1) == 0 ){ // edge case i
	if( j % (nSteps - 1) == 0 ){ midW = 1.0; } // edge case j
	else if( j % 2 == 0 ){ midW = 2.0; } // even j
	else{ midW = 4.0; } // odd j
      }
      else if( i % 2 == 0 ){ // even i
	if( j % (nSteps - 1) == 0 ){ midW = 2.0; } // edge case j
	else if( j % 2 == 0 ){ midW = 4.0; } // even j
	else{ midW = 8.0; } // odd j
      }
      else{ // odd i
	if( j % (nSteps - 1) == 0 ){ midW = 4.0; } // edge case j
	else if( j % 2 == 0 ){ midW = 8.0; } // even j
	else{ midW = 16.0; } // odd j
      }
      // finally, evaluate f at this point
      const double xev  = mPi0 + i * hEPi;
      const double yev  = 0.0 + j * hENu;
      const double fev  = f->Eval( xev, yev );

      // and add to integral
      intNow += std::abs( preSimp * midW * fev );
    }
  }

  delete f;

  intNow *= preFac * bigMats * smallMats;

  return intNow;
}
//----------------------------------------------------------------------------
// formula for N --> pi pi0 ell decay rate
double BRCalculator::PiPi0EllForm( double *x, double *par ){
    double MN = par[0];
    double MMu = par[1];
    double MPi = par[2];
    double MPi0 = par[3];
    
    double Epi = x[0];
    double Emu = x[1];

    double pi0Term = ( MN - Emu - Epi > MPi0 ) ? 
      std::sqrt( std::pow( ( MN - Emu - Epi ), 2.0 ) - MPi0 * MPi0 ) : 0.0;
    
    double ETerm =
      std::sqrt( Emu*Emu - MMu*MMu ) *
      std::sqrt( Epi*Epi - MPi*MPi ) *
      pi0Term / ( MN - Emu - Epi );
    
    double FracNum1 = MN*MN - 2.0*( MN-Emu-Epi )*MN + MPi0*MPi0;
    double FracNum2 = MN*MN - 2.0*Emu*MN + 2.0*MMu*MMu;
    double FracNum3 = MN*MN - MPi0*MPi0;
    double FracNum4 = MN*MN - 2.0*( MN-Emu-Epi )*MN + MPi0*MPi0 + MMu*MMu - MPi*MPi;
    double FracNum = FracNum1*FracNum2 - FracNum3*FracNum4;
    double FracDen = std::pow( MN*MN - 2.0*( MN - Emu - Epi ) * MN + MPi0*MPi0 , 2.0 );
    
    return ETerm * FracNum / FracDen;
}
//----------------------------------------------------------------------------
// formula for N --> pi0 pi0 nu decay rate
double BRCalculator::Pi0Pi0NuForm( double *x, double *par ){
    double MN = par[0];
    double MPi0 = par[1];
    
    double Epi = x[0]; // leading pi-zero energy
    double Enu = x[1];

    double ETerm = 
      std::sqrt( Epi*Epi - MPi0*MPi0 ) *
      (Enu + MN) * Enu * Enu *
      (MN - Enu - Epi);

    double Frac1 = 1.0 / ( Enu * ( MN - Enu - Epi ) + MPi0 * MPi0 - MN * MN );
    double Frac2 = 1.0 / ( Enu * Epi + MPi0 * MPi0 - MN * MN );

    return ETerm * std::pow( ( Frac1 + Frac2 ), 2.0 );
}
