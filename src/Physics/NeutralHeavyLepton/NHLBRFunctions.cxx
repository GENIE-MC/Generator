#include "NHLBRFunctions.h"

using namespace genie;
using namespace genie::NHL;

// initialise the parameters
void NHLSelector::InitParameters() {
  LOG( "NHL", pDEBUG ) << "Initialising parameters from config files. . .";

  wAng = utils::nhl::GetCfgDouble( "Param", "WeakInt", "WeinbergAngle" );
  s2w = std::pow( std::sin( wAng ), 2.0 );

  Vud = utils::nhl::GetCfgDouble( "Param", "CKM", "CKM-Vud" );
  Vud2 = Vud * Vud;

  fpi = utils::nhl::GetCfgDouble( "NHL", "External", "Pion-FFactor" );
  fpi2 = fpi * fpi;

  BR_C1 = 1./4. * ( 1. - 4. * s2w + 8. * s2w * s2w );
  BR_C2 = 1./2. * ( -s2w + 2. * s2w * s2w );

  Ue1 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Ue1" );
  Ue2 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Ue2" );
  Ue3 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Ue3" );
  Um1 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Um1" );
  Um2 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Um2" );
  Um3 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Um3" );
  Ut1 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Ut1" );
  Ut2 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Ut2" );
  Ut3 = genie::utils::nhl::GetCfgDouble( "NHL", "External", "PMNS-Ut3" );

  fParamsInitialised = true;
}

// Get Coloma et al's form factor functions
double NHLSelector::GetColomaF1( double x ) {
  if( !fParamsInitialised ) InitParameters();

  if( x < 0. || x > 0.5 ) { LOG( "NHL", pERROR ) << "BRFunctions::GetColomaF1:: Illegal x = " << x; exit( 3 ); }
  if( x == 0.5 ) return 0.;
  int i = x/NHLSelector::PARTWIDTH;
  if( x - i*NHLSelector::PARTWIDTH ==0 ) return NHLSelector::ColomaF1[i];
  return 1./2. * ( NHLSelector::ColomaF1[i] + NHLSelector::ColomaF1[i+1] );
}

double NHLSelector::GetColomaF2( double x ) {
  if( !fParamsInitialised ) InitParameters();

  if( x < 0. || x > 0.5 ) { LOG( "NHL", pERROR ) << "BRFunctions::GetColomaF2:: Illegal x = " << x; exit( 3 ); }
  if( x == 0.5 ) return 0.;
  int i = x/NHLSelector::PARTWIDTH;
  if( x - i*NHLSelector::PARTWIDTH==0 ) return NHLSelector::ColomaF2[i];
  return 1./2. * ( NHLSelector::ColomaF2[i] + NHLSelector::ColomaF2[i+1] );
}

// total decay widths, various channels
double NHLSelector::DWidth_PiZeroAndNu( const double M, const double Ue42, const double Umu42, const double Ut42 ) {
  if( !fParamsInitialised ) InitParameters();

  const double x       = genie::utils::nhl::MassX( mPi0, M );
  const double preFac  = GF2 * M*M*M / ( 32. * pi );
  const double kinPart = ( 1. - x*x ) * ( 1. - x*x );
  return preFac * ( Ue42 + Umu42 + Ut42 ) * fpi2 * kinPart;
}

double NHLSelector::DWidth_PiAndLepton( const double M, const double Ua42, const double ma ) {
  if( !fParamsInitialised ) InitParameters();

  const double xPi     = genie::utils::nhl::MassX( mPi, M );
  const double xLep    = genie::utils::nhl::MassX( ma, M );
  const double preFac  = GF2 * M*M*M / ( 16. * pi );
  const double kalPart = TMath::Sqrt( genie::utils::nhl::Kallen( 1, xPi*xPi, xLep*xLep ) );
  const double othPart = 1. - xPi*xPi - xLep*xLep * ( 2. + xPi*xPi - xLep*xLep );
  return preFac * fpi2 * Ua42 * Vud2 * kalPart * othPart;
}

double NHLSelector::DWidth_Invisible( const double M, const double Ue42, const double Umu42, const double Ut42 ) {
  if( !fParamsInitialised ) InitParameters();
  
  const double preFac = GF2 * TMath::Power( M, 5. ) / ( 192. * pi*pi*pi );
  return preFac * ( Ue42 + Umu42 + Ut42 );
}

double NHLSelector::DWidth_SameLepton( const double M, const double Ue42, const double Umu42, const double Ut42, const double mb, bool bIsMu ) {
  if( !fParamsInitialised ) InitParameters();

  const double preFac = GF2 * TMath::Power( M, 5. ) / ( 192. * pi*pi*pi );
  const double x      = genie::utils::nhl::MassX( mb, M );
  const double f1     = GetColomaF1( x );
  const double f2     = GetColomaF2( x );
  const double C1Part = ( Ue42 + Umu42 + Ut42 ) * f1 * BR_C1;
  const double C2Part = ( Ue42 + Umu42 + Ut42 ) * f2 * BR_C2;
  const double D1Part = bIsMu ? 2. * s2w * Umu42 * f1 : 2. * s2w * Ue42 * f1;
  const double D2Part = bIsMu ? s2w * Umu42 * f2 : s2w * Ue42 * f2;
  return preFac * ( C1Part + C2Part + D1Part + D2Part );
}

double NHLSelector::DWidth_DiffLepton( const double M, const double Ua42, const double Ub42, const int IsMajorana ) {
  if( !fParamsInitialised ) InitParameters();

  const double preFac = GF2 * TMath::Power( M, 5. ) / ( 192. * pi*pi*pi );
  const double x = genie::utils::nhl::MassX( mMu, M );
  const double kinPol = 1. - 8. * x*x + 8. * TMath::Power( x, 6. ) - TMath::Power( x, 8. );
  const double kinLn  = -12. * TMath::Power( x, 4. ) * TMath::Log( x*x );
  const double kinPart = kinPol + kinLn;
  const double coupPart = IsMajorana ? Ua42 : Ua42 + Ub42; // 2nd diagram in Majorana case!
  return preFac * kinPart * coupPart;
}

// note that these BR are very very tiny.
double NHLSelector::DWidth_PiPi0Ell( const double M, const double ml,
					      const double Ue42, const double Umu42, const double Ut42,
					      const bool isElectron)
{
  if( !fParamsInitialised ) InitParameters();

  // because the actual decay width is very hard to integrate onto a full DWidth,
  // build 2Differential and then integrate numerically
  // using Simpson's method for 2D.

  const double preFac = fpi2 * fpi2 * GF2 * GF2 * Vud2 * M / ( 32.0 * pi*pi*pi );
  const double Ua1 = isElectron ? Ue1 : Um1;
  const double Ua2 = isElectron ? Ue2 : Um2;
  const double Ua3 = isElectron ? Ue3 : Um3;
  const double Ua4 = isElectron ? std::sqrt( Ue42 ) : std::sqrt( Umu42 );

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

// *especially* this channel, there's N4 in the propagator so it emits *both* the pi-zeros!!!
// It is subleading in |U_\ell 4|^2, therefore not important to get this exactly right
double NHLSelector::DWidth_Pi0Pi0Nu( const double M,
					      const double Ue42, const double Umu42, const double Ut42 )
{ 
  if( !fParamsInitialised ) InitParameters();

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

// differential decay width for NHL channels!

void NHLSelector::Diff1Width_PiAndLepton_CosTheta( const double M, const double Ua42,
							    const double ml,
							    double &thePreFac, 
							    double &theCnstPart,
							    double &thePropPart ) {
  if( !fParamsInitialised ) InitParameters();

  const double preFac   = 1. / ( 32.0 * pi * M*M*M );
  const double sqrKal   = std::sqrt( genie::utils::nhl::Kallen( M*M, mPi*mPi, ml*ml ) );
  const double formPart = fpi2 * Ua42 * Vud2 * GF2;
  const double parConst = std::pow( ( M*M - ml*ml ), 2.0 ) - mPi*mPi*( M*M + ml*ml );
  const double parCoeff = -1.0 * ( M*M - ml*ml ) * std::sqrt( genie::utils::nhl::Kallen( M*M, mPi*mPi, ml*ml ) );
  
  thePreFac   = preFac * sqrKal * formPart;
  theCnstPart = parConst;
  thePropPart = parCoeff; // modulo |P| * cos(theta)
}

// formula for N --> pi pi0 ell decay rate
double NHLSelector::PiPi0EllForm( double *x, double *par ){
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

// formula for N --> pi0 pi0 nu decay rate
double NHLSelector::Pi0Pi0NuForm( double *x, double *par ){
  if( !fParamsInitialised ) InitParameters();
  
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
