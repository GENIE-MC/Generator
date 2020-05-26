//____________________________________________________________________________
/*
  Copyright (c) 2003-2019, The GENIE Collaboration
  For the full text of the license visit http://copyright.genie-mc.org
  or see $GENIE/LICENSE

  Author: Marco Roda
  University of Liverpool
  <mroda \at liverpool.ac.uk>

  For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"

#include "Physics/Coherent/XSection/COHDeltaCurrent.h"

using namespace genie;


COHDeltaCurrent::COHDeltaCurrent() :
  COHFormFactorI("genie::COHDeltaCurrent")
{

}
//____________________________________________________________________________
COHDeltaCurrent::COHDeltaCurrent(string config) :
  COHFormFactorI("genie::COHDeltaCurrent", config)
{

}
//____________________________________________________________________________
COHDeltaCurrent::~COHDeltaCurrent()
{

}
//____________________________________________________________________________
GTrace R( const Interaction * i, 
	  const COHFormFactorI * ff ) const {
  
  GTrace t ;
  return t ;
}
//____________________________________________________________________________
void COHDeltaCurrent::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDeltaCurrent::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void COHDeltaCurrent::LoadConfig(void)
{

  bool good_config = true ;

  const Algorithm * algo = SubAlg( "DeltaFormFactor" ) ; 

  if ( ! delta_ff ) {

    good_config = false ;
    LOG("COHDeltaCurrent", pERROR ) << "Configuration does not know DeltaFormFactor field" ;
  }

  
  
  delta_ff = dynamic_cast<const DeltaTransitionFormFactor*>( algo ) ;

  if ( ! delta_ff ) {
    good_config = false ; 
    LOG("COHDeltaCurrent", pERROR ) << "Algo " << algo -> Id().Name() << " is not a DeltaTransitionFormFactor" ;
  }


  if ( ! good_config ) {
    LOG("COHDeltaCurrent", pFATAL ) << "Bad configuration, exit" ;
    exit( 78 ) ;
  }

}
//____________________________________________________________________________
GTrace DirTrace( const Interaction * i,
		 const COHFormFactorI * ff ) const {

  // these calculations expects the interaction to be in the lab frame with the incoming neutrino parallel to z

  TLorentzVector * probe = i -> InitState().GetProbeP4( kRfLab ) ;
  TLorentzVector   out_neutrino = i -> Kine().FSLeptonP4() ; 
  TLorentzVector   t_photon = i -> Kine().HadSystP4() ; 

  TLorentzVector t_q = *probe - out_neutrino ; 
  double Q2 = t_q.M2() ;

  delete probe ; 

    // this requires Q2
  double C3V   = delta_ff -> C3V( Q2 ) ;
  double C3vNC = delta_ff -> C3VNC( Q2 ) ;
  double C5aNC = delta_ff -> C5ANC( Q2 ) ;

  double mDelta = util::res::Mass( this -> Resonance() ) ; 
  double mDelta2 = pow( mDelta, 2 ); 
  
  // the following contractions requires a vector with time coordinate in the 0-th position
  // that is not the case for TLorentzVector so we need to rearrange it
  std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() } ;
  double p0 = out_neutrino.E() ; 
  std::array<double, 4> kg = { t_photon.E(), t_photon.X(), t_photon.Y(), t_photon.Z() } ;
  double mn = kNucleonMass ;
  double mn2 = pow( mn, 2 ) ;


  GTrace tr ;

  tr[0][0] = -(C3v*C3vNC*( mDelta*p0 + mn*(p0 + q[0]))*( 4*mDelta2*q[1]*kg[1] - 4*pow(p0,2)*q[1]*kg[1] 
							 + pow(q[1],3)*kg[1] + q[1]*q[3]*q[3]*kg[1] 
							 + q[1]*q[1]*kg[1]*kg[1] - q[3]*q[3]*kg[1]*kg[1]
							 - q[3]*q[3]*kg[2]*kg[2] 
							 + q[0]*q[0]*( q[1]*q[1] - 3*q[1]*kg[1] + q[3]*( 2*q[3] - 3*kg[3] ) ) 
							 + q[3]*( 4*mDelta2 - 4*pow(p0,2) + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] )*kg[3] 
							 - 8*p0*q[0]*( q[1]*kg[1] + q[3]*kg[3] ) 
							 ) 
	       ) / ( 3.*mDelta2*pow(mn,3) );
  
  tr[0][1] = ( std::complex( 0., 1. )*C3v*C5aNC*mDelta*( mDelta*p0 + mn*( p0 + q[0] ) )*q[3]*kg[2]*( 8*mn2 
												     - 2*q[0]*( p0 + q[0] ) 
												     + q[1]*( q[1] + kg[1] ) + q[3]*( q[3] + kg[3]) ) 
	       - C3v*C3vNC*mn*( 16*mDelta2 * mn * q[0] * ( p0 + q[0] ) * kg[1] 
				+ 4*pow(mDelta,3)*( 4*p0*q[0]*kg[1] + q[3]*( -( q[3]*kg[1]) + q[1]*kg[3]) ) 
				- 4*mn*(p0 + q[0])*( 4*pow(p0,2)*q[0]*kg[1] + 8*p0*q[0] * q[0]*kg[1] 
						     + 4*pow( q[0], 3)*kg[1] - q[0]*q[1] * q[1]*kg[1] 
						     + 2*p0*q[3]*(-(q[3]*kg[1]) + q[1]*kg[3]) 
						     - q[0]*q[1]*( 2*kg[1]*kg[1] + kg[2]*kg[2] + kg[3]*(- q[3] + kg[3]) ) 
						     - q[0]*kg[1]*( q[0]*q[0] + q[3]*(2*q[3] + kg[3] ) ) ) 
				+ mDelta*( -16*pow(p0,3)*q[0]*kg[1] - q[3]*(q[3]*kg[1] - q[1]*kg[3])*( -3*q[0]*q[0] + q[1]*q[1] 
												       + 2*q[1]*kg[1] + q[3]*( q[3] + 2*kg[3] ) ) 
					   - 4*pow(p0,2)*(8*q[0]*q[0]*kg[1] + 3*q[3]*( -( q[3]*kg[1] ) + q[1]*kg[3] ) ) 
					   + 4*p0*q[0]*( q[0]*q[0]*( q[1] - 3*kg[1] ) + q[1]*q[1]*kg[1] + q[3]*kg[1]*( 4*q[3] + kg[3] ) 
							 + q[1]*( kg[1]*kg[1] - 3*q[3]*kg[3] ) ) )
				) 
	       ) / (12.*mDelta2*pow(mn,4));


  tr[0][2] = ( C3v*( std::complex<double>(0,-1)*C5aNC*mDelta*( mDelta*p0 + mn*( p0 + q[0] ) )*( q[3]*kg[1] - q[1]*kg[3] )*( 8*mn2 - 
															    2*q[0]*( p0 + q[0] ) + 
															    q[1]*( q[1] + kg[1] ) + 
															    q[3]*( q[3] + kg[3] ) ) +
		     C3vNC*mn*kg[2]*( -16*mDelta2*mn*q[0]*( p0 + q[0] ) + 
				      4*pow(mDelta,3)*( -4*p0*q[0] + q[1]*q[1] + q[3]*q[3]) +
				      4*mn*( p0 + q[0] )*( 4*pow(p0,2)*q[0] + 8*p0*q[0]*q[0] + 3*pow(q[0],3) - 
							   2*p0*(q[1]*q[1] + q[3]*q[3] ) -
							   q[0]*q[1]*( 2*q[1] + kg[1] ) - 
							   q[0]*q[3]*( 2*q[3] + kg[3] ) ) +
				      mDelta*( 16*pow(p0,3)*q[0] + 4*pow(p0,2)*( 8*q[0]*q[0] - 3*( q[1]*q[1] + q[3]*q[3] ) ) +
					       (q[1]*q[1] + q[3]*q[3])*( -3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) +
					       4*p0*q[0]*( 3*q[0]*q[0] - q[1]*( 4*q[1] + kg[1] ) - 
							   q[3]*( 4*q[3] + kg[3]) ) ) )
		     ) ) / ( 12.*mDelta2*pow(mn,4) );
  
  
  tr[0][3] = ( C3v*( std::complex<double>(0,-1)*C5aNC*mDelta*( mDelta*p0 + mn*( p0 + q[0] ) )*q[1]*kg[2]*( 8*mn2 - 2*q[0]*( p0 + q[0] ) 
													   + q[1]*( q[1] + kg[1] ) 
													   + q[3]*( q[3] + kg[3] ) ) 
		     + C3vNC*mn*( -16*mDelta2*mn*q[0]*( p0 + q[0] )*kg[3] + 4*pow(mDelta,3)*( -(q[1]*q[3]*kg[1]) 
											      - 4*p0*q[0]*kg[3] + q[1]*q[1]*kg[3] ) 
				  + 4*mn*( p0 + q[0] )*( q[0]*q[3]*( kg[1]*( q[1] + kg[1] ) + kg[2]*kg[2] ) 
							 + 8*p0*q[0]*q[0]*kg[3] - q[0]*( -4*pow(p0,2) + q[3]*q[3] + q[1]*( 2*q[1] + kg[1] ) )*kg[3] 
							 + pow(q[0],3)*( -2*q[3] + 3*kg[3] ) + 2*p0*q[1]*( q[3]*kg[1] - q[1]*kg[3] ) ) 
				  + mDelta*( 16*pow(p0,3)*q[0]*kg[3] + 4*pow(p0,2)*( 3*q[1]*q[3]*kg[1] + 8*q[0]*q[0]*kg[3] - 3*q[1]*q[1]*kg[3] ) 
					     + q[1]*( -( q[3]*kg[1] ) + q[1]*kg[3] )*( -3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) 
					     - 4*p0*q[0]*( -( q[3]*(3*q[1]*kg[1] + kg[1]*kg[1] + kg[2]*kg[2] ) ) + q[0]*q[0]*( 2*q[3] - 3*kg[3] ) 
							   + ( q[3]*q[3] + q[1]*( 4*q[1] + kg[1] ) )*kg[3] ) ) ) ) 
	       ) /( 12.*mDelta2*pow(mn,4) );

  
}
//____________________________________________________________________________
GTrace DirTrace( const Interaction * i,
		 const COHFormFactorI * ff ) const {


}
