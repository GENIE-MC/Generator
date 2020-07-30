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
using namespace genie::utils::math;


COHDeltaCurrent::COHDeltaCurrent() :
  COHHadronicCurrentI("genie::COHDeltaCurrent")
{

}

//____________________________________________________________________________
COHDeltaCurrent::COHDeltaCurrent(string config) :
  COHHadronicCurrentI("genie::COHDeltaCurrent", config)
{

}

//____________________________________________________________________________
COHDeltaCurrent::~COHDeltaCurrent()
{

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
GTrace COHDeltaCurrent::R( const Interaction * i,
	  const COHFormFactorI * ff ) const {

  const InitialState & init_state = i -> InitState();

  double k0 = init_state.ProbeE( kRfLab ) ;
  TLorentzVector probe( 0., 0., k0, k0 ) ; 
  TLorentzVector out_neutrino = i -> Kine().FSLeptonP4() ;
  TLorentzVector t_photon = i -> Kine().HadSystP4() ;
  // TODO verify the calculation of t
  TLorentzVector t = 0.5*( (probe - out_neutrino) - t_photon );

  // This is not quite what is used in original code which was
  // the magnitude of the 3-momentum, here we are using 4-momentum
  TLorentzVector p( t.Vect(), 
		    sqrt(constants::kNucleonMass + t.Vect().Mag2() ) );

  double p2 = p.Mag2() ;
  double p0 = out_neutrino.E() ;

  int pdg = init_state.Tgt().Pdg() ;

  // Right now the proton and neutron FF are equal but could change
  double ff_p = ff -> ProtonFF( t.Mag(), pdg );
  double ff_n = ff -> NeutronFF( t.Mag(), pdg );

  std::complex<double> D_prop_dir   = Delta_med -> AverageDirectPropagator( p2, pdg ) ;
  std::complex<double> D_prop_cross = Delta_med -> AverageCrossPropagator( p2 ) ;

  GTrace R = DirTrace(i) ;
  R *= D_prop_dir ;

  GTrace tr_cross = CrsTrace(i) ;
  tr_cross *= D_prop_cross ;

  // Add trace * propagator from direct and crossed diagrams
  R += tr_cross ;
  R *= ( ff_p + ff_n ) / (2*p0) ;

  return R;
}

//____________________________________________________________________________
GTrace COHDeltaCurrent::DirTrace( const Interaction * i ) const {

  // these calculations expects the interaction to be in the lab frame with the incoming neutrino parallel to z
  double k0 = i -> InitState().ProbeE( kRfLab ) ;
  TLorentzVector probe( 0., 0., k0, k0 ) ; 
  TLorentzVector out_neutrino = i -> Kine().FSLeptonP4() ; 
  TLorentzVector t_photon = i -> Kine().HadSystP4() ; 

  TLorentzVector t_q = probe - out_neutrino ; 
  double Q2 = t_q.Mag2() ; // Downstream agnostic to sign

    // this requires Q2
  double C3v   = delta_ff -> C3V( Q2 ) ;
  double C3vNC = delta_ff -> C3VNC( Q2 ) ;
  double C5aNC = delta_ff -> C5ANC( Q2 ) ;


  double mDelta = utils::res::Mass( Resonance() ) ;
  double mDelta2 = pow( mDelta, 2 ); 
  
  // the following contractions requires a vector with time coordinate in the 0-th position
  // that is not the case for TLorentzVector so we need to rearrange it
  std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() } ;
  double p0 = out_neutrino.E() ; 
  std::array<double, 4> kg = { t_photon.E(), t_photon.X(), t_photon.Y(), t_photon.Z() } ;
  double mn = constants::kNucleonMass ;
  double mn2 = pow( mn, 2 ) ;


  GTrace tr;

  tr[0][0] = -( C3v*C3vNC*( mDelta*p0 + mn*(p0 + q[0]) )*( 4*mDelta2*q[1]*kg[1] - 4*pow(p0,2)*q[1]*kg[1]
							    + pow(q[1],3)*kg[1] + q[1]*q[3]*q[3]*kg[1]
							    + q[1]*q[1]*kg[1]*kg[1] - q[3]*q[3]*kg[1]*kg[1]
							    - q[3]*q[3]*kg[2]*kg[2]
							    + q[0]*q[0]*( q[1]*q[1] - 3*q[1]*kg[1] + q[3]*( 2*q[3] - 3*kg[3] ) )
							    + q[3]*( 4*mDelta2 - 4*pow(p0,2) + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] )*kg[3]
							    - 8*p0*q[0]*( q[1]*kg[1] + q[3]*kg[3] ) )
	            ) / ( 3.*mDelta2*pow(mn,3) );

  tr[0][1] = ( std::complex<double>(0.,1.)*C3v*C5aNC*mDelta*( mDelta*p0 + mn*( p0 + q[0] ) )*q[3]*kg[2]*( 8*mn2
                  - 2*q[0]*( p0 + q[0] ) + q[1]*( q[1] + kg[1] ) + q[3]*( q[3] + kg[3]) )
              - C3v*C3vNC*mn*( 16*mDelta2 * mn * q[0] * ( p0 + q[0] ) * kg[1] + 4*pow(mDelta,3)*( 4*p0*q[0]*kg[1]
                  + q[3]*( -( q[3]*kg[1]) + q[1]*kg[3]) ) - 4*mn*(p0 + q[0])*( 4*pow(p0,2)*q[0]*kg[1] + 8*p0*q[0] * q[0]*kg[1]
						      + 4*pow( q[0], 3)*kg[1] - q[0]*q[1] * q[1]*kg[1] + 2*p0*q[3]*(-(q[3]*kg[1]) + q[1]*kg[3])
						      - q[0]*q[1]*( 2*kg[1]*kg[1] + kg[2]*kg[2] + kg[3]*(- q[3] + kg[3]) )
						      - q[0]*kg[1]*( q[0]*q[0] + q[3]*(2*q[3] + kg[3] ) ) )
				        + mDelta*( -16*pow(p0,3)*q[0]*kg[1] - q[3]*(q[3]*kg[1] - q[1]*kg[3])*( -3*q[0]*q[0] + q[1]*q[1]
				          + 2*q[1]*kg[1] + q[3]*( q[3] + 2*kg[3] ) ) - 4*pow(p0,2)*(8*q[0]*q[0]*kg[1] + 3*q[3]*( -( q[3]*kg[1] )
				          + q[1]*kg[3] ) ) + 4*p0*q[0]*( q[0]*q[0]*( q[1] - 3*kg[1] ) + q[1]*q[1]*kg[1] + q[3]*kg[1]*( 4*q[3] + kg[3] )
							  + q[1]*( kg[1]*kg[1] - 3*q[3]*kg[3] ) ) ) )
	            ) / (12.*mDelta2*pow(mn,4));

  tr[0][2] = ( C3v*( std::complex<double>(0.,-1.)*C5aNC*mDelta*( mDelta*p0 + mn*( p0 + q[0] ) )*( q[3]*kg[1]
                  - q[1]*kg[3] )*( 8*mn2 - 2*q[0]*( p0 + q[0] ) + q[1]*( q[1] + kg[1] ) + q[3]*( q[3] + kg[3] ) )
              + C3vNC*mn*kg[2]*( -16*mDelta2*mn*q[0]*( p0 + q[0] ) + 4*pow(mDelta,3)*( -4*p0*q[0] + q[1]*q[1] + q[3]*q[3])
                  + 4*mn*( p0 + q[0] )*( 4*pow(p0,2)*q[0] + 8*p0*q[0]*q[0] + 3*pow(q[0],3) - 2*p0*(q[1]*q[1] + q[3]*q[3] )
                  - q[0]*q[1]*( 2*q[1] + kg[1] ) - q[0]*q[3]*( 2*q[3] + kg[3] ) )
                + mDelta*( 16*pow(p0,3)*q[0] + 4*pow(p0,2)*( 8*q[0]*q[0] - 3*( q[1]*q[1] + q[3]*q[3] ) )
                  + (q[1]*q[1] + q[3]*q[3])*( -3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                  + 4*p0*q[0]*( 3*q[0]*q[0] - q[1]*( 4*q[1] + kg[1] ) - q[3]*( 4*q[3] + kg[3]) ) ) ) )
              ) / ( 12.*mDelta2*pow(mn,4) );

  tr[0][3] = ( C3v*( std::complex<double>(0.,-1.)*C5aNC*mDelta*( mDelta*p0 + mn*( p0 + q[0] ) )*q[1]*kg[2]*( 8*mn2 - 2*q[0]*( p0 + q[0] )
                  + q[1]*( q[1] + kg[1] ) + q[3]*( q[3] + kg[3] ) )
		          + C3vNC*mn*( -16*mDelta2*mn*q[0]*( p0 + q[0] )*kg[3] + 4*pow(mDelta,3)*( -(q[1]*q[3]*kg[1])
                  - 4*p0*q[0]*kg[3] + q[1]*q[1]*kg[3] ) + 4*mn*( p0 + q[0] )*( q[0]*q[3]*( kg[1]*( q[1] + kg[1] ) + kg[2]*kg[2] )
							    + 8*p0*q[0]*q[0]*kg[3] - q[0]*( -4*pow(p0,2) + q[3]*q[3] + q[1]*( 2*q[1] + kg[1] ) )*kg[3]
							    + pow(q[0],3)*( -2*q[3] + 3*kg[3] ) + 2*p0*q[1]*( q[3]*kg[1] - q[1]*kg[3] ) )
				        + mDelta*( 16*pow(p0,3)*q[0]*kg[3] + 4*pow(p0,2)*( 3*q[1]*q[3]*kg[1] + 8*q[0]*q[0]*kg[3] - 3*q[1]*q[1]*kg[3] )
					        + q[1]*( -( q[3]*kg[1] ) + q[1]*kg[3] )*( -3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
					        - 4*p0*q[0]*( -( q[3]*(3*q[1]*kg[1] + kg[1]*kg[1] + kg[2]*kg[2] ) ) + q[0]*q[0]*( 2*q[3] - 3*kg[3] )
							    + ( q[3]*q[3] + q[1]*( 4*q[1] + kg[1] ) )*kg[3] ) ) ) )
	            ) / ( 12.*mDelta2*pow(mn,4) );

  tr[1][0] = (C3v*( std::complex<double>(0.,1.)*C5aNC*q[3]*kg[2]*( 8*pow(mDelta,3)*mn*( p0 + q[0] )
								  - 4*mn2*( p0 + q[0] )*( 2*p0*q[0] + q[0]*q[0] - q[1]*kg[1] - q[3]*kg[3] ) 
								  - 2*mDelta*mn*( p0 + q[0] )*( 8*mn2 + 2*p0*q[0] + q[0]*q[0] - q[1]*kg[1] - q[3]*kg[3] ) 
                + mDelta2*( 8*pow(p0,3) + 12*pow(p0,2)*q[0] + 4*p0*q[0]*q[0] + 8*mn2*( -2*p0 + q[0] )
                  + q[0]*( q[1]*q[1] - q[1]*kg[1] + q[3]*( q[3] - kg[3] ) )
                  - 2*p0*( q[1]*q[1] + 2*q[1]*kg[1] + q[3]*( q[3] + 2*kg[3] ) ) ) )
		          - 2*C3vNC*mn*( 16*mDelta2*mn*q[0]*( p0 + q[0] )*q[1] + 4*pow(mDelta,3)*( q[1]*( 4*p0*q[0]
		              - q[0]*q[0] + kg[1]*kg[1] ) + q[3]*kg[1]*kg[3] )
					      + mDelta*( -16*pow(p0,3)*q[0]*q[1] + 4*p0*q[0]*( pow(q[1],3) + 2*q[1]*q[1]*kg[1] - 3*q[1]*kg[1]*kg[1]
                  + q[3]*kg[1]*( q[3] - 3*kg[3] ) + q[1]*q[3]*( q[3] + kg[3] ) )
						      + ( 3*q[0]*q[0] - q[1]*( q[1] + 2*kg[1] ) - q[3]*( q[3] + 2*kg[3] ) )*( q[0]*q[0]*q[1]
                  - kg[1]*( q[1]*kg[1] + q[3]*kg[3] ) ) - 4*pow(p0,2)*( 5*q[0]*q[0]*q[1] + 3*kg[1]*( q[1]*kg[1] + q[3]*kg[3] ) ) )
					        - 4*mn*( p0 + q[0] )*( 4*pow(p0,2)*q[0]*q[1] + 6*p0*q[0]*q[0]*q[1] + 2*p0*kg[1]*( q[1]*kg[1] + q[3]*kg[3] )
								  + q[0]*( 2*q[0]*q[0]*q[1] - pow(q[1],3) - 2*q[1]*q[1]*kg[1] + q[3]*kg[1]*(-q[3] + kg[3] )
                  + q[1]*( kg[1]*kg[1] - q[3]*(q[3] + kg[3] ) ) ) ) ) )
	            ) / ( 24.*mDelta2*pow(mn,4) );

  tr[1][1] = (C3v*( std::complex<double>(0.,1.)*C5aNC*q[3]*kg[2]*( 4*pow(mDelta,3)*mn*( q[1] + kg[1] )
                  + 2*mn2*( q[1] + kg[1] )*( -( q[0]*( 2*p0 + q[0] ) ) + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*( q[1] + kg[1] )*( -4*p0*q[0] - 3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3]
                  + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                + mDelta2*( q[0]*q[0]*( q[1] - kg[1] ) + 8*mn2*kg[1] - 2*p0*q[0]*kg[1] + q[3]*q[3]*kg[1]
                  + 4*pow(p0,2)*( q[1] + kg[1] ) - q[1]*q[3]*kg[3] ) )
              - 2*C3vNC*mn*( 16*mDelta2*mn*( p0 + q[0] )*( q[0]*q[0] - q[3]*kg[3] )
					        + 4*pow(mDelta,3)*( -( q[0]*( -4*p0*q[0] + q[0]*q[0] + q[3]*q[3] - kg[1]*kg[1] ) )
                  + 2*( -2*p0 + q[0] )*q[3]*kg[3] )
					        - 4*mn*( p0 + q[0] )*( 4*pow(p0,2)*q[0]*q[0] + 2*pow(q[0],4)
								  - q[0]*q[0]*( q[1]*q[1] + 2*q[1]*kg[1] - kg[1]*kg[1] + 3*q[3]*kg[3] )
								  + q[3]*( -( q[3]*( 2*kg[1]*kg[1] + kg[2]*kg[2] ) )
                  +( q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] )*kg[3] )
								  + 2*p0*q[0]*( 3*q[0]*q[0] + kg[1]*kg[1] - q[3]*( q[3] + 2*kg[3] ) ) )
					      - mDelta*( 16*pow(p0,3)*q[0]*q[0] + 4*pow(p0,2)*q[0]*( 5*q[0]*q[0] - 3*q[3]*q[3] + 3*kg[1]*kg[1] - 2*q[3]*kg[3] )
						      + 4*p0*( -2*q[3]*q[3]*kg[1]*kg[1] - q[0]*q[0]*( 2*q[3]*q[3] +( q[1] - kg[1] )*( q[1] + 3*kg[1] ) )
                  + q[3]*( q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] )*kg[3] )
						      - q[0]*( q[0]*q[0] + q[3]*q[3] - kg[1]*kg[1] - 2*q[3]*kg[3] )*( 3*q[0]*q[0] - q[1]*( q[1] + 2*kg[1] )
                  - q[3]*( q[3] + 2*kg[3] ) ) ) ) )
	            ) / ( 24.*mDelta2*pow(mn,4) );

  tr[1][2] = ( C3v*( -2*C3vNC*mn*kg[2]*( 16*mDelta2*mn*( p0 + q[0] )*q[1] + 4*pow(mDelta,3)*( 4*p0*q[1] + q[0]*( -2*q[1] + kg[1] ) ) 
					        + 4*mn*( p0 + q[0] )*( pow(q[1],3) + q[1]*q[3]*q[3] + 2*q[1]*q[1]*kg[1] + q[3]*q[3]*kg[1]
								  - 2*p0*q[0]*( 2*q[1] + kg[1] ) -q[0]*q[0]*( 2*q[1] + kg[1] ) + q[1]*q[3]*kg[3] )
                + mDelta*(-4*pow(p0,2)*q[0]*( 2*q[1] + 3*kg[1] ) + 4*p0*( q[0]*q[0]*( 2*q[1] - 3*kg[1] )
                  + ( q[1]*q[1] + q[3]*q[3] )*( q[1] + 2*kg[1] ) )
						      - q[0]*( 2*q[1] - kg[1] )*( -3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
		          + std::complex<double>(0.,1.)*C5aNC*( 2*mn2*q[3]*kg[2]*kg[2]*( -( q[0]*( 2*p0 + q[0] ) ) + q[1]*kg[1] + q[3]*kg[3] )
							    + 4*pow(mDelta,3)*mn*( -( q[3]*kg[1]*kg[1] ) + q[0]*q[0]*( q[3] - 2*kg[3] ) - 2*p0*q[0]*kg[3]
                  + ( 8*mn2 + q[3]*q[3] + q[1]*( q[1] + kg[1] ) )*kg[3] )
							  + mDelta2*( -2*p0*q[3]*( 2*p0*( -2*q[0]*q[0] + kg[1]*kg[1] )
							    + q[0]*( -3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + q[1]*kg[1] + kg[1]*kg[1] ) )
                  + q[3]*( -( q[0]*( 4*p0 + q[0] ) ) + q[1]*q[1] + q[3]*q[3])*kg[2]*kg[2]
                  + 2*p0*( -4*pow(p0,2)*q[0] - 2*pow(q[0],3) + q[0]*q[1]*( q[1] + kg[1] )
                  + 2*p0*( -3*q[0]*q[0] + q[3]*q[3] + q[1]*( q[1] + kg[1] ) ) )*kg[3]
                  + 8*mn2*( q[3]*( -2*p0*q[0] + kg[2]*kg[2] ) + 2*p0*( 2*p0 + q[0] )*kg[3] ) )
							  + mDelta*mn*( q[3]*( -pow(q[0],4) - kg[1]*kg[1]*( 3*( q[1]*q[1] + q[3]*q[3] ) + 4*q[1]*kg[1] ) 
							    - 2*(q[3]*q[3] + q[1]*(q[1] + kg[1] ) )*kg[2]*kg[2]
                  + q[0]*q[0]*( q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 5*kg[1]*kg[1] + 2*kg[2]*kg[2] ) )
                  + ( 2*pow(q[0],4) + pow((q[1]*q[1] + q[3]*q[3]),2) + 3*q[1]*( q[1]*q[1] + q[3] * q[3])*kg[1]
                  + 2*( q[1] - q[3] )*( q[1] + q[3] )*kg[1]*kg[1]
                  - q[0]*q[0]*( 3*( q[1]*q[1] + q[3]*q[3] ) + 5*q[1]*kg[1] ) )*kg[3] + 4*pow(p0,2)*q[0]*q[0]*( q[3] + kg[3] )
									- 8*mn2*( 2*p0*q[0]*( q[3] + kg[3] ) + q[0]*q[0]*( 2*q[3] + kg[3] ) 
									- kg[3]*( q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
									+ 2*p0*q[0]*( -( q[1]*q[1]*( q[3] + 2*kg[3] ) ) + q[0]*q[0]*( q[3] + 3*kg[3] ) 
									- q[1]*kg[1]*( q[3] + 3*kg[3] ) + q[3]*( 3*kg[1]*kg[1] + kg[2]*kg[2] - q[3]*( q[3] + 3*kg[3] ) ) ) ) ) )
	            ) / ( 24.*mDelta2*pow(mn,4) );

  tr[1][3] = -( C3v*( std::complex<double>(0.,1.)*C5aNC*kg[2]*( 4*pow(mDelta,3)*mn*( 8*mn2 - 2*q[0]*(p0 + q[0]) + q[1]*(q[1] + kg[1]) )
                  - 2*mn2*q[3]*(q[3] + kg[3])*( -(q[0]*(2*p0 + q[0])) + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*( 4*pow(p0,2)*q[0]*q[0] + 6*p0*pow(q[0],3) + 2*pow(q[0],4)
                  - 4*p0*q[0]*q[1]*q[1] - 3*q[0]*q[0]*q[1]*q[1] + pow(q[1],4)
                  + q[1]*q[1]*q[3]*q[3] - 6*p0*q[0]*q[1]*kg[1] - 5*q[0]*q[0]*q[1]*kg[1]
                  + 3*pow(q[1],3)*kg[1] + q[1]*q[3] * q[3]*kg[1] + 2*q[1]*q[1]*kg[1]*kg[1]
                  + 2*q[3]*( -(q[0]*(p0 + q[0])) + q[1]*(q[1] + kg[1]))*kg[3]
                  + 8*mn2*( -( q[0]*(2*p0 + q[0]) ) + q[1]*q[1] + 2*q[1]*kg[1] + q[3]*( q[3] + 2*kg[3] ) ) )
                + mDelta2*( -8*pow(p0,3)*q[0] + 4*pow(p0,2)*( -3*q[0]*q[0] + q[1]*(q[1] + kg[1]) )
                  + 8*mn2*( 2*p0*(2*p0 + q[0]) - q[3]*kg[3] )
                  + 2*p0*q[0]*( -2*q[0] * q[0] + q[1]*(q[1] + kg[1]) + q[3]*(q[3] + 2*kg[3]) )
                  + q[3]*( q[0]*q[0]*(-q[3] + kg[3]) + q[1]*(q[3]*kg[1] - q[1]*kg[3]) ) ) )
              + 2*C3vNC*mn*( 16*mDelta2*mn*(p0 + q[0])*q[1]*kg[3]
                  + 4*pow(mDelta,3)*( q[0]*q[1]*(q[3] - 2*kg[3]) + 4*p0*q[1]*kg[3] + q[0]*kg[1]*kg[3] )
                  - 4*mn*(p0 + q[0])*( q[1]*q[3]*( 2*kg[1]*kg[1] + kg[2]*kg[2] )
                  - q[1]*( q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] )*kg[3]
                  + q[0]*q[0]*( -(q[3]*(q[1] + kg[1])) + (2*q[1] + kg[1])*kg[3] )
                  + 2*p0*q[0]*( kg[1]*kg[3] + q[1]*(q[3] + 2*kg[3]) ) )
                + mDelta*( -4*pow(p0,2)*q[0]*( 3*q[1]*q[3] + 2*q[1]*kg[3] + 3*kg[1]*kg[3] )
                  - q[0]*( q[1]*(q[3] - 2*kg[3]) + kg[1]*kg[3] )*( 3*q[0]*q[0] - q[1]*(q[1]
                  + 2*kg[1]) - q[3]*(q[3] + 2*kg[3]) )
                  + 4*p0*( -2*q[1]*q[3]*kg[1] * kg[1]
                  + q[1]*( q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] )*kg[3]
                  + q[0]*q[0]*( -(q[1]*q[3]) + q[3]*kg[1] + 2*q[1]*kg[3] - 3*kg[1]*kg[3] ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][0] = ( C3v*( -2*C3vNC*mn*kg[2]*( 4*pow(mDelta,3)*(q[1]*kg[1] + q[3]*kg[3])
                  - 4*mn*(p0 + q[0])*( -( q[0]*(q[1]*q[1] + q[3]*q[3] - q[1]*kg[1] - q[3]*kg[3]) )
                  + 2*p0*(q[1]*kg[1] + q[3]*kg[3]) )
                + mDelta*( 4*p0*q[0]*( q[1]*q[1] + q[3]*q[3] - 3*q[1]*kg[1] - 3*q[3]*kg[3] )
                  - 12*pow(p0,2)*(q[1]*kg[1] + q[3]*kg[3])
                  + (q[1]*kg[1] + q[3]*kg[3])*( -3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
              - std::complex<double>(0.,1.)*C5aNC*(q[3]*kg[1] - q[1]*kg[3])*( 8*pow(mDelta,3)*mn*(p0 + q[0])
                  - 4*mn2*(p0 + q[0])*( 2*p0*q[0] + q[0]*q[0] - q[1]*kg[1] - q[3]*kg[3] )
                  - 2*mDelta*mn*(p0 + q[0])*( 8*mn2 + 2*p0*q[0] + q[0]*q[0] - q[1]*kg[1] - q[3]*kg[3] )
                + mDelta2*( 8*pow(p0,3) + 12*pow(p0,2)*q[0] + 4*p0*q[0]*q[0] + 8*mn2*(-2*p0 + q[0])
                  + q[0]*( q[1]*q[1] - q[1]*kg[1] + q[3]*(q[3] - kg[3]) )
                  - 2*p0*( q[1]*q[1] + 2*q[1]*kg[1] + q[3]*(q[3] + 2*kg[3]) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][1] = ( C3v*( -2*C3vNC*mn*kg[2]*( 4*pow(mDelta,3)*q[0]*kg[1] - 4*mn*(p0 + q[0])*( 2*p0*q[0]*kg[1]
                  + q[0] * q[0]*(-q[1] + kg[1]) + q[3]*(-(q[3]*kg[1]) + q[1]*kg[3]) )
                + mDelta*( 4*p0*q[0] * q[0]*(q[1] - 3*kg[1]) - 12*pow(p0,2)*q[0]*kg[1]
                  + 8*p0*q[3]*(q[3]*kg[1] - q[1]*kg[3])
                  + q[0]*kg[1]*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
              - std::complex<double>(0.,1.)*C5aNC*( 2*mn2*(q[1] + kg[1])*(-(q[3]*kg[1])
                  + q[1]*kg[3])*( 2*p0*q[0] + q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] )
                  + 4*pow(mDelta,3)*mn*( q[3]*( q[0] * q[0] + q[1]*kg[1] - kg[2] * kg[2] )
                  + ( 8*mn2 - 2*q[0]*(p0 + q[0]) + q[3] * q[3] )*kg[3] )
                + mDelta2*( -8*pow(p0,3)*q[0]*kg[3] + 2*p0*q[0]*( -( q[3]*(q[1] * q[1] + q[3] * q[3] + q[1]*kg[1]
                  + 2*kg[1] * kg[1] + kg[2] * kg[2] ) ) + q[0] * q[0]*(3*q[3] - 2*kg[3]) + q[1]*(q[1] + 2*kg[1])*kg[3] )
                  + 4*pow(p0,2)*( q[3]*kg[1]*(q[1] + kg[1]) + q[0] * q[0]*(q[3] - 3*kg[3]) + q[3]*kg[3]*(q[3] + kg[3]) )
                  + (q[3]*kg[1] - q[1]*kg[3])*( q[0] * q[0]*(q[1] - kg[1]) + q[3]*(q[3]*kg[1] - q[1]*kg[3]) )
                  + 8*mn2*( 4*pow(p0,2)*kg[3] + 2*p0*q[0]*(-q[3] + kg[3]) + kg[1]*(q[3]*kg[1] - q[1]*kg[3]) ) )
                + mDelta*mn*( q[3]*( -pow(q[0],4) + kg[1]*( pow(q[1],3) + q[1]*q[3] * q[3]
                  + 2*q[1] * q[1]*kg[1] - 2*q[3] * q[3]*kg[1] )
                  - ( q[1] * q[1] + 3*q[3] * q[3] + 2*q[1]*kg[1] )*kg[2] * kg[2]
                  + q[0] * q[0]*( -q[1] * q[1] + q[3] * q[3] - 3*q[1]*kg[1] + 2*kg[1] * kg[1] + 5*kg[2] * kg[2] ) )
                  + ( 2*pow(q[0],4) - q[0] * q[0]*( 3*q[3] * q[3] + 2*q[1]*kg[1] )
                  + q[3] * q[3]*( q[1] * q[1] + q[3] * q[3] + 4*q[1]*kg[1] - 2*kg[2] * kg[2] ) )*kg[3]
                  + 4*pow(p0,2)*q[0] * q[0]*(q[3] + kg[3])
                  - 8*mn2*( 2*p0*q[0]*(q[3] + kg[3]) + q[0] * q[0]*(2*q[3] + kg[3])
                  - kg[3]*( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
                  + 2*p0*q[0]*( -(q[1] * q[1]*q[3]) - q[1]*kg[1]*(3*q[3] + kg[3]) + q[0] * q[0]*(q[3] + 3*kg[3])
                  + q[3]*( kg[1] * kg[1] + 3*kg[2] * kg[2] - q[3]*(q[3] + 3*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][2] = ( C3v*( std::complex<double>(0.,-1.)*C5aNC*kg[2]*(q[3]*kg[1] - q[1]*kg[3])*( 4*pow(mDelta,3)*mn
                + mDelta2*( 8*mn2 + 4*pow(p0,2) - 2*p0*q[0] - q[0]*q[0] + q[1]*q[1] + q[3]*q[3] )
                  + 2*mn2*( -(q[0]*(2*p0 + q[0])) + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*( -4*p0*q[0] - 3*q[0]*q[0] + q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
              + 2*C3vNC*mn*( -16*mDelta2*mn*(p0 + q[0])*(q[0]*q[0] - q[1]*kg[1] - q[3]*kg[3])
                  + 4*pow(mDelta,3)*( q[0]*( q[0]*q[0] + q[1]*q[1] + q[3]*q[3] - 2*q[1]*kg[1] - kg[2]*kg[2] - 2*q[3]*kg[3] )
                  + 4*p0*( -q[0]*q[0] + q[1]*kg[1] + q[3]*kg[3] ) )
                  + 4*mn*(p0 + q[0])*( 4*pow(p0,2)*q[0]*q[0] + 2*pow(q[0],4) + pow(q[1],3)*kg[1] + q[1]*q[3]*q[3]*kg[1]
                  + q[1]*q[1]*kg[1]*kg[1] - q[3]*q[3]*kg[1]*kg[1] - q[1]*q[1]*kg[2]*kg[2] - 2*q[3]*q[3]*kg[2]*kg[2]
                  + q[3]*( q[1]*q[1] + q[3]*q[3] + 2*q[1]*kg[1] )*kg[3]
                  - q[0]*q[0]*( q[1]*q[1] + 3*q[1]*kg[1] - kg[2]*kg[2] + 3*q[3]*kg[3] )
                  +  2*p0*q[0]*( 3*q[0]*q[0] - q[1]*(q[1] + 2*kg[1]) + kg[2]*kg[2] - q[3]*(q[3] + 2*kg[3]) ) )
                + mDelta*( 16*pow(p0,3)*q[0]*q[0] + 4*pow(p0,2)*q[0]*( 5*q[0]*q[0] - 3*q[1]*q[1] - 3*q[3]*q[3]
                  - 2*q[1]*kg[1] + 3*kg[2]*kg[2] - 2*q[3]*kg[3] )
                  - q[0]*( q[0]*q[0] + q[1]*q[1] + q[3]*q[3]
                  - 2*q[1]*kg[1] - kg[2]*kg[2] - 2*q[3]*kg[3] )*( 3*q[0]*q[0] - q[1]*(q[1]
                  + 2*kg[1]) - q[3]*(q[3] + 2*kg[3]) )
                  + 4*p0*( pow(q[1],3)*kg[1] + q[1]*q[3]*q[3]*kg[1] + q[3]*q[3]*( -2*kg[2]*kg[2] + q[3]*kg[3] )
                  + q[0]*q[0]*( q[1]*(-4*q[1] + kg[1]) + 3*kg[2]*kg[2] + q[3]*(-2*q[3] + kg[3]) )
                  + q[1]*q[1]*( 2*kg[1]*kg[1] + kg[3]*(q[3] + 2*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][3] = ( C3v*( -2*C3vNC*mn*kg[2]*( 4*pow(mDelta,3)*q[0]*kg[3] - 4*mn*(p0 + q[0])*( 2*p0*q[0]*kg[3]
                  + q[0] * q[0]*(-q[3] + kg[3]) + q[1]*(q[3]*kg[1] - q[1]*kg[3]) )
                + mDelta*( 4*p0*q[0] * q[0]*(q[3] - 3*kg[3]) - 12*pow(p0,2)*q[0]*kg[3] + 8*p0*q[1]*( -(q[3]*kg[1]) + q[1]*kg[3] )
                  + q[0]*kg[3]*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
              + std::complex<double>(0.,1.)*C5aNC*( -2*mn2*(q[3] + kg[3])*( q[3]*kg[1]
                  - q[1]*kg[3] )*( -( q[0]*(2*p0 + q[0]) ) + q[1]*kg[1] + q[3]*kg[3] )
                  + 4*pow(mDelta,3)*mn*( q[0] * q[0]*(q[1] - 2*kg[1]) + 8*mn2*kg[1]
                  - 2*p0*q[0]*kg[1] + q[1]*( q[1]*kg[1] - kg[2] * kg[2] + q[3]*kg[3] ) )
                + mDelta2*( -8*pow(p0,3)*q[0]*kg[1] + 4*pow(p0,2)*( q[0] * q[0]*(q[1] - 3*kg[1])
                  + q[1]*kg[1]*(q[1] + kg[1]) + q[1]*kg[3]*(q[3] + kg[3]) )
                  + (q[3]*kg[1] - q[1]*kg[3])*( q[0] * q[0]*(-q[3] + kg[3]) + q[1]*(q[3]*kg[1] - q[1]*kg[3]) )
                  + 8*mn2*( 4*pow(p0,2)*kg[1] + 2*p0*q[0]*(-q[1] + kg[1]) + kg[3]*( -(q[3]*kg[1]) + q[1]*kg[3] ) )
                  + 2*p0*q[0]*( -pow(q[1],3) + q[0] * q[0]*(q[1] - 2*kg[1]) + q[3]*kg[1]*(q[3] + 2*kg[3])
                  + q[1]*( 2*kg[1] * kg[1] + kg[2] * kg[2] - q[3]*(q[3] + kg[3]) ) ) )
                + mDelta*mn*( pow(q[0],4)*q[1] - q[0] * q[0]*pow(q[1],3) + q[0] * q[0]*q[1]*q[3] * q[3] + 2*pow(q[0],4)*kg[1]
                  - 3*q[0] * q[0]*q[1] * q[1]*kg[1] + pow(q[1],4)*kg[1] + q[1] * q[1]*q[3] * q[3]*kg[1]
                  - 2*q[0] * q[0]*q[1]*kg[1] * kg[1] + 2*pow(q[1],3)*kg[1] * kg[1] - 2*q[1]*q[3] * q[3]*kg[1] * kg[1]
                  + 4*pow(p0,2)*q[0] * q[0]*(q[1] + kg[1]) + 3*q[0] * q[0]*q[1]*kg[2] * kg[2] - pow(q[1],3)*kg[2] * kg[2]
                  - 3*q[1]*q[3] * q[3]*kg[2] * kg[2] - 2*q[1] * q[1]*kg[1]*kg[2] * kg[2]
                  - 3*q[0] * q[0]*q[1]*q[3]*kg[3] + pow(q[1],3)*q[3]*kg[3] + q[1]*pow(q[3],3)*kg[3]
                  - 2*q[0] * q[0]*q[3]*kg[1]*kg[3] + 4*q[1] * q[1]*q[3]*kg[1]*kg[3] - 2*q[1]*q[3]*kg[2] * kg[2]*kg[3]
                  - 8*mn2*( 2*p0*q[0]*(q[1] + kg[1]) + q[0] * q[0]*(2*q[1] + kg[1])
                  - kg[1]*( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
                  + 2*p0*q[0]*( -pow(q[1],3) - 3*q[1] * q[1]*kg[1] + 2*q[0] * q[0]*(q[1] + 2*kg[1])
                  - q[1]*( q[3] * q[3] + kg[1] * kg[1] - 2*kg[2] * kg[2] + 3*q[3]*kg[3] )
                  - kg[1]*( kg[1] * kg[1] + kg[2] * kg[2] + kg[3]*(q[3] + kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][0] = -( C3v*( std::complex<double>(0.,1.)*C5aNC*q[1]*kg[2]*( 8*pow(mDelta,3)*mn*(p0 + q[0])
                  - 4*mn2*(p0 + q[0])*( 2*p0*q[0] + q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] )
                - 2*mDelta*mn*(p0 + q[0])*( 8*mn2 + 2*p0*q[0] + q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] )
                + mDelta2*( 8*pow(p0,3) + 12*pow(p0,2)*q[0] + 4*p0*q[0] * q[0] + 8*mn2*(-2*p0 + q[0])
                  + q[0]*( q[1] * q[1] - q[1]*kg[1] + q[3]*(q[3] - kg[3]) )
                  - 2*p0*( q[1] * q[1] + 2*q[1]*kg[1] + q[3]*(q[3] + 2*kg[3]) ) ) )
              + 2*C3vNC*mn*( 16*mDelta2*mn*q[0]*(p0 + q[0])*q[3]
                  + 4*pow(mDelta,3)*( 4*p0*q[0]*q[3] - q[3]*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[1]*kg[1]*kg[3] )
                  + 4*mn*(p0 + q[0])*( -4*pow(p0,2)*q[0]*q[3] - 8*p0*q[0] * q[0]*q[3] - 3*pow(q[0],3)*q[3]
                  + 2*p0*q[3]*( kg[1] * kg[1] + kg[2] * kg[2] )
                  + q[0]*q[3]*( q[1] * q[1] + q[3] * q[3] + q[1]*kg[1] + kg[1] * kg[1] + kg[2] * kg[2] )
                  - 2*p0*q[1]*kg[1]*kg[3] + q[0]*( q[1] * q[1] + 2*q[3] * q[3] - q[1]*kg[1] )*kg[3] )
                + mDelta*( -16*pow(p0,3)*q[0]*q[3] - 4*pow(p0,2)*( 8*q[0] * q[0]*q[3]
                  - 3*q[3]*( kg[1] * kg[1] + kg[2] * kg[2] ) + 3*q[1]*kg[1]*kg[3] )
                  - ( q[3]*(kg[1] * kg[1] + kg[2] * kg[2]) - q[1]*kg[1]*kg[3] )*( -3*q[0] * q[0] + q[1] * q[1]
                  + 2*q[1]*kg[1] + q[3]*(q[3] + 2*kg[3]) )
                  + 4*p0*q[0]*( -3*q[0] * q[0]*q[3] + q[1]*kg[1]*(q[3] - 3*kg[3]) + q[1] * q[1]*(q[3] + kg[3])
                  + q[3]*( 3*(kg[1] * kg[1] + kg[2] * kg[2]) + q[3]*(q[3] + 2*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][1] = ( std::complex<double>(0.,1.)*C3v*C5aNC*kg[2]*( -2*mn2*q[1]*(q[1] + kg[1])*( -(q[0]*(2*p0 + q[0]))
                  + q[1]*kg[1] + q[3]*kg[3] ) + 4*pow(mDelta,3)*mn*( 8*mn2 - 2*q[0]*(p0 + q[0]) + q[3]*(q[3] + kg[3]) )
                - mDelta*mn*( -4*pow(p0,2)*q[0] * q[0] - 2*pow(q[0],4) + q[0] * q[0]*q[3] * q[3]
                  - q[1] * q[1]*q[3] * q[3] - pow(q[3],4)
                  + 2*q[0] * q[0]*q[1]*kg[1] - 2*q[1]*q[3] * q[3]*kg[1]
                  + 2*q[3] * q[3]*kg[1] * kg[1] + 2*q[3] * q[3]*kg[2] * kg[2]
                  - q[3]*( -5*q[0] * q[0] + q[1] * q[1] + 3*q[3] * q[3] + 2*q[1]*kg[1] )*kg[3]
                  + 2*p0*q[0]*( -3*q[0] * q[0] + 2*q[3] * q[3] + q[1]*kg[1] + 3*q[3]*kg[3] )
                  + 8*mn2*( 2*p0*q[0] + q[0] * q[0] - q[1]*(q[1] + 2*kg[1]) - q[3]*(q[3] + 2*kg[3]) ) )
                + mDelta2*( -8*pow(p0,3)*q[0] + 8*mn2*( 2*p0*(2*p0 + q[0]) - q[1]*kg[1] )
                  + 4*pow(p0,2)*( -3*q[0] * q[0] + q[3]*(q[3] + kg[3]) )
                  + 2*p0*q[0]*( -2*q[0] * q[0] + q[1]*(q[1] + 2*kg[1]) + q[3]*(q[3] + kg[3]) )
                  + q[1]*( q[0] * q[0]*(-q[1] + kg[1]) + q[3]*( -(q[3]*kg[1]) + q[1]*kg[3] ) ) ) )
              - 2*C3v*C3vNC*mn*( 16*mDelta2*mn*(p0 + q[0])*q[3]*kg[1]
                  + 4*pow(mDelta,3)*( q[0]*q[1]*q[3] + 4*p0*q[3]*kg[1] + q[0]*kg[1]*(-2*q[3] + kg[3]) )
                  - 4*mn*(p0 + q[0])*( 2*p0*q[0]*( q[1]*q[3] + kg[1]*(2*q[3] + kg[3]) )
                  + q[0] * q[0]*( q[1]*(q[3] - kg[3]) + kg[1]*(2*q[3] + kg[3]) )
                  - q[3]*( q[1] * q[1]*kg[1] + q[1]*( 2*kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) )
                + mDelta*( -4*pow(p0,2)*q[0]*( 3*q[1]*q[3] + 2*q[3]*kg[1] + 3*kg[1]*kg[3] )
                  + q[0]*( q[1]*q[3] + kg[1]*(-2*q[3] + kg[3]) )*( -3*q[0] * q[0] + q[1] * q[1] + 2*q[1]*kg[1] + q[3]*(q[3] + 2*kg[3]) )
                  + 4*p0*( q[0] * q[0]*( -3*q[1]*q[3] + 2*q[3]*kg[1] + q[1]*kg[3] - 3*kg[1]*kg[3] )
                  + q[3]*( q[1] * q[1]*kg[1] + 2*q[1]*(kg[1] * kg[1] + kg[2] * kg[2]) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][2] = ( C3v*( -2*C3vNC*mn*kg[2]*( 16*mDelta2*mn*(p0 + q[0])*q[3]
                  - 4*mn*(p0 + q[0])*( -( q[3]*( -2*q[0]*(2*p0 + q[0]) + q[3] * q[3] + q[1]*(q[1] + kg[1]) ) )
                  + ( 2*p0*q[0] + q[0] * q[0] - q[1] * q[1] - 2*q[3] * q[3])*kg[3] )
                  + 4*pow(mDelta,3)*( 4*p0*q[3] + q[0]*(-2*q[3] + kg[3]) )
                + mDelta*( -4*pow(p0,2)*q[0]*(2*q[3] + 3*kg[3])
                  - q[0]*(2*q[3] - kg[3])*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                  + 4*p0*( q[0] * q[0]*(2*q[3] - 3*kg[3]) + ( q[1] * q[1] + q[3] * q[3] )*(q[3] + 2*kg[3]) ) ) )
              - std::complex<double>(0.,1.)*C5aNC*( 2*mn2*q[1]*kg[2] * kg[2]*( -(q[0]*(2*p0 + q[0]))
                  + q[1]*kg[1] + q[3]*kg[3]) + 4*pow(mDelta,3)*mn*( 8*mn2*kg[1] - 2*p0*q[0]*kg[1] - 2*q[0] * q[0]*kg[1]
                  + q[1] * q[1]*kg[1] + q[3] * q[3]*kg[1] + q[1]*kg[1] * kg[1] + q[1]*kg[2] * kg[2] + q[3]*kg[1]*kg[3])
                + mDelta2*( -8*pow(p0,3)*q[0]*kg[1] + q[1]*( -q[0] * q[0] + q[1] * q[1] + q[3] * q[3] )*kg[2] * kg[2]
                  + 8*mn2*( 4*pow(p0,2)*kg[1] + 2*p0*q[0]*(-q[1] + kg[1]) + q[1]*kg[2] * kg[2] )
                  + 4*pow(p0,2)*( q[0] * q[0]*(q[1] - 3*kg[1]) + q[1] * q[1]*kg[1] + q[1]*( kg[1] * kg[1] + kg[2] * kg[2] )
                  + q[3]*kg[1]*(q[3] + kg[3])) + 2*p0*q[0]*( -pow(q[1],3) + 2*q[0] * q[0]*(q[1] - kg[1])
                  + q[3]*kg[1]*(q[3] + kg[3]) - q[1]*( -kg[1] * kg[1] + kg[2] * kg[2] + q[3]*(q[3] + kg[3]) ) ) )
                + mDelta*mn*( 4*pow(q[0],4)*q[1] - 2*q[0] * q[0]*pow(q[1],3) - 2*q[0] * q[0]*q[1]*q[3] * q[3] + 2*pow(q[0],4)*kg[1]
                  - 5*q[0] * q[0]*q[1] * q[1]*kg[1] + pow(q[1],4)*kg[1] - q[0] * q[0]*q[3] * q[3]*kg[1]
                  + 2*q[1] * q[1]*q[3] * q[3]*kg[1] + pow(q[3],4)*kg[1] - 5*q[0] * q[0]*q[1]*kg[1] * kg[1]
                  + 3*pow(q[1],3)*kg[1] * kg[1] + 3*q[1]*q[3] * q[3]*kg[1] * kg[1]
                  + 2*q[1] * q[1]*pow(kg[1],3) - 2*q[3] * q[3]*pow(kg[1],3) + 4*pow(p0,2)*q[0] * q[0]*(q[1] + kg[1])
                  - 3*q[0] * q[0]*q[1]*kg[2] * kg[2] + pow(q[1],3)*kg[2] * kg[2] + q[1]*q[3] * q[3]*kg[2] * kg[2]
                  + 2*q[1] * q[1]*kg[1]*kg[2] * kg[2] - 2*q[3] * q[3]*kg[1]*kg[2] * kg[2] - 2*q[0] * q[0]*q[1]*q[3]*kg[3]
                  - 5*q[0] * q[0]*q[3]*kg[1]*kg[3] + 3*q[1] * q[1]*q[3]*kg[1]*kg[3] + 3*pow(q[3],3)*kg[1]*kg[3]
                  + 4*q[1]*q[3]*kg[1] * kg[1]*kg[3] + 2*q[1]*q[3]*kg[2] * kg[2]*kg[3]
                  - 8*mn2*( 2*p0*q[0]*(q[1] + kg[1]) + q[0] * q[0]*(2*q[1] + kg[1])
                  - kg[1]*( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
                  + 2*p0*q[0]*( -pow(q[1],3) - 3*q[1] * q[1]*kg[1] + q[0] * q[0]*(4*q[1] + 3*kg[1]) - q[3]*kg[1]*(2*q[3] + 3*kg[3])
                  - q[1]*( 3*kg[1] * kg[1] + 2*kg[2] * kg[2] + q[3]*(q[3] + kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][3] = -( C3v*( std::complex<double>(0.,1.)*C5aNC*q[1]*kg[2]*( 4*pow(mDelta,3)*mn*(q[3] + kg[3])
                  + 2*mn2*(q[3] + kg[3])*( -( q[0]*(2*p0 + q[0]) ) + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*(q[3] + kg[3])*( -4*p0*q[0] - 3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                + mDelta2*( -(q[1]*q[3]*kg[1]) + q[0] * q[0]*(q[3] - kg[3]) + 8*mn2*kg[3] - 2*p0*q[0]*kg[3] + q[1] * q[1]*kg[3]
                  + 4*pow(p0,2)*(q[3] + kg[3]) ) )
              + 2*C3vNC*mn*( 16*mDelta2*mn*(p0 + q[0])*( q[0] * q[0] - q[1]*kg[1] )
                  + 4*pow(mDelta,3)*( 4*p0*(q[0] * q[0] - q[1]*kg[1]) - q[0]*( pow(q[1] - kg[1],2) + kg[2] * kg[2] ) )
                  - 4*mn*(p0 + q[0])*( 4*pow(p0,2)*q[0] * q[0] + 3*pow(q[0],4)
                  + 2*p0*q[0]*( 4*q[0] * q[0] - pow(q[1] + kg[1],2) - kg[2] * kg[2] )
                  - q[0] * q[0]*( q[3] * q[3] + (q[1] + kg[1])*(2*q[1] + kg[1]) + kg[2] * kg[2] + 2*q[3]*kg[3] )
                  + q[1]*( q[1] * q[1]*kg[1] + q[1]*( 2*kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) )
                - mDelta*( 16*pow(p0,3)*q[0] * q[0] + 4*pow(p0,2)*q[0]*( 8*q[0] * q[0] - 3*q[1] * q[1]
                  - 2*q[1]*kg[1] - 3*( kg[1] * kg[1] + kg[2] * kg[2] ) )
                  + q[0]*( pow(q[1] - kg[1],2) + kg[2] * kg[2] )*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3]
                  + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) + 4*p0*( 3*pow(q[0],4) - q[0] * q[0]*( 4*q[1] * q[1]
                  - q[1]*kg[1] + 3*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*(q[3] + 2*kg[3]) )
                  + q[1]*( q[1] * q[1]*kg[1] + 2*q[1]*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  return tr;

}

//____________________________________________________________________________
GTrace COHDeltaCurrent::CrsTrace( const Interaction * i ) const {


  // these calculations expects the interaction to be in the lab frame with the incoming neutrino parallel to z
  double k0 = i->InitState().ProbeE( kRfLab ) ;
  TLorentzVector probe( 0., 0., k0, k0 ) ; 
  TLorentzVector out_neutrino = i -> Kine().FSLeptonP4() ;
  TLorentzVector t_photon = i -> Kine().HadSystP4() ;

  TLorentzVector t_q = probe - out_neutrino ;
  double Q2 = t_q.Mag2() ; // Downstream agnostic to sign

  // this requires Q2
  double C3v   = delta_ff -> C3V( Q2 ) ;
  double C3vNC = delta_ff -> C3VNC( Q2 ) ;
  double C5aNC = delta_ff -> C5ANC( Q2 ) ;

  double mDelta = utils::res::Mass( Resonance() ) ;
  double mDelta2 = pow( mDelta, 2 );

  // the following contractions requires a vector with time coordinate in the 0-th position
  // that is not the case for TLorentzVector so we need to rearrange it
  std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() } ;
  double p0 = out_neutrino.E() ;
  std::array<double, 4> kg = { t_photon.E(), t_photon.X(), t_photon.Y(), t_photon.Z() } ;

  double mn = constants::kNucleonMass ;
  double mn2 = pow( mn, 2 ) ;

  GTrace tr;

  tr[0][0] = -( C3v*C3vNC*( (mDelta + mn)*p0 - mn*q[0] )*( 4*mDelta2*q[1]*kg[1] - 4*pow(p0,2)*q[1]*kg[1] + pow(q[1],3)*kg[1]
                  + q[1]*q[3] * q[3]*kg[1] + q[1] * q[1]*kg[1] * kg[1] - q[3] * q[3]*kg[1] * kg[1] - q[3] * q[3]*kg[2] * kg[2]
                  + q[0] * q[0]*( q[1] * q[1] - 3*q[1]*kg[1] + q[3]*(2*q[3] - 3*kg[3]) ) + q[3]*( 4*mDelta2 - 4*pow(p0,2)
                  + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] )*kg[3] + 8*p0*q[0]*(q[1]*kg[1] + q[3]*kg[3]) )
              ) / (3.*mDelta2*pow(mn,3));

  tr[0][1] = -( C3v*( std::complex<double>(0.,1.)*C5aNC*mDelta*( (mDelta + mn)*p0 - mn*q[0] )*q[3]*kg[2]*( 8*mn2
                  + 2*(p0 - q[0])*q[0] + q[1]*(q[1] + kg[1]) + q[3]*(q[3] + kg[3]) )
              + C3vNC*mn*( 16*mDelta2*mn*(p0 - q[0])*q[0]*kg[1] + 4*pow(mDelta,3)*( 4*p0*q[0]*kg[1] + q[3]*(q[3]*kg[1]
                  - q[1]*kg[3]) ) + 4*mn*(p0 - q[0])*( pow(q[0],3)*(q[1] - 3*kg[1]) + 8*p0*q[0] * q[0]*kg[1]
                  + 2*p0*q[3]*(-(q[3]*kg[1]) + q[1]*kg[3]) + q[0]*( -4*pow(p0,2)*kg[1] + 2*q[3] * q[3]*kg[1]
                  + q[1]*kg[1]*(q[1] + kg[1]) + q[3]*(-q[1] + kg[1])*kg[3] ) )
                + mDelta*( -16*pow(p0,3)*q[0]*kg[1] + q[3]*(q[3]*kg[1] - q[1]*kg[3])*( -3*q[0] * q[0] + q[1] * q[1]
                  + 2*q[1]*kg[1] + q[3]*(q[3] + 2*kg[3]) ) + 4*pow(p0,2)*( 8*q[0] * q[0]*kg[1]
                  + 3*q[3]*( -(q[3]*kg[1]) + q[1]*kg[3] ) ) +
                    4*p0*q[0]*( q[0] * q[0]*(q[1] - 3*kg[1]) + q[1] * q[1]*kg[1] + q[3]*kg[1]*(4*q[3] + kg[3]) + q[1]*(kg[1] * kg[1] - 3*q[3]*kg[3]) ) ) ) )
              ) / (12.*mDelta2*pow(mn,4));

  tr[0][2] = ( std::complex<double>(0.,1.)*C3v*C5aNC*mDelta*( (mDelta + mn)*p0 - mn*q[0] )*( q[3]*kg[1]
                  - q[1]*kg[3] )*( 8*mn2 + 2*(p0 - q[0])*q[0] + q[1]*(q[1] + kg[1]) + q[3]*(q[3] + kg[3]) )
              - C3v*C3vNC*mn*kg[2]*( 16*mDelta2*mn*(p0 - q[0])*q[0] + 4*pow(mDelta,3)*( 4*p0*q[0] + q[1] * q[1] + q[3] * q[3] )
                  + 4*mn*(p0 - q[0])*( -4*pow(p0,2)*q[0] + 8*p0*q[0] * q[0] - 3*pow(q[0],3) - 2*p0*( q[1] * q[1] + q[3] * q[3] )
                  + q[0]*q[1]*(2*q[1] + kg[1]) + q[0]*q[3]*(2*q[3] + kg[3]) )
                + mDelta*( -16*pow(p0,3)*q[0] + 4*pow(p0,2)*( 8*q[0] * q[0] - 3*( q[1] * q[1] + q[3] * q[3] ) )
                  + ( q[1] * q[1] + q[3] * q[3] )*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                  + 4*p0*q[0]*( -3*q[0] * q[0] + q[1]*(4*q[1] + kg[1]) + q[3]*(4*q[3] + kg[3]) ) ) )
              ) / (12.*mDelta2*pow(mn,4));

  tr[0][3] = ( std::complex<double>(0.,1.)*C3v*C5aNC*mDelta*( (mDelta + mn)*p0 - mn*q[0] )*q[1]*kg[2]*( 8*mn2
                  + 2*(p0 - q[0])*q[0] + q[1]*(q[1] + kg[1]) + q[3]*(q[3] + kg[3]) )
              - C3v*C3vNC*mn*( 16*mDelta2*mn*(p0 - q[0])*q[0]*kg[3] + 4*pow(mDelta,3)*( -(q[1]*q[3]*kg[1])
                  + 4*p0*q[0]*kg[3] + q[1] * q[1]*kg[3] ) +4*mn*(p0 - q[0])*( -( q[0]*q[3]*(kg[1]*(q[1] + kg[1]) + kg[2] * kg[2] ) )
                  + pow(q[0],3)*(2*q[3] - 3*kg[3]) + 8*p0*q[0] * q[0]*kg[3] + q[0]*( -4*pow(p0,2) + q[3] * q[3]
                  + q[1]*(2*q[1] + kg[1]) )*kg[3] + 2*p0*q[1]*( q[3]*kg[1] - q[1]*kg[3]) )
                + mDelta*( -16*pow(p0,3)*q[0]*kg[3] + 4*pow(p0,2)*( 3*q[1]*q[3]*kg[1] + 8*q[0] * q[0]*kg[3] - 3*q[1] * q[1]*kg[3] )
                  + q[1]*( -(q[3]*kg[1]) + q[1]*kg[3] )*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                  + 4*p0*q[0]*( -( q[3]*(3*q[1]*kg[1] + kg[1] * kg[1] + kg[2] * kg[2] ) ) + q[0] * q[0]*(2*q[3] - 3*kg[3])
                  + ( q[3] * q[3] + q[1]*(4*q[1] + kg[1]) )*kg[3] ) ) )
              ) / (12.*mDelta2*pow(mn,4));

  tr[1][0] = ( std::complex<double>(0.,1.)*C3v*C5aNC*q[3]*kg[2]*( 8*pow(mDelta,3)*mn*(-p0 + q[0])
                  + 2*mDelta*mn*(p0 - q[0])*( 8*mn2 - 2*p0*q[0] + q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] )
                  - 4*mn2*(p0 - q[0])*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta2*( -8*pow(p0,3) + 12*pow(p0,2)*q[0] + 8*mn2*(2*p0 + q[0])
                  + q[0]*( q[1] * q[1] + q[3] * q[3] - q[1]*kg[1] - q[3]*kg[3] )
                  + 2*p0*( -2*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
              - 2*C3v*C3vNC*mn*( 16*mDelta2*mn*(p0 - q[0])*q[0]*q[1]
                  + 4*pow(mDelta,3)*( 4*p0*q[0]*q[1] + q[0] * q[0]*q[1] - kg[1]*(q[1]*kg[1] + q[3]*kg[3]) )
                + mDelta*( -16*pow(p0,3)*q[0]*q[1] + 4*p0*q[0]*( pow(q[1],3) + 2*q[1] * q[1]*kg[1]
                  - 3*q[1]*kg[1] * kg[1] + q[3]*kg[1]*(q[3] - 3*kg[3]) + q[1]*q[3]*(q[3] + kg[3]) )
                  - ( 3*q[0] * q[0] - q[1]*(q[1] + 2*kg[1]) - q[3]*(q[3] + 2*kg[3]) )*( q[0] * q[0]*q[1]
                  - kg[1]*(q[1]*kg[1] + q[3]*kg[3]) ) + 4*pow(p0,2)*( 5*q[0] * q[0]*q[1] + 3*kg[1]*(q[1]*kg[1] + q[3]*kg[3]) ) )
                - 4*mn*(p0 - q[0])*( 4*pow(p0,2)*q[0]*q[1] - 2*p0*( 3*q[0] * q[0]*q[1] + kg[1]*(q[1]*kg[1] + q[3]*kg[3]) )
                  + q[0]*( 2*q[0] * q[0]*q[1] - pow(q[1],3) - 2*q[1] * q[1]*kg[1] + q[3]*kg[1]*(-q[3] + kg[3])
                  + q[1]*( kg[1] * kg[1] - q[3]*(q[3] + kg[3]) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[1][1] = ( std::complex<double>(0.,1.)*C3v*C5aNC*q[3]*kg[2]*( 4*pow(mDelta,3)*mn*(q[1] + kg[1])
                  + 2*mn2*(q[1] + kg[1])*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*(q[1] + kg[1])*( 4*p0*q[0] - 3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                + mDelta2*( q[0] * q[0]*(q[1] - kg[1]) + 8*mn2*kg[1] + 2*p0*q[0]*kg[1] + q[3] * q[3]*kg[1]
                  + 4*pow(p0,2)*(q[1] + kg[1]) - q[1]*q[3]*kg[3] ) )
              - 2*C3v*C3vNC*mn*( 16*mDelta2*mn*(p0 - q[0])*(q[0] * q[0] - q[3]*kg[3]) + 4*pow(mDelta,3)*( q[0]*( 4*p0*q[0]
                  + q[0] * q[0] + q[3] * q[3] - kg[1] * kg[1] ) - 2*(2*p0 + q[0])*q[3]*kg[3] )
                  - 4*mn*(p0 - q[0])*( 4*pow(p0,2)*q[0] * q[0] + 2*pow(q[0],4) + 2*p0*q[0]*( -3*q[0] * q[0]
                  + q[3] * q[3] - kg[1] * kg[1] + 2*q[3]*kg[3] ) - q[0] * q[0]*( q[1] * q[1] + 2*q[1]*kg[1]
                  - kg[1] * kg[1] + 3*q[3]*kg[3] ) + q[3]*( -( q[3]*(2*kg[1] * kg[1] + kg[2] * kg[2] ) )
                  + ( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1])*kg[3] ) )
                + mDelta*( -16*pow(p0,3)*q[0] * q[0] + 4*pow(p0,2)*q[0]*( 5*q[0] * q[0] - 3*q[3] * q[3]
                  + 3*kg[1] * kg[1] - 2*q[3]*kg[3] ) - q[0]*( q[0] * q[0] + q[3] * q[3] - kg[1] * kg[1]
                  - 2*q[3]*kg[3] )*( 3*q[0] * q[0] - q[1]*(q[1] + 2*kg[1]) - q[3]*(q[3] + 2*kg[3]) )
                  + 4*p0*( q[0] * q[0]*( (q[1] - kg[1])*(q[1] + 3*kg[1]) + q[3]*(2*q[3] - kg[3]) )
                  + q[3]*( 2*q[3]*kg[1] * kg[1] - ( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1])*kg[3] ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[1][2] = ( C3v*( -2*C3vNC*mn*kg[2]*( 16*mDelta2*mn*(p0 - q[0])*q[1] + 4*pow(mDelta,3)*( 4*p0*q[1] + 2*q[0]*q[1]
                  - q[0]*kg[1] ) + 4*mn*(p0 - q[0])*( pow(q[1],3) + q[1]*q[3] * q[3] + 2*q[1] * q[1]*kg[1]
                  + q[3] * q[3]*kg[1] + 2*p0*q[0]*(2*q[1] + kg[1]) - q[0] * q[0]*(2*q[1] + kg[1]) + q[1]*q[3]*kg[3] )
                + mDelta*( 4*pow(p0,2)*q[0]*(2*q[1] + 3*kg[1]) + 4*p0*( q[0] * q[0]*(2*q[1] - 3*kg[1])
                  + ( q[1] * q[1] + q[3] * q[3] )*(q[1] + 2*kg[1]) ) + q[0]*(2*q[1] - kg[1])*( -3*q[0] * q[0]
                  + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
              + std::complex<double>(0.,1.)*C5aNC*( 2*mn2*q[3]*kg[2] * kg[2]*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1]
                  + q[3]*kg[3] ) + 4*pow(mDelta,3)*mn*( -(q[3]*kg[1] * kg[1]) + q[0] * q[0]*(q[3] - 2*kg[3]) + 2*p0*q[0]*kg[3]
                  + ( 8*mn2 + q[3] * q[3] + q[1]*(q[1] + kg[1]) )*kg[3] )
                + mDelta2*( 2*p0*q[3]*( q[0]*( 4*p0*q[0] - 3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] ) + q[0]*q[1]*kg[1]
                  + (-2*p0 + q[0])*kg[1] * kg[1] ) + q[3]*( 4*p0*q[0] - q[0] * q[0] + q[1] * q[1]
                  + q[3] * q[3] )*kg[2] * kg[2] + 2*p0*( 4*pow(p0,2)*q[0] + 2*pow(q[0],3) - q[0]*q[1]*(q[1] + kg[1])
                  + 2*p0*( -3*q[0] * q[0] + q[3] * q[3] + q[1]*(q[1] + kg[1]) ) )*kg[3] + 8*mn2*( q[3]*kg[2] * kg[2]
                  + 2*p0*q[0]*(q[3] - kg[3]) + 4*pow(p0,2)*kg[3] ) )
                + mDelta*mn*( q[3]*( 4*pow(p0,2)*q[0] * q[0] - pow(q[0],4) - kg[1] * kg[1]*( 3*( q[1] * q[1]
                  + q[3] * q[3] ) + 4*q[1]*kg[1] ) - 2*( q[3] * q[3] + q[1]*(q[1] + kg[1]) )*kg[2] * kg[2]
                  + 2*p0*q[0]*( -q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + q[1]*kg[1] - 3*kg[1] * kg[1] - kg[2] * kg[2] )
                  + q[0] * q[0]*( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 5*kg[1] * kg[1] + 2*kg[2] * kg[2] ) )
                  + ( 4*pow(p0,2)*q[0] * q[0] + 2*pow(q[0],4) + pow(q[1] * q[1] + q[3] * q[3],2)
                  + 3*q[1]*( q[1] * q[1] + q[3] * q[3] )*kg[1] + 2*(q[1] - q[3])*(q[1] + q[3])*kg[1] * kg[1]
                  + 2*p0*q[0]*( -3*q[0] * q[0] + 2*q[1] * q[1] + 3*q[3] * q[3] + 3*q[1]*kg[1] )
                  - q[0] * q[0]*( 3*( q[1] * q[1] + q[3] * q[3] ) + 5*q[1]*kg[1] ) )*kg[3]
                  + 8*mn2*( 2*p0*q[0]*(q[3] + kg[3]) - q[0] * q[0]*(2*q[3] + kg[3]) + kg[3]*( q[1] * q[1]
                  + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[1][3] = ( C3v*( 2*C3vNC*mn*( 16*mDelta2*mn*(-p0 + q[0])*q[1]*kg[3] + 4*pow(mDelta,3)*( q[0]*q[1]*(q[3] - 2*kg[3])
                  - 4*p0*q[1]*kg[3] + q[0]*kg[1]*kg[3] ) - 4*mn*(p0 - q[0])*( -( q[1]*q[3]*(2*kg[1] * kg[1]
                  + kg[2] * kg[2] ) ) + q[1]*( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] )*kg[3] + q[0] * q[0]*( q[3]*(q[1]
                  + kg[1]) - (2*q[1] + kg[1])*kg[3] ) + 2*p0*q[0]*( kg[1]*kg[3] + q[1]*(q[3] + 2*kg[3]) ) )
                + mDelta*( 4*p0*q[3]*( q[0] * q[0]*(q[1] - kg[1]) + 2*q[1]*kg[1] * kg[1] ) - 4*p0*( q[0] * q[0]*(2*q[1]
                  - 3*kg[1]) + q[1]*( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] ) )*kg[3] - 4*pow(p0,2)*q[0]*( 3*q[1]*q[3]
                  + 2*q[1]*kg[3] + 3*kg[1]*kg[3] ) + q[0]*( q[1]*(q[3] - 2*kg[3])
                  + kg[1]*kg[3] )*( -3*q[0] * q[0] + q[1] * q[1] + 2*q[1]*kg[1] + q[3]*(q[3] + 2*kg[3]) ) ) )
              - std::complex<double>(0.,1.)*C5aNC*kg[2]*( 4*pow(mDelta,3)*mn*( 8*mn2 + 2*(p0 - q[0])*q[0] + q[1]*(q[1] + kg[1]) )
                  - 2*mn2*q[3]*(q[3] + kg[3])*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*( 4*pow(p0,2)*q[0] * q[0] - 6*p0*pow(q[0],3) + 2*pow(q[0],4) + 4*p0*q[0]*q[1] * q[1]
                  - 3*q[0] * q[0]*q[1] * q[1] + pow(q[1],4) + q[1] * q[1]*q[3] * q[3] + 6*p0*q[0]*q[1]*kg[1]
                  - 5*q[0] * q[0]*q[1]*kg[1] + 3*pow(q[1],3)*kg[1] + q[1]*q[3] * q[3]*kg[1] + 2*q[1] * q[1]*kg[1] * kg[1]
                  + 2*q[3]*( (p0 - q[0])*q[0] + q[1]*(q[1] + kg[1]) )*kg[3] + 8*mn2*( 2*p0*q[0] - q[0] * q[0]
                  + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
                + mDelta2*( 8*pow(p0,3)*q[0] + 4*pow(p0,2)*( -3*q[0] * q[0] + q[1]*(q[1] + kg[1]) )
                  + 8*mn2*( 4*pow(p0,2) - 2*p0*q[0] - q[3]*kg[3] ) + 2*p0*q[0]*( 2*q[0] * q[0] - q[1]*(q[1] + kg[1])
                  - q[3]*(q[3] + 2*kg[3]) ) + q[3]*( q[0] * q[0]*(-q[3] + kg[3]) + q[1]*(q[3]*kg[1] - q[1]*kg[3]) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][0] = ( C3v*( 2*C3vNC*mn*kg[2]*( 4*pow(mDelta,3)*(q[1]*kg[1] + q[3]*kg[3]) - 4*mn*(p0 - q[0])*( q[0]*( q[1] * q[1]
                  + q[3] * q[3] - q[1]*kg[1] - q[3]*kg[3] ) + 2*p0*(q[1]*kg[1] + q[3]*kg[3]) )
                + mDelta*( -4*p0*q[0]*( q[1] * q[1] + q[3] * q[3] - 3*q[1]*kg[1] - 3*q[3]*kg[3] ) - 12*pow(p0,2)*(q[1]*kg[1]
                  + q[3]*kg[3]) + (q[1]*kg[1] + q[3]*kg[3])*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
              - std::complex<double>(0.,1.)*C5aNC*(q[3]*kg[1] - q[1]*kg[3])*( 8*pow(mDelta,3)*mn*(-p0 + q[0])
                + 2*mDelta*mn*(p0 - q[0])*( 8*mn2 - 2*p0*q[0] + q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] )
                  - 4*mn2*(p0 - q[0])*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta2*( -8*pow(p0,3) + 12*pow(p0,2)*q[0] + 8*mn2*(2*p0 + q[0]) + q[0]*( q[1] * q[1] + q[3] * q[3]
                  - q[1]*kg[1] - q[3]*kg[3] ) + 2*p0*( -2*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][1] = ( C3v*( 2*C3vNC*mn*kg[2]*( 4*pow(mDelta,3)*q[0]*kg[1] - 4*mn*(p0 - q[0])*( q[0] * q[0]*(q[1]
                  - kg[1]) + 2*p0*q[0]*kg[1] + q[3]*(q[3]*kg[1] - q[1]*kg[3]) )
                + mDelta*( -12*pow(p0,2)*q[0]*kg[1] + q[0]*kg[1]*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3]
                  + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) - 4*p0*( q[0] * q[0]*(q[1] - 3*kg[1]) + 2*q[3]*(q[3]*kg[1] - q[1]*kg[3]) ) ) )
              - std::complex<double>(0.,1.)*C5aNC*( -2*mn2*(q[1] + kg[1])*(-(q[3]*kg[1]) + q[1]*kg[3])*( 2*p0*q[0]
                  - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] ) + 4*pow(mDelta,3)*mn*( q[3]*( q[0] * q[0] + q[1]*kg[1]
                  - kg[2] * kg[2] ) + ( 8*mn2 + 2*(p0 - q[0])*q[0] + q[3] * q[3] )*kg[3])
                + mDelta2*( 8*pow(p0,3)*q[0]*kg[3] + 4*pow(p0,2)*( q[3]*kg[1]*(q[1] + kg[1]) + q[0] * q[0]*(q[3]
                  - 3*kg[3]) + q[3]*kg[3]*(q[3] + kg[3]) ) + 2*p0*q[0]*( q[3]*( q[1] * q[1] + q[3] * q[3] + q[1]*kg[1]
                  + 2*kg[1] * kg[1] + kg[2] * kg[2] ) - q[1]*(q[1] + 2*kg[1])*kg[3] + q[0] * q[0]*(-3*q[3] + 2*kg[3]) )
                  + (q[3]*kg[1] - q[1]*kg[3])*( q[0] * q[0]*(q[1] - kg[1]) + q[3]*(q[3]*kg[1] - q[1]*kg[3]) )
                  + 8*mn2*( 2*p0*q[0]*(q[3] - kg[3]) + 4*pow(p0,2)*kg[3] + kg[1]*(q[3]*kg[1] - q[1]*kg[3]) ) )
                + mDelta*mn*( q[3]*( -pow(q[0],4) + kg[1]*( pow(q[1],3) + q[1]*q[3] * q[3] + 2*q[1] * q[1]*kg[1]
                  - 2*q[3] * q[3]*kg[1] ) - ( q[1] * q[1] + 3*q[3] * q[3] + 2*q[1]*kg[1] )*kg[2] * kg[2]
                  + q[0] * q[0]*( -q[1] * q[1] + q[3] * q[3] - 3*q[1]*kg[1] + 2*kg[1] * kg[1] + 5*kg[2] * kg[2] ) )
                  + ( 2*pow(q[0],4) - q[0] * q[0]*( 3*q[3] * q[3] + 2*q[1]*kg[1] ) + q[3] * q[3]*( q[1] * q[1]
                  + q[3] * q[3] + 4*q[1]*kg[1] - 2*kg[2] * kg[2] ) )*kg[3] + 4*pow(p0,2)*q[0] * q[0]*(q[3] + kg[3])
                  + 8*mn2*( 2*p0*q[0]*(q[3] + kg[3]) - q[0] * q[0]*(2*q[3] + kg[3]) + kg[3]*( q[1] * q[1] + q[3] * q[3]
                  + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) + 2*p0*q[0]*( q[1] * q[1]*q[3] + q[1]*kg[1]*(3*q[3] + kg[3])
                  - q[0] * q[0]*(q[3] + 3*kg[3]) + q[3]*( q[3] * q[3] - kg[1] * kg[1] - 3*kg[2] * kg[2] + 3*q[3]*kg[3] ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][2] = -( C3v*( std::complex<double>(0.,1.)*C5aNC*kg[2]*(q[3]*kg[1] - q[1]*kg[3])*( 4*pow(mDelta,3)*mn
                  + mDelta2*( 8*mn2 + 4*pow(p0,2) + 2*p0*q[0] - q[0] * q[0] + q[1] * q[1] + q[3] * q[3] )
                  + 2*mn2*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*( 4*p0*q[0] - 3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
              + 2*C3vNC*mn*( 16*mDelta2*mn*(p0 - q[0])*( q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] )
                  + 4*pow(mDelta,3)*( q[0]*( q[0] * q[0] + q[1] * q[1] + q[3] * q[3] - 2*q[1]*kg[1] - kg[2] * kg[2]
                  - 2*q[3]*kg[3] ) + 4*p0*( q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] ) ) - 4*mn*(p0 - q[0])*( 4*pow(p0,2)*q[0] * q[0]
                  + 2*pow(q[0],4) + pow(q[1],3)*kg[1] + q[1]*q[3] * q[3]*kg[1] + q[1] * q[1]*kg[1] * kg[1]
                  - q[3] * q[3]*kg[1] * kg[1] - q[1] * q[1]*kg[2] * kg[2] - 2*q[3] * q[3]*kg[2] * kg[2]
                  + q[3]*( q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] )*kg[3] + 2*p0*q[0]*( -3*q[0] * q[0] + q[1] * q[1]
                  + q[3] * q[3] + 2*q[1]*kg[1] - kg[2] * kg[2] + 2*q[3]*kg[3] ) - q[0] * q[0]*( q[1] * q[1]
                  + 3*q[1]*kg[1] - kg[2] * kg[2] + 3*q[3]*kg[3] ) )
                - mDelta*( 16*pow(p0,3)*q[0] * q[0] + 4*pow(p0,2)*q[0]*( -5*q[0] * q[0] + 3*q[1] * q[1] + 3*q[3] * q[3]
                  + 2*q[1]*kg[1] - 3*kg[2] * kg[2] + 2*q[3]*kg[3] ) + q[0]*( q[0] * q[0] + q[1] * q[1] + q[3] * q[3]
                  - 2*q[1]*kg[1] - kg[2] * kg[2] - 2*q[3]*kg[3] )*( 3*q[0] * q[0] - q[1]*(q[1] + 2*kg[1]) - q[3]*(q[3] + 2*kg[3]) )
                  + 4*p0*( pow(q[1],3)*kg[1] + q[1]*q[3] * q[3]*kg[1] + q[3] * q[3]*( -2*kg[2] * kg[2] + q[3]*kg[3] )
                  + q[0] * q[0]*( q[1]*(-4*q[1] + kg[1]) + 3*kg[2] * kg[2] + q[3]*(-2*q[3] + kg[3]) )
                  + q[1] * q[1]*( 2*kg[1] * kg[1] + kg[3]*(q[3] + 2*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[2][3] = ( C3v*( 2*C3vNC*mn*kg[2]*( 4*pow(mDelta,3)*q[0]*kg[3] - 4*mn*(p0 - q[0])*( q[0] * q[0]*(q[3] - kg[3])
                  + 2*p0*q[0]*kg[3] + q[1]*(-(q[3]*kg[1]) + q[1]*kg[3]) )
                + mDelta*( -12*pow(p0,2)*q[0]*kg[3] + q[0]*kg[3]*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3]
                  + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) - 4*p0*( q[0] * q[0]*(q[3] - 3*kg[3]) + 2*q[1]*( -(q[3]*kg[1]) + q[1]*kg[3] ) ) ) )
              + std::complex<double>(0.,1.)*C5aNC*( -2*mn2*(q[3] + kg[3])*(q[3]*kg[1] - q[1]*kg[3])*( 2*p0*q[0]
                  - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] ) + 4*pow(mDelta,3)*mn*( q[0] * q[0]*(q[1] - 2*kg[1])
                  + 8*mn2*kg[1] + 2*p0*q[0]*kg[1] + q[1]*( q[1]*kg[1] - kg[2] * kg[2] + q[3]*kg[3] ) )
                + mDelta2*( 8*pow(p0,3)*q[0]*kg[1] + 4*pow(p0,2)*( q[0] * q[0]*(q[1] - 3*kg[1]) + q[1]*kg[1]*(q[1] + kg[1])
                  + q[1]*kg[3]*(q[3] + kg[3]) ) + (q[3]*kg[1] - q[1]*kg[3])*( q[0] * q[0]*(-q[3] + kg[3])
                  + q[1]*(q[3]*kg[1] - q[1]*kg[3]) ) + 8*mn2*( 2*p0*q[0]*(q[1] - kg[1]) + 4*pow(p0,2)*kg[1]
                  + kg[3]*( -(q[3]*kg[1]) + q[1]*kg[3]) ) - 2*p0*q[0]*( -pow(q[1],3) + q[0] * q[0]*(q[1] - 2*kg[1])
                  + q[3]*kg[1]*(q[3] + 2*kg[3]) + q[1]*( 2*kg[1] * kg[1] + kg[2] * kg[2] - q[3]*(q[3] + kg[3]) ) ) )
                + mDelta*mn*( pow(q[0],4)*q[1] - q[0] * q[0]*pow(q[1],3) + q[0] * q[0]*q[1]*q[3] * q[3] + 2*pow(q[0],4)*kg[1]
                  - 3*q[0] * q[0]*q[1] * q[1]*kg[1] + pow(q[1],4)*kg[1] + q[1] * q[1]*q[3] * q[3]*kg[1] - 2*q[0] * q[0]*q[1]*kg[1] * kg[1]
                  + 2*pow(q[1],3)*kg[1] * kg[1] - 2*q[1]*q[3] * q[3]*kg[1] * kg[1] + 4*pow(p0,2)*q[0] * q[0]*(q[1] + kg[1])
                  + 3*q[0] * q[0]*q[1]*kg[2] * kg[2] - pow(q[1],3)*kg[2] * kg[2] - 3*q[1]*q[3] * q[3]*kg[2] * kg[2]
                  - 2*q[1] * q[1]*kg[1]*kg[2] * kg[2] - 3*q[0] * q[0]*q[1]*q[3]*kg[3] + pow(q[1],3)*q[3]*kg[3] + q[1]*pow(q[3],3)*kg[3]
                  - 2*q[0] * q[0]*q[3]*kg[1]*kg[3] + 4*q[1] * q[1]*q[3]*kg[1]*kg[3] - 2*q[1]*q[3]*kg[2] * kg[2]*kg[3]
                  + 8*mn2*( 2*p0*q[0]*(q[1] + kg[1]) - q[0] * q[0]*(2*q[1] + kg[1]) + kg[1]*( q[1] * q[1] + q[3] * q[3]
                  + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) + 2*p0*q[0]*( pow(q[1],3) + 3*q[1] * q[1]*kg[1]
                  - 2*q[0] * q[0]*(q[1] + 2*kg[1]) + q[1]*( q[3] * q[3] + kg[1] * kg[1] - 2*kg[2] * kg[2]
                  + 3*q[3]*kg[3] ) + kg[1]*( kg[1] * kg[1] + kg[2] * kg[2] + kg[3]*(q[3] + kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][0] = -( C3v*( std::complex<double>(0.,1.)*C5aNC*q[1]*kg[2]*( 8*pow(mDelta,3)*mn*(-p0 + q[0])
                  + 2*mDelta*mn*(p0 - q[0])*( 8*mn2 - 2*p0*q[0] + q[0] * q[0] - q[1]*kg[1] - q[3]*kg[3] )
                  - 4*mn2*(p0 - q[0])*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta2*( -8*pow(p0,3) + 12*pow(p0,2)*q[0] + 8*mn2*(2*p0 + q[0]) + q[0]*( q[1] * q[1] + q[3] * q[3]
                  - q[1]*kg[1] - q[3]*kg[3] ) + 2*p0*( -2*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) )
              + 2*C3vNC*mn*( 16*mDelta2*mn*(p0 - q[0])*q[0]*q[3] + 4*pow(mDelta,3)*( q[3]*( 4*p0*q[0] + kg[1] * kg[1]
                  + kg[2] * kg[2] ) - q[1]*kg[1]*kg[3] ) + 4*mn*(p0 - q[0])*( -4*pow(p0,2)*q[0]*q[3] + 8*p0*q[0] * q[0]*q[3]
                  - 3*pow(q[0],3)*q[3] - 2*p0*q[3]*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[0]*q[3]*( q[1] * q[1]
                  + q[3] * q[3] + q[1]*kg[1] + kg[1] * kg[1] + kg[2] * kg[2] ) + 2*p0*q[1]*kg[1]*kg[3] + q[0]*(q[1] * q[1]
                  + 2*q[3] * q[3] - q[1]*kg[1])*kg[3] )
                + mDelta*( -16*pow(p0,3)*q[0]*q[3] + 4*pow(p0,2)*( 8*q[0] * q[0]*q[3] - 3*q[3]*(kg[1] * kg[1]
                  + kg[2] * kg[2]) + 3*q[1]*kg[1]*kg[3] ) + ( q[3]*( kg[1] * kg[1] + kg[2] * kg[2] )
                  - q[1]*kg[1]*kg[3] )*( -3*q[0] * q[0] + q[1] * q[1] + 2*q[1]*kg[1] + q[3]*(q[3] + 2*kg[3]) )
                  + 4*p0*q[0]*( -3*q[0] * q[0]*q[3] + q[1]*kg[1]*(q[3] - 3*kg[3]) + q[1] * q[1]*(q[3] + kg[3])
                  + q[3]*( 3*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*(q[3] + 2*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][1] = ( C3v*( std::complex<double>(0.,1.)*C5aNC*kg[2]*( -2*mn2*q[1]*(q[1] + kg[1])*( 2*p0*q[0] - q[0] * q[0]
                  + q[1]*kg[1] + q[3]*kg[3] ) + 4*pow(mDelta,3)*mn*( 8*mn2 + 2*(p0 - q[0])*q[0] + q[3]*(q[3] + kg[3]) )
                + mDelta*mn*( 4*pow(p0,2)*q[0] * q[0] - 6*p0*pow(q[0],3) + 2*pow(q[0],4) + 4*p0*q[0]*q[3] * q[3]
                  - q[0] * q[0]*q[3] * q[3] + q[1] * q[1]*q[3] * q[3] + pow(q[3],4) + 2*p0*q[0]*q[1]*kg[1] - 2*q[0] * q[0]*q[1]*kg[1]
                  + 2*q[1]*q[3] * q[3]*kg[1] - 2*q[3] * q[3]*kg[1] * kg[1] - 2*q[3] * q[3]*kg[2] * kg[2]
                  + q[3]*( 6*p0*q[0] - 5*q[0] * q[0] + q[1] * q[1] + 3*q[3] * q[3] + 2*q[1]*kg[1] )*kg[3]
                  + 8*mn2*( 2*p0*q[0] - q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) )
                + mDelta2*( 8*pow(p0,3)*q[0] + 8*mn2*( 4*pow(p0,2) - 2*p0*q[0] - q[1]*kg[1] )
                  + 2*p0*q[0]*( 2*q[0] * q[0] - q[1]*(q[1] + 2*kg[1]) - q[3]*(q[3] + kg[3]) ) + 4*pow(p0,2)*( -3*q[0] * q[0]
                  + q[3]*(q[3] + kg[3]) ) + q[1]*( q[0] * q[0]*(-q[1] + kg[1]) + q[3]*( -(q[3]*kg[1]) + q[1]*kg[3] ) ) ) )
              + 2*C3vNC*mn*( 16*mDelta2*mn*(-p0 + q[0])*q[3]*kg[1] + 4*pow(mDelta,3)*( q[0]*q[1]*q[3] - 4*p0*q[3]*kg[1]
                  + q[0]*kg[1]*(-2*q[3] + kg[3]) ) - 4*mn*(p0 - q[0])*( 2*p0*q[0]*( q[1]*q[3] + kg[1]*(2*q[3] + kg[3]) )
                  - q[0] * q[0]*( q[1]*(q[3] - kg[3]) + kg[1]*(2*q[3] + kg[3]) ) + q[3]*( q[1] * q[1]*kg[1]
                  + q[1]*( 2*kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) )
                + mDelta*( -4*pow(p0,2)*q[0]*( 3*q[1]*q[3] + 2*q[3]*kg[1] + 3*kg[1]*kg[3] )
                  + q[0]*( -3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3])*( q[1]*q[3] + kg[1]*(-2*q[3] + kg[3]) )
                  - 4*p0*( q[0] * q[0]*( -3*q[1]*q[3] + 2*q[3]*kg[1] + q[1]*kg[3] - 3*kg[1]*kg[3] )
                  + q[3]*( q[1] * q[1]*kg[1] + 2*q[1]*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][2] = ( C3v*( -2*C3vNC*mn*kg[2]*( 16*mDelta2*mn*(p0 - q[0])*q[3] + 4*pow(mDelta,3)*( 4*p0*q[3] + 2*q[0]*q[3] - q[0]*kg[3] )
                  + 4*mn*(p0 - q[0])*( q[3]*( 4*p0*q[0] - 2*q[0] * q[0] + q[3] * q[3] + q[1]*(q[1] + kg[1]) )
                  + ( 2*p0*q[0] - q[0] * q[0] + q[1] * q[1] + 2*q[3] * q[3])*kg[3] )
                + mDelta*( 4*pow(p0,2)*q[0]*(2*q[3] + 3*kg[3]) + q[0]*(2*q[3] - kg[3])*( -3*q[0] * q[0] + q[1] * q[1]
                  + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) + 4*p0*( q[0] * q[0]*(2*q[3] - 3*kg[3]) + ( q[1] * q[1]
                  + q[3] * q[3] )*(q[3] + 2*kg[3]) ) ) )
              - std::complex<double>(0.,1.)*C5aNC*( 2*mn2*q[1]*kg[2] * kg[2]*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                  + 4*pow(mDelta,3)*mn*( 8*mn2*kg[1] + 2*p0*q[0]*kg[1] - 2*q[0] * q[0]*kg[1] + q[1] * q[1]*kg[1] + q[3] * q[3]*kg[1]
                  + q[1]*kg[1] * kg[1] + q[1]*kg[2] * kg[2] + q[3]*kg[1]*kg[3])
                + mDelta2*( 8*pow(p0,3)*q[0]*kg[1] + q[1]*( -q[0] * q[0] + q[1] * q[1] + q[3] * q[3] )*kg[2] * kg[2]
                  + 8*mn2*( 2*p0*q[0]*(q[1] - kg[1]) + 4*pow(p0,2)*kg[1] + q[1]*kg[2] * kg[2] ) + 4*pow(p0,2)*( q[0] * q[0]*(q[1]
                  - 3*kg[1]) + q[1] * q[1]*kg[1] + q[1]*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + kg[3]) )
                  - 2*p0*q[0]*( -pow(q[1],3) + 2*q[0] * q[0]*(q[1] - kg[1]) + q[3]*kg[1]*(q[3] + kg[3])
                  - q[1]*( -kg[1] * kg[1] + kg[2] * kg[2] + q[3]*(q[3] + kg[3]) ) ) )
                + mDelta*mn*( 4*pow(q[0],4)*q[1] - 2*q[0] * q[0]*pow(q[1],3) - 2*q[0] * q[0]*q[1]*q[3] * q[3] + 2*pow(q[0],4)*kg[1]
                  - 5*q[0] * q[0]*q[1] * q[1]*kg[1] + pow(q[1],4)*kg[1] - q[0] * q[0]*q[3] * q[3]*kg[1] + 2*q[1] * q[1]*q[3] * q[3]*kg[1]
                  + pow(q[3],4)*kg[1] - 5*q[0] * q[0]*q[1]*kg[1] * kg[1] + 3*pow(q[1],3)*kg[1] * kg[1] + 3*q[1]*q[3] * q[3]*kg[1] * kg[1]
                  + 2*q[1] * q[1]*pow(kg[1],3) - 2*q[3] * q[3]*pow(kg[1],3) + 4*pow(p0,2)*q[0] * q[0]*(q[1] + kg[1])
                  - 3*q[0] * q[0]*q[1]*kg[2] * kg[2] + pow(q[1],3)*kg[2] * kg[2] + q[1]*q[3] * q[3]*kg[2] * kg[2]
                  + 2*q[1] * q[1]*kg[1]*kg[2] * kg[2] - 2*q[3] * q[3]*kg[1]*kg[2] * kg[2] + q[3]*( -( q[0] * q[0]*(2*q[1]
                  + 5*kg[1]) ) + kg[1]*( 3*( q[1] * q[1] + q[3] * q[3] ) + 4*q[1]*kg[1] ) + 2*q[1]*kg[2] * kg[2] )*kg[3]
                  + 8*mn2*( 2*p0*q[0]*(q[1] + kg[1]) - q[0] * q[0]*(2*q[1] + kg[1]) + kg[1]*( q[1] * q[1] + q[3] * q[3]
                  + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) ) + 2*p0*q[0]*( pow(q[1],3) + 3*q[1] * q[1]*kg[1] - q[0] * q[0]*(4*q[1]
                  + 3*kg[1]) + q[3]*kg[1]*(2*q[3] + 3*kg[3]) + q[1]*( 3*kg[1] * kg[1] + 2*kg[2] * kg[2] + q[3]*(q[3] + kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  tr[3][3] = -( C3v*( std::complex<double>(0.,1.)*C5aNC*q[1]*kg[2]*( 4*pow(mDelta,3)*mn*(q[3] + kg[3])
                  + 2*mn2*(q[3] + kg[3])*( 2*p0*q[0] - q[0] * q[0] + q[1]*kg[1] + q[3]*kg[3] )
                + mDelta*mn*(q[3] + kg[3])*( 4*p0*q[0] - 3*q[0] * q[0] + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] )
                + mDelta2*( -(q[1]*q[3]*kg[1]) + q[0] * q[0]*(q[3] - kg[3]) + 8*mn2*kg[3] + 2*p0*q[0]*kg[3] + q[1] * q[1]*kg[3]
                  + 4*pow(p0,2)*(q[3] + kg[3]) ) )
              + 2*C3vNC*mn*( 16*mDelta2*mn*(p0 - q[0])*(q[0] * q[0] - q[1]*kg[1]) + 4*pow(mDelta,3)*( 4*p0*(q[0] * q[0]
                  - q[1]*kg[1]) + q[0]*(pow(q[1] - kg[1],2) + kg[2] * kg[2]) ) - 4*mn*(p0 - q[0])*( 4*pow(p0,2)*q[0] * q[0]
                  + 3*pow(q[0],4) + 2*p0*q[0]*( -4*q[0] * q[0] + pow(q[1] + kg[1],2) + kg[2] * kg[2])
                  - q[0] * q[0]*( q[3] * q[3] + (q[1] + kg[1])*(2*q[1] + kg[1]) + kg[2] * kg[2] + 2*q[3]*kg[3] )
                  + q[1]*( q[1] * q[1]*kg[1] + q[1]*( 2*kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) )
                + mDelta*( -16*pow(p0,3)*q[0] * q[0] + 4*pow(p0,2)*q[0]*( 8*q[0] * q[0] - 3*q[1] * q[1] - 2*q[1]*kg[1]
                  - 3*(kg[1] * kg[1] + kg[2] * kg[2])) + q[0]*( pow(q[1] - kg[1],2) + kg[2] * kg[2] )*( -3*q[0] * q[0]
                  + q[1] * q[1] + q[3] * q[3] + 2*q[1]*kg[1] + 2*q[3]*kg[3] ) - 4*p0*( 3*pow(q[0],4)
                  - q[0] * q[0]*( 4*q[1] * q[1] - q[1]*kg[1] + 3*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*(q[3] + 2*kg[3]) )
                  + q[1]*( q[1] * q[1]*kg[1] + 2*q[1]*( kg[1] * kg[1] + kg[2] * kg[2] ) + q[3]*kg[1]*(q[3] + 2*kg[3]) ) ) ) ) )
              ) / (24.*mDelta2*pow(mn,4));

  return tr;

}
