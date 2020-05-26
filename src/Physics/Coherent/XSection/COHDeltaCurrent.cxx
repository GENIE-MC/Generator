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
  TLorentzVecotr   out_neutrino = i -> Kine().FSLeptonP4() ; 

  TLorentzVector t_q = *probe - out_neutrino ; 
  double Q2 = t_q.M2() ;

  delete probe ; 

    // this requires Q2
  double C3V   = delta_ff -> C3V( Q2 ) ;
  double C3vNC = delta_ff -> C3VNC( Q2 ) ;

  double mDelta = util::res::Mass( this -> Resonance() ) ; 
  double mDelta2 = pow( mDelta, 2 ); 
  
  // the following contractions requires a vector with time coordinate in the 0-th position
  // that is not the case for TLorentzVector so we need to rearrange it
  std::array<double, 4> q = { t_q.E(), t_q.X(), t_q.Y(), t_q.Z() } ;



  GTrace tr ;

  tr[0][0] = -(C3v*C3vNC*( mDelta*p0 + 
			   this->mn*(p0 + this->q[0]))*(4*mDelta2*this->q[1]*this->kg[1] - 4*pow(p0,2)*this->q[1]*this->kg[1] + 
							pow(this->q[1],3)*this->kg[1] + 
							this->q[1]*this->q[3] * this->q[3]*this->kg[1] +
							this->q[1] * this->q[1]*this->kg[1] * this->kg[1] - 
							this->q[3] * this->q[3]*this->kg[1] * this->kg[1] - 
							this->q[3] * this->q[3]*this->kg[2] * this->kg[2] +
							this->q[0] * this->q[0]*(this->q[1] * this->q[1] - 
										 3*this->q[1]*this->kg[1] + this->q[3]*(2*this->q[3] - 3*this->kg[3])) +
							this->q[3]*(4*mDelta2 - 4*pow(p0,2) + this->q[1] * this->q[1] + 
								    this->q[3] * this->q[3] + 2*this->q[1]*this->kg[1])*this->kg[3] - 
							8*p0*this->q[0]*(this->q[1]*this->kg[1] + this->q[3]*this->kg[3])))/
    (3.*mDelta2*pow(this->mn,3));


}
//____________________________________________________________________________
GTrace DirTrace( const Interaction * i,
		 const COHFormFactorI * ff ) const {


}
