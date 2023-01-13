#include <iostream>
#include <memory>
#include <cmath>

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

#include "Framework/Interaction/Interaction.h"

#include "Math/Vector4D.h"
#include "Math/Boost.h" 


using namespace genie;
using namespace std;


double find_threshold( const Interaction & i, double prec = 0.00001, double max = 10. ) {

  double min = 0. ;

  Interaction in(i);

  // max - min > prec * 0.5 * ( max + min)
  while ( (2-prec)*max > (2+prec)*min ) {

    double centre = 0.5*(max+min);
    in.InitStatePtr()->SetProbeE(centre);
    
    if ( in.PhaseSpace().IsAboveThreshold() ) {
      max = centre;
    }
    else { 
      min = centre; 
    }
    
  }
  
  return 0.5*(min+max);

}


double new_threshold( const Interaction & i ) {

  int finalNucPDG = i.InitState().Tgt().HitNucPdg();
  if ( i.ProcInfo().IsWeakCC() ) finalNucPDG = pdg::SwitchProtonNeutron( finalNucPDG );

  double m_l = i.FSPrimLepton()->Mass();
  double m_p = i.InitState().Probe()->Mass();
  double E_T = i.InitState().Tgt().HitNucP4().E();
  double p_T = i.InitState().Tgt().HitNucP4().P();
  double W = PDGLibrary::Instance()->Find( finalNucPDG )->Mass();
  
  double alpha = m_l*m_l + W*W + 2*W*m_l - m_p*m_p - E_T*E_T + p_T*p_T ;

  double p_T_p = i.InitState().Tgt().HitNucP4().Pz();

  if ( E_T > abs(p_T_p) ) {
    auto a2 = alpha*alpha;
    auto Et2 = E_T*E_T;
    auto ptp2 = p_T_p*p_T_p;
    auto mp2 = m_p*m_p;

    auto delta = a2*Et2 - 4*(Et2-ptp2)*(ptp2*mp2+a2/4);

    if( delta < 0. ) return m_p ;

    auto threshold = alpha*E_T;
    auto sqrt_delta = sqrt( delta );

    if ( p_T_p >= 0. ) {
      threshold += sqrt_delta;
    }
    else {
      threshold -= sqrt_delta;
    }
    threshold /= 2*(Et2 - ptp2);
    
    return threshold;

  }
  else {

    if ( p_T_p < 0 ) {
      if ( alpha < 2*m_p*E_T ) return m_p ;
      else {
	auto a2 = alpha*alpha;
	auto Et2 = E_T*E_T;
	auto ptp2 = p_T_p*p_T_p;
	auto mp2 = m_p*m_p;

	auto delta = a2*Et2 - 4*(Et2-ptp2)*(ptp2*mp2+a2/4);
	auto threshold = 0.5*(alpha*E_T + sqrt( delta ))/(Et2 - ptp2);
	return threshold ;
      }
    }
    else { 
      return m_p;
    }
  }
  
}






void threshold_test() {


  unique_ptr<Interaction> i_p {Interaction::QELCC( 1000180400, kPdgNeutron, kPdgNuMu, 2. )} ;

  // auto th = i_p -> PhaseSpace().Threshold();

  // cout << th << endl;

  TLorentzVector on_shell_target;
  on_shell_target.SetXYZM( 0.1, -0.2, 0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.2 );
    
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( on_shell_target ) ;

  i_p -> InitState().Tgt().HitNucP4().Print();

  cout << *i_p << endl;

  auto th = i_p -> PhaseSpace().Threshold();

  cout << th << endl;

  auto lab_th = find_threshold(*i_p);

  cout << lab_th << endl ;

  auto new_th = new_threshold(*i_p);
  cout << new_th << endl ;

}
