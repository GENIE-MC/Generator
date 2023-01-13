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


double find_threshold( const Interaction & i, double prec = 0.000001, double max = 10. ) {

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

    cout << "Simple case" << endl;
    
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

    cout << "Weird case, alpha = " << alpha << endl;

    if ( p_T_p < 0 ) {
      if ( alpha < 2*m_p*E_T ) {
	cout << "No intersection" << endl;
	return m_p ;
      }
      else {
	auto a2 = alpha*alpha;
	auto Et2 = E_T*E_T;
	auto ptp2 = p_T_p*p_T_p;
	auto mp2 = m_p*m_p;

	auto delta = a2*Et2 - 4*(Et2-ptp2)*(ptp2*mp2+a2/4);
	auto threshold = 0.5*(alpha*E_T - sqrt( delta ))/(Et2 - ptp2);
	return threshold ;
      }
    }
    else { 
      return m_p;
    }
  }
  
}


void print_initial_conditions( TLorentzVector tgt ) {

  cout << "Target: ( " << tgt.X() << ", " << tgt.Y() << ", " << tgt.Z() << ", " << tgt.E() << " ) with beta: " << tgt.Beta() << endl ;

}


void compare_threshold( double old_th, double new_th ) {

  cout << "Old threshold: " << old_th << " GeV" << endl;
  cout << "New threshold: " << new_th << " GeV" << endl;
  cout << "Relative difference:  " << (new_th - old_th)/old_th*100 << " %" << endl;


}



void threshold_test() {


  unique_ptr<Interaction> i_p {Interaction::QELCC( 1000180400, kPdgNeutron, kPdgNuMu, 2. )} ;

  cout << "--------- Neutrino on Target on shell, forward nucleon ------------" << endl;

  TLorentzVector on_shell_target;
  on_shell_target.SetXYZM( 0.1, -0.2, 0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( on_shell_target ) ;

  print_initial_conditions( on_shell_target );
  
  auto old_th = find_threshold(*i_p);
  auto new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );
  

  cout << "--------- Neutrino on Target on shell, backward nucleon ------------" << endl;

  on_shell_target.SetXYZM( 0.1, -0.2, -0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( on_shell_target ) ;

  print_initial_conditions( on_shell_target );
  
  old_th = find_threshold(*i_p);
  new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );

 
  cout << "--------- Neutrino on Target off shell, small beta, forward nucleon ------------" << endl;
  TLorentzVector off_shell_target;
  off_shell_target.SetXYZM( 0.1, -0.2, 0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.2 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );
  
  old_th = find_threshold(*i_p);
  new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );

  cout << "--------- Neutrino on Target off shell, small beta, backward nucleon ------------" << endl;
  off_shell_target.SetXYZM( 0.1, -0.2, -0.25, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.2 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );
  
  old_th = find_threshold(*i_p);
  new_th = new_threshold(*i_p); 

  compare_threshold( old_th, new_th );
  
  cout << endl << "--------- Neutrino on Target off shell, beta > 1, forward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.8, 0.1, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;

  cout << endl << "--------- Neutrino on Target off shell, beta > 1, backward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.8, -0.1, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;

  cout << endl << "--------- Neutrino on Target off shell, beta > 1, forward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.1, 0.8, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;

  cout << endl << "--------- Neutrino on Target off shell, beta > 1, backward nucleon --------" << endl;
  off_shell_target.SetXYZT( 0.2, 0.1, -0.8, PDGLibrary::Instance()->Find(kPdgNeutron)->Mass() - 0.25 );
  i_p -> InitStatePtr() -> TgtPtr() -> SetHitNucP4( off_shell_target ) ;

  print_initial_conditions( off_shell_target );

  new_th = new_threshold(*i_p); 

  cout << "New threshold: " << new_th << " GeV" << endl;



}
