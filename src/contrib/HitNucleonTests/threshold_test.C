#include <iostream>
#include <memory>

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

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

  double m_l = 
  double m_p = i.InitState().Probe()->Mass();
  double E_T = i.InitState().Tgt().HitNucP4().E();
  double p_T = i.InitState().Tgt().HitNucP4().p();
  double W = 
  

  double alpha = 

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
  

}
