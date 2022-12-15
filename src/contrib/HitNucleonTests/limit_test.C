#include <iostream>

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "Framework/Utils/KineUtils.h"

#include "Framework/Utils/Range1.h"

#include "Math/Vector4D.h"
#include "Math/Boost.h" 


ostream & operator << ( ostream & out, const genie::Range1D_t & r ) {
  
  out << '[' << r.min << ", " << r.max << ']' ;
  return out ;
}

template< class A, class B>
void print_initial_conditions( A & probe, 
			       B & target ) {

  cout << "Probe: (" << probe.x() << ", " << probe.y() << ", " << probe.z() << ", " << probe.e() << " )" << endl ;
  
  cout << "Target: (" << target.x() << ", " << target.y() << ", " << target.z() << ", " << target.e() << " ) with beta: " << target.Beta() << endl ;
  
  auto total = probe + target ;
  auto sqrt_s = total.M() ;

  cout << "System CoM energy: " << sqrt_s << endl ;
}

bool compare_ranges( const genie::Range1D_t & old, const genie::Range1D_t & new_r ) {

  cout << "Old: " << old << endl;
  cout << "New: " << new_r << endl;
  
  if ( old.min == new_r.min && old.max == new_r.max ) {
    cout << "Same range" << endl;
    return true ;
  }
  else {
    auto min_dev = (new_r.min - old.min)/old.min;
    auto max_dev = (new_r.max - old.max)/old.max;
    cout << "Different resutls: " << endl 
	 << "Min deviation "  << min_dev*100 
	 << " % , Max deviation " << max_dev*100 << " %" << endl ;
    
    return false;
  }

}



void Wlimit( ) {

  using namespace genie;
  using namespace std;

  cout << "--------- Neutrino on Target on shell ------------" << endl;
  
  ROOT::Math::PxPyPzEVector probe( 0., 0., 2., 2. );

  double proton_mass = PDGLibrary::Instance()->Find(kPdgProton)->Mass();
  ROOT::Math::PxPyPzMVector nucleon( 0.1, 0.2, 0.1, proton_mass );

  print_initial_conditions( probe, nucleon );
  
  auto total = probe + nucleon ;
  auto sqrt_s = total.M() ;

  double mu_mass =  PDGLibrary::Instance()->Find(kPdgMuon)->Mass();
  auto new_range = utils::kinematics::InelWLim(sqrt_s, mu_mass);

  // old range
  auto boost = nucleon.BoostToCM();
  ROOT::Math::Boost t(boost);
  auto old_probe = t(probe);
  auto old_range =  utils::kinematics::InelWLim(old_probe.E(), proton_mass, mu_mass);  

  compare_ranges( old_range, new_range );


  cout << endl << "--------- Neutrino on Target off shell, small beta --------" << endl;

  nucleon.SetPxPyPzE( 0.1, 0.2, 0.1, proton_mass - 0.1 );

  print_initial_conditions( probe, nucleon);

  total = probe + nucleon ;
  sqrt_s = total.M() ;
  new_range = utils::kinematics::InelWLim(sqrt_s, mu_mass);

  boost = nucleon.BoostToCM();
  t.SetComponents(boost);
  old_probe = t(probe);

  old_range =  utils::kinematics::InelWLim(old_probe.E(), nucleon.M(), mu_mass);
  compare_ranges( old_range, new_range );

  // this si the old range calculated using the on shell mass, obvioulsy wrong
  old_range =  utils::kinematics::InelWLim(old_probe.E(), proton_mass, mu_mass);
  cout << "Range with wrong inputs: " << old_range << " for comparison only" << endl;


  cout << endl << "--------- Neutrino on Target off shell, beta > 1 --------" << endl;

  nucleon.SetPxPyPzE( 0.2, 0.8, 0.1, proton_mass - 0.2 );

  print_initial_conditions( probe, nucleon);

  total = probe + nucleon ;
  sqrt_s = total.M() ;

  new_range = utils::kinematics::InelWLim(sqrt_s, mu_mass);

  std::cout << "Range: " << new_range << std::endl;

}


void Q2_W_limit( ) {

  using namespace genie;
  using namespace std;

  cout << "--------- Neutrino on Target on shell ------------" << endl;

  ROOT::Math::PxPyPzEVector probe( 0., 0., 2., 2. );

  double proton_mass = PDGLibrary::Instance()->Find(kPdgProton)->Mass();

  ROOT::Math::PxPyPzMVector nucleon( 0.1, 0.2, 0.1, proton_mass );

  print_initial_conditions( probe, nucleon);
 
  auto total = probe + nucleon ;
  auto sqrt_s = total.M() ;
  auto com_boost = total.BoostToCM();

  ROOT::Math::Boost com_t(com_boost);

  auto com_probe = com_t(probe);

  double mu_mass =  PDGLibrary::Instance()->Find(kPdgMuon)->Mass();

  double W = 1.3 ;

  auto new_range = utils::kinematics::NewInelQ2Lim_W(sqrt_s, com_probe.E(), mu_mass, W );

    // old range

  auto boost = nucleon.BoostToCM();
  
  ROOT::Math::Boost t(boost);
  
  auto old_probe = t(probe);

  auto old_range =  utils::kinematics::InelQ2Lim_W(old_probe.E(), proton_mass, mu_mass, W);  

  compare_ranges( old_range, new_range);


  cout << endl << "--------- Neutrino on Target off shell, small beta --------" << endl;

  nucleon.SetPxPyPzE( 0.1, 0.2, 0.1, proton_mass - 0.1 );

  print_initial_conditions( probe, nucleon);

  total = probe + nucleon ;
  sqrt_s = total.M() ;
  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);

  com_probe = com_t(probe);
  

  new_range = utils::kinematics::NewInelQ2Lim_W(sqrt_s, com_probe.E(), mu_mass, W );

  boost = nucleon.BoostToCM();
  t.SetComponents(boost);
  old_probe = t(probe);



  old_range =  utils::kinematics::InelQ2Lim_W(old_probe.E(), nucleon.M(), mu_mass, W);

  compare_ranges( old_range, new_range );

  auto wrong_range =  utils::kinematics::InelQ2Lim_W(old_probe.E(), proton_mass, mu_mass, W);
  cout << "Range with wrong inputs: " << wrong_range << " for comparison only" << endl;


  cout << endl << "--------- Neutrino on Target off shell, beta > 1 --------" << endl;
  nucleon.SetPxPyPzE( 0.2, 0.8, 0.1, proton_mass - 0.2 );

  print_initial_conditions( probe, nucleon);
  
  total = probe + nucleon ;
  sqrt_s = total.M() ;
  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);
  com_probe = com_t(probe);
  
  
  new_range = utils::kinematics::NewInelQ2Lim_W(sqrt_s, com_probe.E(), mu_mass, W );


  std::cout << "Range: " << new_range << std::endl;


}


void Q2_limit( ) {

  using namespace genie;
  using namespace std;

  cout << "--------- Neutrino on Target on shell ------------" << endl;
  ROOT::Math::PxPyPzEVector probe( 0., 0., 2., 2. );

  double proton_mass = PDGLibrary::Instance()->Find(kPdgProton)->Mass();

  ROOT::Math::PxPyPzMVector nucleon( 0.1, 0.2, 0.1, proton_mass );
  
  print_initial_conditions( probe, nucleon );
  
  auto total = probe + nucleon ;
  auto sqrt_s = total.M() ;
  auto com_boost = total.BoostToCM();

  ROOT::Math::Boost com_t(com_boost);

  auto com_probe = com_t(probe);

  double mu_mass =  PDGLibrary::Instance()->Find(kPdgMuon)->Mass();

  auto new_range = utils::kinematics::NewInelQ2Lim(sqrt_s, com_probe.E(), mu_mass );

  // old range

  auto boost = nucleon.BoostToCM();
  
  ROOT::Math::Boost t(boost);
  
  auto old_probe = t(probe);

  auto old_range =  utils::kinematics::InelQ2Lim(old_probe.E(), proton_mass, mu_mass );  

  compare_ranges( old_range, new_range );

  

  cout << endl << "--------- Neutrino on Target off shell, small beta --------" << endl;

  nucleon.SetPxPyPzE( 0.1, 0.2, 0.1, proton_mass - 0.1 );
  
  print_initial_conditions(probe, nucleon);

  total = probe + nucleon ;
  sqrt_s = total.M() ;
  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);

  com_probe = com_t(probe);
 
  new_range = utils::kinematics::NewInelQ2Lim(sqrt_s, com_probe.E(), mu_mass );

  boost = nucleon.BoostToCM();
  t.SetComponents(boost);
  old_probe = t(probe);

  old_range =  utils::kinematics::InelQ2Lim(old_probe.E(), nucleon.M(), mu_mass);

  compare_ranges( old_range, new_range );

  old_range =  utils::kinematics::InelQ2Lim(old_probe.E(), proton_mass, mu_mass);

  cout << "Range with wrong inputs: " << old_range << " for comparison only" << endl;


  cout << endl << "--------- Neutrino on Target off shell, beta > 1 --------" << endl;

  nucleon.SetPxPyPzE( 0.2, 0.8, 0.1, proton_mass - 0.2 );

  print_initial_conditions( probe, nucleon);

  total = probe + nucleon ;
 
  sqrt_s = total.M() ;

  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);

  com_probe = com_t(probe);
  
  new_range = utils::kinematics::NewInelQ2Lim(sqrt_s, com_probe.E(), mu_mass );

  std::cout << "Range: " << new_range << std::endl;


}
