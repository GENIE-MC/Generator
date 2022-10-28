#include <iostream>

#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "Framework/Utils/KineUtils.h"

#include "Framework/Utils/Range1.h"

#include "Math/Vector4D.h"
#include "Math/Boost.h" 


ostream & operator << ( ostream & out, const genie::Range1D_t r ) {
  
  out << '[' << r.min << ", " << r.max << ']' ;
  return out ;
}



void Wlimit( ) {

  using namespace genie;

  // neutrino case on shell
  ROOT::Math::PxPyPzEVector probe( 0., 0., 2., 2. );

  double proton_mass = PDGLibrary::Instance()->Find(kPdgProton)->Mass();

  ROOT::Math::PxPyPzMVector nucleon( 0.1, 0.2, 0.1, proton_mass );
  
  std::cout << probe.E() << std::endl ;
  std::cout << nucleon.E() << std::endl ;
  
  auto total = probe + nucleon ;
  auto sqrt_s = total.M() ;

  std::cout << sqrt_s << std::endl;

  double mu_mass =  PDGLibrary::Instance()->Find(kPdgMuon)->Mass();

  auto new_range = utils::kinematics::InelWLim(sqrt_s, mu_mass);

  std::cout << "New: " << new_range << std::endl;

  // old range

  auto boost = nucleon.BoostToCM();
  
  ROOT::Math::Boost t(boost);
  
  auto old_probe = t(probe);

  std::cout << old_probe << std::endl;


  auto old_range =  utils::kinematics::InelWLim(old_probe.E(), proton_mass, mu_mass);  

  std::cout << "Old: " << old_range << std::endl;
  

  // neutrino case off shell nucleon

  nucleon.SetPxPyPzE( 0.1, 0.2, 0.1, proton_mass - 0.1 );
  std::cout << "Off shell target: " << nucleon << std::endl ;

  total = probe + nucleon ;
  sqrt_s = total.M() ;

  std::cout << sqrt_s << std::endl;
  new_range = utils::kinematics::InelWLim(sqrt_s, mu_mass);

  std::cout << "New: " << new_range << std::endl;

  boost = nucleon.BoostToCM();
  t.SetComponents(boost);
  old_probe = t(probe);

  old_range =  utils::kinematics::InelWLim(old_probe.E(), proton_mass, mu_mass);

  std::cout << "Old: " << old_range << std::endl;

  old_range =  utils::kinematics::InelWLim(old_probe.E(), nucleon.M(), mu_mass);
  std::cout << "Old with correction: " << old_range << std::endl;


  // neutrino case very off shell (beta > 1)

  nucleon.SetPxPyPzE( 0.2, 0.8, 0.1, proton_mass - 0.2 );
  std::cout << "very Off shell target: " << nucleon << std::endl ;

  std::cout << "Beta: " << nucleon.Beta() << std::endl;

  total = probe + nucleon ;
 
  sqrt_s = total.M() ;

  std::cout << sqrt_s << std::endl;
  new_range = utils::kinematics::InelWLim(sqrt_s, mu_mass);

  std::cout << "New: " << new_range << std::endl;

}


void Q2_W_limit( ) {

  using namespace genie;

  // neutrino case on shell
  ROOT::Math::PxPyPzEVector probe( 0., 0., 2., 2. );

  double proton_mass = PDGLibrary::Instance()->Find(kPdgProton)->Mass();

  ROOT::Math::PxPyPzMVector nucleon( 0.1, 0.2, 0.1, proton_mass );
  
  std::cout << probe.E() << std::endl ;
  std::cout << nucleon.E() << std::endl ;
  
  auto total = probe + nucleon ;
  auto sqrt_s = total.M() ;
  auto com_boost = total.BoostToCM();

  ROOT::Math::Boost com_t(com_boost);

  auto com_probe = com_t(probe);

  std::cout << sqrt_s << std::endl;

  double mu_mass =  PDGLibrary::Instance()->Find(kPdgMuon)->Mass();

  double W = 1.3 ;

  auto new_range = utils::kinematics::NewInelQ2Lim_W(sqrt_s, com_probe.E(), mu_mass, W );

  std::cout << "New: " << new_range << std::endl;

  // old range

  auto boost = nucleon.BoostToCM();
  
  ROOT::Math::Boost t(boost);
  
  auto old_probe = t(probe);

  std::cout << old_probe << std::endl;

  auto old_range =  utils::kinematics::InelQ2Lim_W(old_probe.E(), proton_mass, mu_mass, W);  

  std::cout << "Old: " << old_range << std::endl;
  

  // neutrino case off shell nucleon

  nucleon.SetPxPyPzE( 0.1, 0.2, 0.1, proton_mass - 0.1 );
  std::cout << "Off shell target: " << nucleon << std::endl ;

  total = probe + nucleon ;
  sqrt_s = total.M() ;
  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);

  com_probe = com_t(probe);
  

  std::cout << sqrt_s << std::endl;
  new_range = utils::kinematics::NewInelQ2Lim_W(sqrt_s, com_probe.E(), mu_mass, W );

  std::cout << "New: " << new_range << std::endl;

  boost = nucleon.BoostToCM();
  t.SetComponents(boost);
  old_probe = t(probe);

  old_range =  utils::kinematics::InelQ2Lim_W(old_probe.E(), proton_mass, mu_mass, W);

  std::cout << "Old: " << old_range << std::endl;

  old_range =  utils::kinematics::InelQ2Lim_W(old_probe.E(), nucleon.M(), mu_mass, W);
  std::cout << "Old with correction: " << old_range << std::endl;

  // neutrino case very off shell (beta > 1)

  nucleon.SetPxPyPzE( 0.2, 0.8, 0.1, proton_mass - 0.2 );
  std::cout << "very Off shell target: " << nucleon << std::endl ;

  std::cout << "Beta: " << nucleon.Beta() << std::endl;

  total = probe + nucleon ;
 
  sqrt_s = total.M() ;

  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);

  com_probe = com_t(probe);
  

  std::cout << sqrt_s << std::endl;
  new_range = utils::kinematics::NewInelQ2Lim_W(sqrt_s, com_probe.E(), mu_mass, W );


  std::cout << "New: " << new_range << std::endl;


}


void Q2_limit( ) {

  using namespace genie;

  // neutrino case on shell
  ROOT::Math::PxPyPzEVector probe( 0., 0., 2., 2. );

  double proton_mass = PDGLibrary::Instance()->Find(kPdgProton)->Mass();

  ROOT::Math::PxPyPzMVector nucleon( 0.1, 0.2, 0.1, proton_mass );
  
  std::cout << probe.E() << std::endl ;
  std::cout << nucleon.E() << std::endl ;
  
  auto total = probe + nucleon ;
  auto sqrt_s = total.M() ;
  auto com_boost = total.BoostToCM();

  ROOT::Math::Boost com_t(com_boost);

  auto com_probe = com_t(probe);

  std::cout << sqrt_s << std::endl;

  double mu_mass =  PDGLibrary::Instance()->Find(kPdgMuon)->Mass();

  double W = 1.3 ;

  auto new_range = utils::kinematics::NewInelQ2Lim(sqrt_s, com_probe.E(), mu_mass );

  std::cout << "New: " << new_range << std::endl;

  // old range

  auto boost = nucleon.BoostToCM();
  
  ROOT::Math::Boost t(boost);
  
  auto old_probe = t(probe);

  std::cout << old_probe << std::endl;

  auto old_range =  utils::kinematics::InelQ2Lim(old_probe.E(), proton_mass, mu_mass );  

  std::cout << "Old: " << old_range << std::endl;
  

  // neutrino case off shell nucleon

  nucleon.SetPxPyPzE( 0.1, 0.2, 0.1, proton_mass - 0.1 );
  std::cout << "Off shell target: " << nucleon << std::endl ;

  total = probe + nucleon ;
  sqrt_s = total.M() ;
  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);

  com_probe = com_t(probe);
  

  std::cout << sqrt_s << std::endl;
  new_range = utils::kinematics::NewInelQ2Lim(sqrt_s, com_probe.E(), mu_mass );

  std::cout << "New: " << new_range << std::endl;

  boost = nucleon.BoostToCM();
  t.SetComponents(boost);
  old_probe = t(probe);

  old_range =  utils::kinematics::InelQ2Lim(old_probe.E(), proton_mass, mu_mass);

  std::cout << "Old: " << old_range << std::endl;

  old_range =  utils::kinematics::InelQ2Lim(old_probe.E(), nucleon.M(), mu_mass);
  std::cout << "Old with correction: " << old_range << std::endl;

  // neutrino case very off shell (beta > 1)

  nucleon.SetPxPyPzE( 0.2, 0.8, 0.1, proton_mass - 0.2 );
  std::cout << "very Off shell target: " << nucleon << std::endl ;

  std::cout << "Beta: " << nucleon.Beta() << std::endl;

  total = probe + nucleon ;
 
  sqrt_s = total.M() ;

  com_boost = total.BoostToCM();

  com_t.SetComponents(com_boost);

  com_probe = com_t(probe);
  

  std::cout << sqrt_s << std::endl;
  new_range = utils::kinematics::NewInelQ2Lim(sqrt_s, com_probe.E(), mu_mass );


  std::cout << "New: " << new_range << std::endl;


}
